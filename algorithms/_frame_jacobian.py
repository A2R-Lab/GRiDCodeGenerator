"""General-frame geometric Jacobian + operational-space inertia codegen (E2).

NEW, additive algorithm family. Nothing in the existing emit path is touched:
the family is only emitted when the ``frame_jacobian`` key is selected, mirroring
how the centroidal quick-wins (`_centroidal.py`) gate on their grid:: deps.

Lives in the kinematics (s_XmatsHom) domain: one inner pass builds the world
homogeneous transform of every joint (BFS chain-up, identical to
`end_effector_pose` / `centroidal_inner`), then assembles the geometric Jacobian
of a chosen target joint at runtime for one of the three pinocchio reference
frames:

    LOCAL (0)               -- twist in the frame's own body axes
    WORLD (1)               -- spatial Jacobian at the world origin
    LOCAL_WORLD_ALIGNED (2) -- at the frame origin, world-aligned axes

Output J is 6 x NUM_VEL column-major, rows ordered [linear(3); angular(3)] to
match pinocchio's `getFrameJacobian` / `getJointJacobian`. A second device
function composes Lambda = (J Minv J^T)^{-1}, the 6x6 operational-space inertia.

This is a direct CUDA transcription of the RBDReference numpy oracle
(`RBDReference.frame_jacobian` / `.osc_inertia`), which agrees with pinocchio to
float64 rounding. Correctness-first (single-block, serial inner assembly); the
per-column independence is left for a future perf pass.
"""

import numpy as np


__all__ = [
    "gen_frame_jacobian_inner",
    "gen_frame_jacobian_device",
    "gen_frame_jacobian",
]


# Reference-frame enum (matches pin.ReferenceFrame ordering used by the host).
_REF_LOCAL = 0
_REF_WORLD = 1
_REF_LWA = 2


def _frame_jacobian_inner_temp_mem_size(self):
    NJ = self.robot.get_num_joints()
    # world homogeneous transform per joint (16 each).
    return 16 * NJ


def gen_frame_jacobian_inner(self):
    """Emit frame_jacobian_inner: build per-joint world transforms then the
    geometric Jacobian (6 x NV, [linear; angular]) of `target_jid` in
    `reference_frame`. Correctness-first single-block assembly."""
    NJ = self.robot.get_num_joints()
    nv = self.robot.get_num_vel()
    n_bfs_levels = self.robot.get_max_bfs_level() + 1

    func_params = [
        "s_J is the output 6 x NUM_VEL geometric Jacobian (column-major, [linear; angular])",
        "target_jid is the joint id whose frame Jacobian is requested",
        "reference_frame is 0=LOCAL, 1=WORLD, 2=LOCAL_WORLD_ALIGNED",
        "s_q is the vector of joint positions (unused; baked into s_Xhom)",
        "s_Xhom is the per-joint LOCAL homogeneous transforms",
        "d_robotModel is the GPU model helpers",
        "s_temp is scratch of size " + str(_frame_jacobian_inner_temp_mem_size(self))]
    func_def_middle = ("T *s_J, const int target_jid, const int reference_frame, "
                       "const T *s_q, const T *s_Xhom, const robotModel<T> *d_robotModel, ")
    func_def = "void frame_jacobian_inner(" + func_def_middle + "T *s_temp) {"
    self.gen_add_func_doc("Compute a general-frame geometric Jacobian", [], func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    self.gen_add_code_line("(void)s_q;")

    self.gen_add_code_line("T *s_Xworld = &s_temp[0];")

    # ---- Step 1: world homogeneous transforms by BFS level (chain-up) ----
    self.gen_add_code_line("// Step 1: world homogeneous transforms (chain-up of local s_Xhom)")
    for level in range(n_bfs_levels):
        ids_at_level = self.robot.get_ids_by_bfs_level(level)
        if not ids_at_level:
            continue
        njs = len(ids_at_level)
        self.gen_add_parallel_loop("ind", str(16 * njs))
        self.gen_add_code_line("int slot = ind / 16; int ele = ind % 16;")
        self.gen_add_code_line("int row = ele & 3; int col = ele >> 2;")
        jid_list = [str(j) for j in ids_at_level]
        par_list = [str(self.robot.get_parent_id(j)) for j in ids_at_level]
        if njs > 1:
            self.gen_add_multi_threaded_select("slot", "<", [str(i + 1) for i in range(njs)],
                                               [("int", "jid", jid_list), ("int", "par", par_list)])
        else:
            self.gen_add_code_line("const int jid = " + jid_list[0] + "; const int par = " + par_list[0] + ";")
        self.gen_add_code_line("if (par == -1) { s_Xworld[16*jid + ele] = s_Xhom[16*jid + ele]; }")
        self.gen_add_code_line("else { s_Xworld[16*jid + ele] = dot_prod<T,4,4,1>(&s_Xworld[16*par + row], &s_Xhom[16*jid + 4*col]); }")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    # ---- Step 2: zero the Jacobian ----
    self.gen_add_parallel_loop("ind", str(6 * nv))
    self.gen_add_code_line("s_J[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # ---- Step 3: assemble J at the frame ORIGIN with WORLD axes ----
    # Bake the per-target-joint chain jobs at codegen time, behind a runtime
    # `target_jid` switch. For each (target jid, ancestor-or-self jj, S-col c)
    # the world contribution to J[:, vi] is the screw of joint jj evaluated at
    # the frame origin p_f (= s_Xworld[16*target_jid + 12..14]).
    self.gen_add_code_line("// Step 3: geometric Jacobian at frame origin, world axes")
    jobs = []  # (target_jid, jj, vi, ang_local[3], lin_local[3])
    for jid in range(NJ):
        if self.robot.get_parent_id(jid) == -1 and jid != 0:
            pass
        chain = sorted(self.robot.get_ancestors_by_id(jid)) + [jid]
        for jj in chain:
            S = np.asarray(self.robot.get_S_by_id(jj), dtype=np.float64)
            if S.ndim == 1:
                S = S.reshape(-1, 1)
            vinds = self.robot.get_joint_index_v(jj)
            if not isinstance(vinds, (list, tuple, np.ndarray)):
                vinds = [vinds]
            else:
                vinds = list(vinds)
            for c in range(S.shape[1]):
                vi = vinds[c] if c < len(vinds) else vinds[-1]
                ang = [float(S[0, c]), float(S[1, c]), float(S[2, c])]
                lin = [float(S[3, c]), float(S[4, c]), float(S[5, c])]
                jobs.append((jid, jj, vi, ang, lin))

    njobs = len(jobs)
    if njobs > 0:
        tj_arr = ", ".join(str(j[0]) for j in jobs)
        jj_arr = ", ".join(str(j[1]) for j in jobs)
        vi_arr = ", ".join(str(j[2]) for j in jobs)
        ax_arr = ", ".join("static_cast<T>({:.17g})".format(v) for j in jobs for v in j[3])
        lx_arr = ", ".join("static_cast<T>({:.17g})".format(v) for j in jobs for v in j[4])
        self.gen_add_code_line("const int fj_target[" + str(njobs) + "] = {" + tj_arr + "};")
        self.gen_add_code_line("const int fj_jj[" + str(njobs) + "] = {" + jj_arr + "};")
        self.gen_add_code_line("const int fj_vi[" + str(njobs) + "] = {" + vi_arr + "};")
        self.gen_add_code_line("const T fj_ang[" + str(3 * njobs) + "] = {" + ax_arr + "};")
        self.gen_add_code_line("const T fj_lin[" + str(3 * njobs) + "] = {" + lx_arr + "};")
        # Serial accumulation (correctness-first; columns may repeat vi).
        self.gen_add_serial_ops()
        # frame origin p_f in world.
        self.gen_add_code_line("const T *Xf = &s_Xworld[16*target_jid];")
        self.gen_add_code_line("T pfx = Xf[12], pfy = Xf[13], pfz = Xf[14];")
        self.gen_add_code_line("for (int t = 0; t < " + str(njobs) + "; ++t) {", True)
        self.gen_add_code_line("if (fj_target[t] != target_jid) continue;")
        self.gen_add_code_line("int jj = fj_jj[t]; int vi = fj_vi[t];")
        self.gen_add_code_line("const T *Xj = &s_Xworld[16*jj];")
        self.gen_add_code_line("T a0=fj_ang[3*t], a1=fj_ang[3*t+1], a2=fj_ang[3*t+2];")
        self.gen_add_code_line("T l0=fj_lin[3*t], l1=fj_lin[3*t+1], l2=fj_lin[3*t+2];")
        # world axes: aw = R_jj * ang_local ; lw = R_jj * lin_local
        self.gen_add_code_line("T aw0 = Xj[0]*a0 + Xj[4]*a1 + Xj[8]*a2;")
        self.gen_add_code_line("T aw1 = Xj[1]*a0 + Xj[5]*a1 + Xj[9]*a2;")
        self.gen_add_code_line("T aw2 = Xj[2]*a0 + Xj[6]*a1 + Xj[10]*a2;")
        self.gen_add_code_line("T lw0 = Xj[0]*l0 + Xj[4]*l1 + Xj[8]*l2;")
        self.gen_add_code_line("T lw1 = Xj[1]*l0 + Xj[5]*l1 + Xj[9]*l2;")
        self.gen_add_code_line("T lw2 = Xj[2]*l0 + Xj[6]*l1 + Xj[10]*l2;")
        # p_jj (world origin of joint jj)
        self.gen_add_code_line("T pjx = Xj[12], pjy = Xj[13], pjz = Xj[14];")
        # linear at the FRAME origin = lw + aw x (p_f - p_jj)
        self.gen_add_code_line("T dx = pfx - pjx, dy = pfy - pjy, dz = pfz - pjz;")
        self.gen_add_code_line("T linf0 = lw0 + (aw1*dz - aw2*dy);")
        self.gen_add_code_line("T linf1 = lw1 + (aw2*dx - aw0*dz);")
        self.gen_add_code_line("T linf2 = lw2 + (aw0*dy - aw1*dx);")
        # accumulate into J[:, vi] ([linear; angular] col-major 6 x NV)
        self.gen_add_code_line("T *Jc = &s_J[6*vi];")
        self.gen_add_code_line("Jc[0]+=linf0; Jc[1]+=linf1; Jc[2]+=linf2; Jc[3]+=aw0; Jc[4]+=aw1; Jc[5]+=aw2;")
        self.gen_add_end_control_flow()  # for t
        self.gen_add_end_control_flow()  # serial
        self.gen_add_sync()

    # ---- Step 4: apply the reference-frame transform in place ----
    # s_J currently holds the LOCAL_WORLD_ALIGNED Jacobian (frame origin, world
    # axes). Convert to WORLD or LOCAL when requested.
    self.gen_add_code_line("// Step 4: reference-frame transform (in place per column)")
    self.gen_add_serial_ops()
    self.gen_add_code_line("const T *Xf2 = &s_Xworld[16*target_jid];")
    self.gen_add_code_line("T Rf[9]; for (int c=0;c<3;++c) for (int r=0;r<3;++r) Rf[r+3*c] = Xf2[r + 4*c];")
    self.gen_add_code_line("T pfx2 = Xf2[12], pfy2 = Xf2[13], pfz2 = Xf2[14];")
    self.gen_add_code_line("for (int vi = 0; vi < " + str(nv) + "; ++vi) {", True)
    self.gen_add_code_line("T *Jc = &s_J[6*vi];")
    self.gen_add_code_line("T v0=Jc[0], v1=Jc[1], v2=Jc[2], w0=Jc[3], w1=Jc[4], w2=Jc[5];")
    # if/else-if chain emitted as raw lines (no auto control-flow bookkeeping).
    self.gen_add_code_lines([
        "if (reference_frame == " + str(_REF_WORLD) + ") {",
        # spatial Jacobian at world origin: v_world = v_lwa + p_f x w ; w unchanged.
        "  Jc[0] = v0 + (pfy2*w2 - pfz2*w1);",
        "  Jc[1] = v1 + (pfz2*w0 - pfx2*w2);",
        "  Jc[2] = v2 + (pfx2*w1 - pfy2*w0);",
        "} else if (reference_frame == " + str(_REF_LOCAL) + ") {",
        # rotate both blocks into the frame body axes: Jc <- Rf^T Jc
        "  Jc[0] = Rf[0]*v0 + Rf[1]*v1 + Rf[2]*v2;",
        "  Jc[1] = Rf[3]*v0 + Rf[4]*v1 + Rf[5]*v2;",
        "  Jc[2] = Rf[6]*v0 + Rf[7]*v1 + Rf[8]*v2;",
        "  Jc[3] = Rf[0]*w0 + Rf[1]*w1 + Rf[2]*w2;",
        "  Jc[4] = Rf[3]*w0 + Rf[4]*w1 + Rf[5]*w2;",
        "  Jc[5] = Rf[6]*w0 + Rf[7]*w1 + Rf[8]*w2;",
        "}",
    ])
    self.gen_add_end_control_flow()  # for vi
    self.gen_add_end_control_flow()  # serial
    self.gen_add_sync()
    self.gen_add_end_function()


def gen_frame_jacobian_device(self):
    """Auto-smem device wrapper around frame_jacobian_inner."""
    nv = self.robot.get_num_vel()
    func_def = ("void frame_jacobian_device(T *s_J, const int target_jid, const int reference_frame, "
                "const T *s_q, const robotModel<T> *d_robotModel) {")
    func_params = ["s_J holds the 6 x NUM_VEL geometric Jacobian (column-major, [linear; angular])",
                   "target_jid is the joint id of the frame",
                   "reference_frame is 0=LOCAL, 1=WORLD, 2=LOCAL_WORLD_ALIGNED",
                   "s_q is the joint position vector",
                   "d_robotModel is the GPU model helpers"]
    self.gen_add_func_doc("Compute a general-frame geometric Jacobian", [], func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    # s_J is caller-provided (function param); the arena only holds the XmatsHom
    # world-transform machinery + inner scratch.
    self.gen_XmatsHom_helpers_temp_shared_memory_code(
        _frame_jacobian_inner_temp_mem_size(self),
        include_linalg_scratch=True, linalg_scratch_bytes="GRID_EE_LINALG_SHARED_BYTES<T>()")
    self.gen_load_update_XmatsHom_helpers_function_call()
    self.gen_add_code_line("frame_jacobian_inner<T>(s_J, target_jid, reference_frame, s_q, s_XmatsHom, d_robotModel, s_temp);")
    self.gen_add_sync()
    self.gen_add_end_function()


def gen_frame_jacobian(self):
    self.gen_frame_jacobian_inner()
    self.gen_frame_jacobian_device()
