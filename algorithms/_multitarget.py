"""General batched multi-target kinematics emitter (W1b live; W2a gradient prototype).

STATUS:
  * POSITION (build_target_batch + gen_multi_target_position{,_inner,_device,...}) is
    LIVE — wired in algorithms/__init__.py + the GRiDCodeGenerator class import list,
    dispatched from gen_all_code behind the opt-in `multi_target_batch` kwarg (default
    None -> not emitted, all existing robots byte-identical). Validated on iiwa14 +
    baxter (test/cuda_equivalents/test_cuda_multi_target_position.py): world positions
    vs a NumPy FK oracle, offset==0 == end_effector_pose, thread-invariant, sanitizers.
    Kernel/host/gridData-buffer registration + algo_registry descriptor = W1b.3.
  * GRADIENT (gen_multi_target_gradient_phaseB) is still a PROTOTYPE for W2a (the offset
    epilogue dpos = Jv + Jw x (R.r), no FK re-walk). Not wired.

Design: docs/open-tasks/design_W1b_batched_multitarget_position_2026-07-07.md.

A TARGET = (anchor_jid, offset[3]) — a fixed point off a link (named grasp point,
gripper, or collision sphere). A BATCH = an ordered target list + named sub-ranges
("all" + groups). World position pos[t] = X_world[anchor(t)] @ [offset,1], computed
by ONE shared FK (all link world transforms) + a cheap parallel-over-targets
extraction driven by a baked (anchor, offset) table — same table-driven idiom as the
W1a hessian collapse. Subsumes backlog D (multi-named-EE-target).
"""


# ---------------------------------------------------------------------------
# Build-time: assemble the target batch (anchors + baked offsets + named groups)
# ---------------------------------------------------------------------------
def build_target_batch(self, targets):
    """Given `targets` = ordered list of dicts {"anchor_jid": int, "offset": (x,y,z),
    "group": str (optional)}, return a flat batch descriptor:
        {
          "n": N,
          "anchor": [anchor_jid ...]          # len N, GRiD joint/frame index
          "offset": [x,y,z, x,y,z, ...]       # len 3N, LOCAL frame, near-zero snapped
          "groups": {name: (lo, hi)}          # contiguous [lo,hi) index ranges + "all"
        }
    Offsets are snapped |c|<1e-15 -> 0.0 (bit-identical world-axis emission, matches the
    axis-literal snap in _eepose_gradient_hessian). Targets are grouped contiguously so a
    caller can request "all" (whole [0,N)) or a named group's slice.
    NOTE (W3/FLANGE): for foam spheres, `anchor_jid` MUST be resolved through GRiD's own
    URDFParser link->frame id, NOT copied from foam's actuated-joint-count table.
    """
    def _snap(v):
        return [float(c) if abs(c) >= 1e-15 else 0.0 for c in v]

    # stable-sort by group so each group is a contiguous slice; ungrouped -> "_" first.
    ordered = sorted(range(len(targets)), key=lambda i: targets[i].get("group", ""))
    anchor, offset, groups = [], [], {}
    for new_idx, i in enumerate(ordered):
        t = targets[i]
        anchor.append(int(t["anchor_jid"]))
        offset.extend(_snap(t["offset"]))
        g = t.get("group", "")
        if g:
            lo, hi = groups.get(g, (new_idx, new_idx))
            groups[g] = (min(lo, new_idx), new_idx + 1)
    n = len(anchor)
    groups["all"] = (0, n)
    return {"n": n, "anchor": anchor, "offset": offset, "groups": groups}


# NOTE: the shared BFS world-FK chain-up lives in _eepose_gradient_hessian.py as
# emit_world_fk_chainup (committed, byte-identical refactor of the gradient inner's
# Steps 1+1b). The position/gradient emitters below import and call THAT — there is
# no separate copy here (the earlier prototype emit_shared_world_fk was superseded).


# ---------------------------------------------------------------------------
# Batched multi-target POSITION: scratch sizing.
# ---------------------------------------------------------------------------
def gen_multi_target_position_inner_temp_mem_size(self, batch=None):
    """Helper scratch for the position inner = the s_Xworld world-transform arena,
    16 per joint (+ appended fixed-anchor slots when welded targets are baked in).
    Mirrors _eepose_xworld_slot_count in _eepose_gradient_hessian.py so the two
    kinematics paths size the shared FK arena identically. The batch's output buffer
    (s_out_pos = 3*N) is allocated by the kernel/device wrapper, NOT counted here
    (same split as end_effector_pose: s_temp holds only the chain scratch)."""
    from ._eepose_gradient_hessian import _eepose_xworld_slot_count
    return 16 * _eepose_xworld_slot_count(self)


# ---------------------------------------------------------------------------
# Batched multi-target POSITION inner.
# ---------------------------------------------------------------------------
def gen_multi_target_position_inner(self, batch):
    """Emit `multi_target_position_inner<T>`: compute ALL targets' world positions in one
    call. s_out_pos is 3*N (xyz per target). One shared FK (s_Xworld) + a parallel
    extraction over the baked (anchor, offset) table. `batch` = build_target_batch(...) output.

    pos[t][row] = R_world[anchor]·offset + p_world
                = X[row]*o0 + X[row+4]*o1 + X[row+8]*o2 + X[row+12]   (X col-major 4x4)
    """
    n = batch["n"]
    n_joints = self.robot.get_num_joints()  # s_Xworld slot count (extend for fixed anchors as needed)

    # --- function boilerplate (models gen_end_effector_pose_inner) ---
    func_params = [
        "s_out_pos is shared memory of size 3*N_TARGETS (xyz per target), N_TARGETS = " + str(n),
        "s_q is the vector of joint positions",
        "s_Xhom is the per-joint local homogeneous transforms (already updated for q)",
        "s_temp is helper shared memory (holds s_Xworld = 16*NUM_JOINTS)",
        "d_workspace is the global-memory scratch used when !TEMP_IN_SMEM",
    ]
    func_notes = [
        "Computes world positions of a baked batch of fixed-offset targets (grasp points / spheres).",
        "One shared FK (world transforms) + parallel-over-targets offset extraction.",
    ]
    func_def_start = "void multi_target_position_inner("
    func_def_middle = "T *s_out_pos, const T *s_q, const T *s_Xhom, "
    func_def_end = "T *s_temp, T *d_workspace, unsigned char *s_linalg_smem) {"
    func_def_middle, func_params = self.gen_insert_helpers_func_def_params(
        func_def_middle, func_params, -1, NO_XI_FLAG=True)
    self.gen_add_func_doc("Batched multi-target world positions", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T, bool TEMP_IN_SMEM = true>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def_start + func_def_middle + func_def_end, True)
    self.gen_add_code_line("if constexpr (!TEMP_IN_SMEM) { s_temp = d_workspace; } else { (void)d_workspace; }")
    self.gen_add_code_line("(void)s_q; (void)s_linalg_smem;")
    self.gen_add_code_line("T *s_Xworld = s_temp;   // 16 * " + str(n_joints))

    # --- shared FK (the committed emit_world_fk_chainup in _eepose_gradient_hessian) ---
    from ._eepose_gradient_hessian import emit_world_fk_chainup
    emit_world_fk_chainup(
        self,
        header_lines=["//",
                      "// Build world transforms for every joint via BFS-level chain-up",
                      "//"],
        fixed_anchors=None)  # welded (fixed-joint) anchors: pass their (anchor,parent) here (W3)

    # --- baked batch tables ---
    self.gen_add_code_line("// baked target batch: anchor frame id + LOCAL offset per target")
    self.gen_add_code_line("static const int mt_anchor[" + str(n) + "] = {" +
                           ", ".join(str(a) for a in batch["anchor"]) + "};")
    self.gen_add_code_line("const T mt_offset[" + str(3 * n) + "] = {" +
                           ", ".join("static_cast<T>({:.17g})".format(v) for v in batch["offset"]) + "};")

    # --- parallel extraction: one thread per (target, xyz) ---
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Extract each target's world position = R_world[anchor] @ offset + p_world")
    self.gen_add_code_line("//")
    self.gen_add_parallel_loop("ind", str(3 * n))
    self.gen_add_code_line("int row = ind % 3; int t = ind / 3;")
    self.gen_add_code_line("const T *X = &s_Xworld[16 * mt_anchor[t]];")
    self.gen_add_code_line("const T *o = &mt_offset[3 * t];")
    self.gen_add_code_line("s_out_pos[3*t + row] = X[row]*o[0] + X[row + 4]*o[1] + X[row + 8]*o[2] + X[row + 12];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# Inner function-call helper (mirrors gen_end_effector_pose_inner_function_call).
# ---------------------------------------------------------------------------
def gen_multi_target_position_inner_function_call(self, updated_var_names=None, temp_in_smem_expr="true"):
    var_names = dict(
        s_Xhom_name="s_XmatsHom",
        s_out_pos_name="s_out_pos",
        s_q_name="s_q",
        s_topology_helpers_name="s_topology_helpers",
        s_temp_name="s_temp",
        d_workspace_name="nullptr",
        s_linalg_smem_name="s_linalg_smem",
    )
    if updated_var_names is not None:
        for key, value in updated_var_names.items():
            var_names[key] = value
    code_start = ("multi_target_position_inner<T, " + temp_in_smem_expr + ">(" +
                  var_names["s_out_pos_name"] + ", " + var_names["s_q_name"] + ", ")
    code_middle = var_names["s_Xhom_name"] + ", "
    code_end = var_names["s_temp_name"] + ", " + var_names["d_workspace_name"] + ", " + var_names["s_linalg_smem_name"] + ");"
    # NO_XI: this family takes s_Xhom (not s_XImats); mirror the def's NO_XI_FLAG.
    code_middle += self.gen_insert_helpers_function_call(updated_var_names=var_names, NO_XI_FLAG=True)
    self.gen_add_code_line(code_start + code_middle + code_end)


# ---------------------------------------------------------------------------
# Device wrapper (models gen_end_effector_pose_device): allocates the shared
# XmatsHom/temp arena, builds the local per-joint transforms for q, then calls
# the batched inner. Caller supplies the s_out_pos output buffer (3*N shared),
# exactly like end_effector_pose_device(s_pose, d_q, m). No gridData dependency.
# ---------------------------------------------------------------------------
def gen_multi_target_position_device(self, batch):
    n = batch["n"]
    func_params = [
        "s_out_pos is a pointer to shared memory of size 3*N_TARGETS where N_TARGETS = " + str(n),
        "s_q is the vector of joint positions",
        "d_robotModel is the pointer to the initialized model specific helpers on the GPU (XImats, topology_helpers, etc.)"]
    func_notes = ["Computes world positions of a baked batch of fixed-offset targets (grasp points / spheres)."]
    func_def_start = "void multi_target_position_device("
    func_def_middle = "T *s_out_pos, const T *s_q, "
    func_def_end = "const robotModel<T> *d_robotModel) {"
    func_def = func_def_start + func_def_middle + func_def_end
    self.gen_add_func_doc("Computes batched multi-target world positions", func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    shared_mem_size = self.gen_multi_target_position_inner_temp_mem_size(batch)
    self.gen_XmatsHom_helpers_temp_shared_memory_code(shared_mem_size, include_linalg_scratch=True,
                                                      linalg_scratch_bytes="GRID_EE_LINALG_SHARED_BYTES<T>()")
    self.gen_load_update_XmatsHom_helpers_function_call()
    self.gen_multi_target_position_inner_function_call()
    self.gen_add_end_function()


# ---------------------------------------------------------------------------
# Dispatcher — W1b.2 scope: inner + device (validated via a focused runner that
# calls multi_target_position_device directly). Kernel + host + gridData buffer
# registration = W1b.3.
# ---------------------------------------------------------------------------
def gen_multi_target_position(self, batch):
    n = batch["n"]
    XHom_size, _dXhom, _d2Xhom = self.gen_get_Xhom_size()
    total_t = XHom_size + self.gen_multi_target_position_inner_temp_mem_size(batch)
    self.gen_add_code_lines([
        "// W1b batched multi-target world positions (opt-in via multi_target_batch); NUM_MULTI_TARGETS = " + str(n),
        "const int NUM_MULTI_TARGETS = " + str(n) + ";",
        "template <typename T> __host__ __device__ inline size_t MULTI_TARGET_POSITION_DYNAMIC_SHARED_MEM_BYTES() "
        "{ return grid_shared_arena_bytes<T>(" + str(total_t) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
    ])
    self.gen_multi_target_position_inner(batch)
    self.gen_multi_target_position_device(batch)


# ---------------------------------------------------------------------------
# Batched multi-target GRADIENT — Phase B (the novel offset epilogue).
# Design: docs/open-tasks/design_W2a_batched_multitarget_gradient_2026-07-07.md
#
# Phase A (NOT prototyped here) = build s_Jv, s_Jw (3 x nv per DISTINCT anchor) by
# reusing gen_end_effector_pose_gradient_inner Steps 1-3 with all_ees := distinct
# anchors (ee-index -> anchor_idx). That already bakes the per-(anchor, S-col) job
# table (eeg_job_*) and fills Jv = axw x (p_anchor - p_j), Jw = axw (revolute) etc.,
# incl. the mimic alpha-fold. A __syncthreads follows (Jv/Jw complete).
#
# Phase B (below) = per-(target, vi, row) offset-corrected writeback, driven by the
# baked batch tables + a target->anchor_idx map. NO FK re-walk:
#   dpos[t][:,j] = Jv[anchor,:,j] + Jw[anchor,:,j] x (R_world[anchor] . r_local)
# ---------------------------------------------------------------------------
def gen_multi_target_gradient_phaseB(self, batch, anchor_index_of_target):
    """Emit the offset epilogue. Assumes s_Jv/s_Jw (3 x nv per distinct anchor) are built
    (Phase A) and s_Xworld holds anchor world transforms. `anchor_index_of_target[t]` maps
    a target to its slot in the deduped anchor set; `batch` from build_target_batch.
    Output s_out_grad is 3 * nv * N_targets (position gradient rows only)."""
    n = batch["n"]
    nv = self.robot.get_num_vel()

    # baked: target -> anchor slot (into s_Jv/s_Jw) and target -> local offset.
    self.gen_add_code_line("static const int mt_anchor_idx[" + str(n) + "] = {" +
                           ", ".join(str(a) for a in anchor_index_of_target) + "};")
    self.gen_add_code_line("static const int mt_anchor[" + str(n) + "] = {" +
                           ", ".join(str(a) for a in batch["anchor"]) + "};")
    self.gen_add_code_line("const T mt_offset[" + str(3 * n) + "] = {" +
                           ", ".join("static_cast<T>({:.17g})".format(v) for v in batch["offset"]) + "};")

    # Pre-pass: rotate each target's LOCAL offset into world once (vi-independent).
    # s_ro is 3*N scratch. ro = R_world[anchor] @ r_local.
    self.gen_add_code_line("// ro[t] = R_world[anchor(t)] @ offset(t)  (once per target, vi-independent)")
    self.gen_add_parallel_loop("t", str(n))
    self.gen_add_code_line("const T *X = &s_Xworld[16 * mt_anchor[t]];")
    self.gen_add_code_line("const T *o = &mt_offset[3 * t];")
    self.gen_add_code_line("s_ro[3*t + 0] = X[0]*o[0] + X[4]*o[1] + X[8]*o[2];")
    self.gen_add_code_line("s_ro[3*t + 1] = X[1]*o[0] + X[5]*o[1] + X[9]*o[2];")
    self.gen_add_code_line("s_ro[3*t + 2] = X[2]*o[0] + X[6]*o[1] + X[10]*o[2];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # Main: per (target, vi) -> corrected 3-vector. dpos = Jv + Jw x ro.
    self.gen_add_code_line("// dpos[t][:,vi] = Jv[anchor,:,vi] + Jw[anchor,:,vi] x ro[t]")
    self.gen_add_parallel_loop("ind", str(n * nv))
    self.gen_add_code_line("int vi = ind % " + str(nv) + "; int t = ind / " + str(nv) + ";")
    self.gen_add_code_line("int jb = 3 * (" + str(nv) + " * mt_anchor_idx[t] + vi);")
    self.gen_add_code_line("T Jv0 = s_Jv[jb+0], Jv1 = s_Jv[jb+1], Jv2 = s_Jv[jb+2];")
    self.gen_add_code_line("T Jw0 = s_Jw[jb+0], Jw1 = s_Jw[jb+1], Jw2 = s_Jw[jb+2];")
    self.gen_add_code_line("T r0 = s_ro[3*t+0], r1 = s_ro[3*t+1], r2 = s_ro[3*t+2];")
    self.gen_add_code_line("int ob = 3 * (" + str(nv) + " * t + vi);")
    self.gen_add_code_line("s_out_grad[ob + 0] = Jv0 + (Jw1*r2 - Jw2*r1);")
    self.gen_add_code_line("s_out_grad[ob + 1] = Jv1 + (Jw2*r0 - Jw0*r2);")
    self.gen_add_code_line("s_out_grad[ob + 2] = Jv2 + (Jw0*r1 - Jw1*r0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()


# ---------------------------------------------------------------------------
# INTEGRATION STATUS
# ---------------------------------------------------------------------------
# [x] gen_multi_target_position_inner_temp_mem_size -> 16 * n_xworld_slots (+ fixed anchors)
# [x] gen_multi_target_position_device wrapper (models gen_end_effector_pose_device);
#     tier-awareness deferred to W2b.
# [x] shared world-FK: gradient inner Steps 1+1b factored into emit_world_fk_chainup
#     (_eepose_gradient_hessian.py), byte-identical gate PASSED (GCG cb73296); this
#     emitter imports and calls it.
# [x] wired into GRiDCodeGenerator.gen_all_code behind the opt-in multi_target_batch kwarg
#     (default None -> not emitted; existing robots byte-identical).
# [x] test: baxter (multi-anchor) + iiwa14 (offset==0==ee_pose) positions vs NumPy FK oracle;
#     thread-invariance 1/32/256 (bit-identical); synccheck/racecheck/memcheck clean.
# W1b.3 (remaining): _kernel/_host wrappers + gridData d_/h_ buffer + algo_registry
#     AlgoEntry/AlgoDescriptor rows + KERNEL_OVERLOADS + bench GRID_HAS_* wrappers.
# W2a (remaining): wire gen_multi_target_gradient_phaseB (offset epilogue, no FK re-walk).
