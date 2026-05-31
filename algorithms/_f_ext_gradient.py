"""External-force gradient column (section A of the differentiability plan).

Emits the three f_ext-gradient outputs (all q-only, f_ext-value independent —
f_ext enters RNEA additively & linearly so the Jacobian carries no f_ext value):

  dtau/dfext      = -J(q)^T          (section A.1)  stacked body-Jacobian transpose
  dqdd/dfext      =  M^{-1} J^T      (section A.2)  operational-space inverse-inertia
  d(id_du)/dfext  = -dJ^T/dq         (section A.3)  q-derivative of the body Jacobian

J^T is the stacked SPATIAL body-Jacobian transpose in each link's LOCAL frame:

  -J^T[v_j, 6*i + k] = -( S_j^T  P_{i,j} )_k   for j on path(root->i), else 0
        P_{i,j} = X[j+1]^T X[j+2]^T ... X[i]^T   (composed local 6x6 spatial
                                                  transforms, [angular;linear])

This is EXACTLY the matrix the RNEA backward force sweep applies (f[parent] +=
X^T f), so feeding a unit local wrench at link i and running the back-prop yields
column block i of -J^T. We emit it as the explicit per-(body, chain-joint) build
from s_XImats (the same 6x6 local transforms RNEA loads), one running 6x6 product
per body chained root-ward. The convention (LOCAL frame, SUBTRACTED -> sign -)
is locked to the T4 forward path (apply_external_forces: f[:,i] -= f_ext[i]).

The output layout is body-major: s_dtau_dfext is nv x (6*NB), column-major in the
[v_row + nv*col] sense used by the rest of GRiD's dense gradient outputs.

dqdd/dfext is -s_Minv @ s_dtau_dfext (a single nv x nv * nv x 6NB GEMM reusing the
direct_minv s_Minv). dJ^T/dq is central-FD of the analytic -J^T over each
generalized coordinate (the same FD-on-Jacobian strategy the d2ee GPU path uses);
the q-dot block is identically zero (J^T is q-only) and is not stored.
"""


def _f_ext_grad_chain_jobs(self):
    """Bake the per-(body i, chain-joint j, S-column) fill jobs for -J^T.

    Returns (NB, nv, jobs) where each job is a dict:
      { 'i': body id, 'j': chain joint id, 'vi': velocity slot, 'Scol': the 6-vec
        motion subspace column, 'tf_chain': the ordered joint ids [j+1, ..., i]
        whose local 6x6 motion transforms X[m] are applied (left-fold) to Scol to
        push it from joint j's frame down to body i's frame. }

    The body-Jacobian column of body i for chain joint j is
      col = X[i] X[i-1] ... X[j+1] S_j   (Featherstone motion transforms),
    written to row v_j, column-block i of -J^T (negated). Out-of-chain columns are
    absent here (left zero by the inner's init).
    """
    import numpy as _np
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    jobs = []
    for i in range(NB):
        chain = sorted(self.robot.get_ancestors_by_id(i)) + [i]
        for j in chain:
            S = _np.asarray(self.robot.get_S_by_id(j), dtype=_np.float64)
            if S.ndim == 1:
                S = S.reshape(-1, 1)
            try:
                vinds = self.robot.get_joint_index_v(j)
            except Exception:
                vinds = self.robot.get_joint_index_q(j)
            if not isinstance(vinds, (list, tuple, _np.ndarray)):
                vinds = [vinds]
            vinds = list(vinds)
            # transform chain: the body-Jacobian column of body i for chain joint
            # j is the motion subspace S_j transformed from joint j's frame DOWN to
            # body i's frame: col = X[i] X[i-1] ... X[j+1] S_j  (Featherstone motion
            # transforms, [angular;linear]). Applying as a left-fold over a running
            # 6-vector means apply X[m] (NO transpose) for m = j+1, j+2, ..., i. We
            # bake the chain in that (root-ward-reversed) order. This is exactly the
            # transpose of the RNEA backward force sweep's P_{i,j} = X[j+1]^T..X[i]^T
            # (verified bit-exact vs the rnea_bpass unit-wrench oracle).
            tf_chain = []
            m = i
            while m != j:
                tf_chain.append(int(m))
                m = self.robot.get_parent_id(m)
            tf_chain = list(reversed(tf_chain))  # j+1, j+2, ..., i
            for c in range(S.shape[1]):
                vi = vinds[c] if c < len(vinds) else vinds[-1]
                Scol = [float(x) for x in S[:6, c]]
                jobs.append({
                    "i": int(i), "j": int(j), "vi": int(vi),
                    "Scol": Scol, "tf_chain": tf_chain,
                })
    return NB, nv, jobs


def gen_f_ext_gradient_inner_temp_mem_size(self):
    # scratch: per work-item 6x6 product buffer is built serially via two 36-slot
    # double buffers shared across the block. We allocate 2 * 36 * NB so each body
    # i has its own running product (parallel across bodies). Plus a 6-vec staging
    # per job is folded into the output directly.
    NB = self.robot.get_num_bodies()
    return 2 * 36 * NB


def gen_f_ext_gradient_inner_function_call(self, updated_var_names=None):
    var_names = dict(
        s_dtau_dfext_name="s_dtau_dfext",
        s_q_name="s_q",
        s_temp_name="s_temp",
    )
    if updated_var_names is not None:
        for k, v in updated_var_names.items():
            var_names[k] = v
    code_start = "f_ext_gradient_jacobianT_inner<T>(" + var_names["s_dtau_dfext_name"] + ", " + var_names["s_q_name"] + ", "
    code_mid = self.gen_insert_helpers_function_call()
    code_end = var_names["s_temp_name"] + ");"
    self.gen_add_code_line(code_start + code_mid + code_end)


def gen_f_ext_gradient_jacobianT_inner(self):
    """Emit f_ext_gradient_jacobianT_inner: builds -J^T into s_dtau_dfext.

    s_dtau_dfext is nv x (6*NB), zeroed then filled per chain job. Assumes
    s_XImats holds the per-joint LOCAL 6x6 spatial transforms for the current q.
    """
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    _, _, jobs = _f_ext_grad_chain_jobs(self)

    func_params = [
        "s_dtau_dfext is the output dtau/dfext = -J^T, size NV*(6*NB) = " + str(nv * 6 * NB),
        "s_q is the vector of joint positions (used only via s_XImats)",
        "s_temp is helper shared memory of size " + str(self.gen_f_ext_gradient_inner_temp_mem_size()),
    ]
    func_notes = [
        "Assumes s_XImats is updated already for the current s_q.",
        "Output is the LOCAL-frame stacked body-Jacobian transpose, negated.",
        "Column block i (6 cols) is the joint-torque response to a unit local",
        "wrench on body i; nonzero only on rows v_j for j on path(root->i).",
    ]
    func_def_start = "void f_ext_gradient_jacobianT_inner(T *s_dtau_dfext, const T *s_q, "
    func_def_end = "T *s_temp) {"
    func_def_start, func_params = self.gen_insert_helpers_func_def_params(func_def_start, func_params, -1)
    func_def = func_def_start + func_def_end

    self.gen_add_func_doc("Computes dtau/dfext = -J(q)^T (stacked local body-Jacobian transpose)",
                          func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    self.gen_add_code_line("(void)s_q;")

    out_size = nv * 6 * NB
    # zero the output
    self.gen_add_code_line("// zero the full nv x 6*NB output (out-of-chain cols stay zero)")
    self.gen_add_parallel_loop("ind", str(out_size))
    self.gen_add_code_line("s_dtau_dfext[ind] = static_cast<T>(0);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()

    # For each job (body i, chain joint j): the -J^T column block is -col, where
    #   col = X[i] X[i-1] ... X[j+1] S_j   (push the motion subspace S_j from joint
    # j's frame down to body i's frame via the local 6x6 motion transforms X[m]).
    # We build col as a left-fold over a running 6-vector: start col = S_j, then
    # apply col := X[m] @ col for m = j+1, j+2, ..., i (tf_chain order). This is the
    # transpose of the RNEA backward force sweep's S_j^T (X[j+1]^T...X[i]^T) and was
    # verified bit-exact vs the rnea_bpass unit-wrench oracle (and the GPU emit vs
    # the numpy + pinocchio f_ext_gradient oracle on iiwa14-fixed and go2-floating).
    #
    # The result -col is written to output row v_j, column block i. Jobs that share
    # (i, v_j) accumulate (+=) -- matching the reference's S-column / mimic v-slot
    # fold. To keep the parallel slab-fill race-free the final += reduction is run
    # serially on lane 0 (the slab itself is filled fully in parallel).
    self.gen_add_code_line("//")
    self.gen_add_code_line("// Per chain job: col(v_j, body i) = -(X[i]..X[j+1] S_j) (motion-transform pushdown)")
    self.gen_add_code_line("//")

    if len(jobs) > 0:
        # Bake the per-job chains as a flat const int array with [start,len]
        # offsets, plus the per-job (body i, output row v_j) and the 6-vec S column.
        flat_chain = []
        job_off = []
        job_len = []
        job_i = []
        job_vrow = []
        job_S = []
        for job in jobs:
            job_off.append(len(flat_chain))
            job_len.append(len(job["tf_chain"]))
            flat_chain.extend(job["tf_chain"])
            job_i.append(job["i"])
            job_vrow.append(job["vi"])
            job_S.append(job["Scol"])

        def _ints(vals):
            return "{ " + ", ".join(str(v) for v in vals) + " }"

        def _floats(vals):
            return "{ " + ", ".join("static_cast<T>({:.17g})".format(v) for v in vals) + " }"

        if len(flat_chain) == 0:
            flat_chain = [0]  # avoid zero-size array
        self.gen_add_code_line("static const int feg_chain[]   = " + _ints(flat_chain) + ";")
        self.gen_add_code_line("static const int feg_job_off[]  = " + _ints(job_off) + ";")
        self.gen_add_code_line("static const int feg_job_len[]  = " + _ints(job_len) + ";")
        self.gen_add_code_line("static const int feg_job_i[]    = " + _ints(job_i) + ";")
        self.gen_add_code_line("static const int feg_job_vrow[] = " + _ints(job_vrow) + ";")
        flatS = [s for job in job_S for s in job]
        self.gen_add_code_line("const T feg_job_S[] = " + _floats(flatS) + ";")

        njobs = len(jobs)
        # Two-phase to avoid += races on shared (i, v_j) destinations: (1) each job
        # (one work-item) computes its -col 6-vector chain in parallel into a per-job
        # scratch slab (njobs*6 in s_temp); (2) a lane-0 serial reduce sums each
        # job's slab into the output row v_j / column-block i (folding shared v-slots
        # in deterministic order). The chain matvec is a cheap short serial loop.
        slab = njobs * 6
        self.gen_add_code_line("// per-job contribution slab in s_temp (njobs*6)")
        self.gen_add_code_line("T *s_feg_slab = s_temp;   // size " + str(slab))
        self.gen_add_parallel_loop("jb", str(njobs))
        self.gen_add_code_line("int off = feg_job_off[jb]; int len = feg_job_len[jb];")
        self.gen_add_code_line("T col[6];")
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int r = 0; r < 6; ++r) { col[r] = feg_job_S[6*jb + r]; }")
        # chain: for each m in feg_chain[off..off+len): col := X[m] @ col
        self.gen_add_code_line("for (int s = 0; s < len; ++s) {", True)
        self.gen_add_code_line("int m = feg_chain[off + s];")
        self.gen_add_code_line("const T *X = &s_XImats[36*m];")
        self.gen_add_code_line("T tmp[6];")
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int r = 0; r < 6; ++r) {", True)
        # col := X @ col (motion transform, NO transpose). X is column-major 6x6:
        # element (row r, col c) = X[6*c + r]. So (X @ col)[r] = sum_c X[6*c + r] col[c]
        # = dot_prod with stride 6 on X starting at r, stride 1 on col.
        self.gen_add_code_line("tmp[r] = dot_prod<T,6,6,1>(&X[r], col);")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int r = 0; r < 6; ++r) { col[r] = tmp[r]; }")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int r = 0; r < 6; ++r) { s_feg_slab[6*jb + r] = -col[r]; }")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

        # reduce slab -> output: for each job add its 6 cells into
        # s_dtau_dfext[vrow + nv*(6*i + r)]  (column-major nv-row layout)
        self.gen_add_code_line("// reduce per-job contributions into the output (serial over jobs to fold shared v-slots)")
        self.gen_add_serial_ops()
        self.gen_add_code_line("for (int jb = 0; jb < " + str(njobs) + "; ++jb) {", True)
        self.gen_add_code_line("int i = feg_job_i[jb]; int vrow = feg_job_vrow[jb];")
        self.gen_add_code_line("#pragma unroll")
        self.gen_add_code_line("for (int r = 0; r < 6; ++r) {", True)
        self.gen_add_code_line("s_dtau_dfext[vrow + " + str(nv) + "*(6*i + r)] += s_feg_slab[6*jb + r];")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow()
        self.gen_add_sync()
    self.gen_add_end_function()


def gen_f_ext_gradient_output_size(self):
    """Number of T elements in EACH of the two first-order f_ext-grad outputs:
    dtau_dfext and dqdd_dfext are nv x (6*NB)."""
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    return nv * 6 * NB


def _f_ext_gradient_dq_smem_count(self):
    """T-element shared count for the -dJ^T/dq kernel arena (fixed base).

    Layout in s_temp: s_qpert[n_pos] | s_JTp[nv*6NB] | s_JTm[nv*6NB] |
    s_jt_temp[jt_inner] | s_xi_scratch[xi]. The XImats buffer + s_q live in their
    own arena regions (declared via gen_XImats_helpers_temp_shared_memory_code)."""
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    n_pos = self.robot.get_num_pos()
    out6 = 6 * NB
    jt_temp = self.gen_f_ext_gradient_inner_temp_mem_size()
    xi_scratch = self.gen_load_update_XImats_helpers_temp_mem_size()
    return n_pos + 2 * nv * out6 + jt_temp + xi_scratch


def _emit_f_ext_gradient_dq_body(self, out_ptr_expr):
    """Emit the per-timestep -dJ^T/dq FD body. Assumes s_q (smem), s_XImats, and
    s_temp arena are already declared/loaded. Writes into `out_ptr_expr` (a global
    or shared pointer to the nv*6NB*nv output for this timestep)."""
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    n_pos = self.robot.get_num_pos()
    out6 = 6 * NB
    jt_temp = self.gen_f_ext_gradient_inner_temp_mem_size()
    self.gen_add_code_line("T *s_did_du_dfext = " + out_ptr_expr + ";")
    self.gen_add_code_line("const T fd_h = static_cast<T>(1e-3);")
    self.gen_add_code_line("T *s_qpert = s_temp;")
    self.gen_add_code_line("T *s_JTp = &s_temp[" + str(n_pos) + "];")
    self.gen_add_code_line("T *s_JTm = &s_temp[" + str(n_pos + nv * out6) + "];")
    self.gen_add_code_line("T *s_jt_temp = &s_temp[" + str(n_pos + 2 * nv * out6) + "];")
    self.gen_add_code_line("T *s_xi_scratch = &s_temp[" + str(n_pos + 2 * nv * out6 + jt_temp) + "];")
    # loop over each q coordinate qi in [0, nv)
    self.gen_add_code_line("for (int qi = 0; qi < " + str(nv) + "; ++qi) {", True)
    # s_qpert = s_q (copy) for +h
    self.gen_add_parallel_loop("ind", str(n_pos))
    self.gen_add_code_line("s_qpert[ind] = s_q[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("s_qpert[qi] += fd_h;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_load_update_XImats_helpers_function_call(updated_var_names={"s_q_name": "s_qpert", "s_temp_name": "s_xi_scratch"})
    self.gen_add_sync()
    self.gen_f_ext_gradient_inner_function_call(updated_var_names={
        "s_dtau_dfext_name": "s_JTp", "s_q_name": "s_qpert", "s_temp_name": "s_jt_temp"})
    self.gen_add_sync()
    # re-copy s_qpert = s_q then -h (the XImats helper may have written into
    # s_xi_scratch only, but be defensive and re-seed s_qpert from s_q).
    self.gen_add_parallel_loop("ind", str(n_pos))
    self.gen_add_code_line("s_qpert[ind] = s_q[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_serial_ops()
    self.gen_add_code_line("s_qpert[qi] -= fd_h;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_load_update_XImats_helpers_function_call(updated_var_names={"s_q_name": "s_qpert", "s_temp_name": "s_xi_scratch"})
    self.gen_add_sync()
    self.gen_f_ext_gradient_inner_function_call(updated_var_names={
        "s_dtau_dfext_name": "s_JTm", "s_q_name": "s_qpert", "s_temp_name": "s_jt_temp"})
    self.gen_add_sync()
    # central diff into output column qi: out[...][qi] = (JTp - JTm)/(2h).
    # s_JTp/s_JTm already hold -J^T (the inner emits -J^T), so this is -dJ^T/dq.
    # output layout: [ (row v_j) + nv*(6NB col) + nv*6NB*qi ]
    self.gen_add_parallel_loop("ind", str(nv * out6))
    self.gen_add_code_line("s_did_du_dfext[ind + " + str(nv * out6) + "*qi] = "
                           "(s_JTp[ind] - s_JTm[ind]) / (static_cast<T>(2)*fd_h);")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_control_flow()  # for qi
    # Restore s_XImats / s_q-state for the ORIGINAL q so any later use is correct.
    self.gen_load_update_XImats_helpers_function_call(updated_var_names={"s_temp_name": "s_xi_scratch"})
    self.gen_add_sync()


def gen_f_ext_gradient_dq_kernel(self, single_call_timing=False):
    """Emit f_ext_gradient_dq_kernel: the mixed second-order block
    d(id_du)/dfext = -dJ^T/dq  (section A.3), size nv x (6*NB) x nv, FIXED-BASE.

    Central finite-difference of the analytic A.1 -J^T over each generalized
    coordinate (the same FD-on-Jacobian approach the d2ee GPU path uses for the
    kinematic Hessian). A velocity-coordinate perturbation equals q[i] += h
    directly on a fixed base. Floating-base A.3 needs the SE(3) Lie integrator to
    perturb along the root twist and is DEFERRED (the numpy + pinocchio oracle
    ships -dJ^T/dq for BOTH modes, so the math is validated). The q-dot block is
    identically zero (J^T is q-only) and is not emitted."""
    if self.robot.floating_base:
        return
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    n_pos = self.robot.get_num_pos()
    out6 = 6 * NB
    out_each = nv * out6 * nv

    func_params = [
        "d_did_du_dfext is the output -dJ^T/dq, size NV*(6*NB)*NV = " + str(out_each) + " per timestep",
        "d_q is the joint positions, stride_q the per-timestep stride",
        "d_robotModel is the initialized model helpers on the GPU",
        "NUM_TIMESTEPS is the trajectory length (or timing reps)",
    ]
    func_def_start = ("void f_ext_gradient_dq_kernel(T *d_did_du_dfext, "
                      "const T *d_q, const int stride_q, ")
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    self.gen_add_func_doc("Compute -dJ^T/dq = d(id_du)/dfext (section A.3, fixed base, batched kernel)",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    shared_extra = _f_ext_gradient_dq_smem_count(self)
    self.gen_XImats_helpers_temp_shared_memory_code(
        shared_extra, extra_t_buffers=[("s_q", n_pos)], include_linalg_scratch=True)
    if not single_call_timing:
        self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
        self.gen_kernel_load_inputs("q", str(n_pos), stride="stride_q")
        self.gen_add_code_line("// compute")
        _emit_f_ext_gradient_dq_body(self, "&d_did_du_dfext[k*" + str(out_each) + "]")
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q", str(n_pos))
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q", str(n_pos), feedback_from="did_du_dfext")
        _emit_f_ext_gradient_dq_body(self, "d_did_du_dfext")
        self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_f_ext_gradient_dq_host(self, mode=0):
    """Host wrapper for the -dJ^T/dq kernel (fixed base only)."""
    if self.robot.floating_base:
        return
    single_call_timing = (mode == 1)
    compute_only = (mode == 2)
    func_params = [
        "hd_data is the packaged input and output pointers",
        "d_robotModel is the initialized model helpers on the GPU",
        "num_timesteps is the trajectory length (or timing reps)",
        "streams are CUDA streams for async transfers",
    ]
    func_def_start = ("void f_ext_gradient_dq(gridData<T, KIND> *hd_data, "
                      "const robotModel<T> *d_robotModel, const int num_timesteps,")
    func_def_end = "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    self.gen_add_func_doc("Compute -dJ^T/dq = d(id_du)/dfext (host wrapper, fixed base)", [], func_params, None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"f_ext_gradient_dq requires all-data or dynamics gridData\");")
    out_each = "NUM_VEL*6*NUM_BODIES*NUM_VEL"
    func_call_start = ("f_ext_gradient_dq_kernel<T><<<block_dimms,thread_dimms,F_EXT_GRAD_DQ_DYNAMIC_SHARED_MEM_BYTES<T>()>>>("
                       "hd_data->d_did_du_dfext,hd_data->d_q,stride_q,")
    func_call_end = "d_robotModel,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>", "kernel_single_timing<T>")
    if not compute_only:
        self.gen_add_code_lines([
            "// start code with memory transfer",
            "int stride_q;",
            "if (USE_COMPRESSED_MEM) {stride_q = NUM_JOINTS; gpuErrchk(cudaMemcpyAsync(hd_data->d_q,hd_data->h_q,stride_q*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}",
            "else {stride_q = 3*NUM_JOINTS; gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}",
            "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q = USE_COMPRESSED_MEM ? NUM_JOINTS: 3*NUM_JOINTS;")
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q", "hd_data->d_q_qd_u") + "}"
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    if single_call_timing:
        func_call_code.insert(0, "struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"f_ext_gradient_dq\", F_EXT_GRAD_DQ_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_did_du_dfext,hd_data->d_did_du_dfext," + out_each + "*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();"])
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("f_ext_gradient_dq"))
    self.gen_add_end_function()


def gen_f_ext_gradient_device(self):
    """Emit f_ext_gradient_device: computes dtau/dfext = -J^T and
    dqdd/dfext = M^{-1} J^T into caller-provided shared buffers.

    Reuses direct_minv_inner for s_Minv (the same inverse-inertia buffer fd_du
    consumes) and the f_ext_gradient_jacobianT_inner for -J^T, then one
    nv x nv * nv x 6NB GEMM (dqdd = -Minv @ dtau). Both outputs are q-only and
    f_ext-VALUE independent (so this device takes q, not f_ext)."""
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    out6 = 6 * NB

    func_params = [
        "s_dtau_dfext is the output dtau/dfext = -J^T, size NV*6*NB = " + str(nv * out6),
        "s_dqdd_dfext is the output dqdd/dfext = M^{-1} J^T, size NV*6*NB = " + str(nv * out6),
        "s_q is the vector of joint positions",
        "d_robotModel is the initialized model helpers on the GPU",
    ]
    func_notes = [
        "Both outputs are q-only (f_ext enters RNEA additively & linearly).",
        "dqdd = -Minv @ dtau (since dtau = -J^T, M^{-1} J^T = -Minv @ dtau).",
    ]
    func_def = ("void f_ext_gradient_device(T *s_dtau_dfext, T *s_dqdd_dfext, "
                "const T *s_q, const robotModel<T> *d_robotModel) {")
    # scratch: max of the J^T inner temp and the direct_minv inner temp, plus an
    # nv*nv s_Minv buffer.
    jt_temp = self.gen_f_ext_gradient_inner_temp_mem_size()
    minv_temp = self.gen_direct_minv_inner_temp_mem_size()
    shared_extra = nv * nv + max(jt_temp, minv_temp)

    self.gen_add_func_doc("Compute the f_ext gradient (dtau/dfext, dqdd/dfext)",
                          func_notes, func_params, None)
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line(func_def, True)
    self.gen_XImats_helpers_temp_shared_memory_code(
        shared_extra, extra_t_buffers=None, include_linalg_scratch=True)
    self.gen_load_update_XImats_helpers_function_call()
    # s_Minv lives at head of s_temp; the inner scratch follows.
    self.gen_add_code_line("T *s_Minv = s_temp;")
    self.gen_add_code_line("T *s_fext_temp = &s_temp[" + str(nv * nv) + "];")
    # build -J^T
    self.gen_f_ext_gradient_inner_function_call(
        updated_var_names={"s_temp_name": "s_fext_temp"})
    self.gen_add_sync()
    # Minv into s_Minv (F kept in smem; inner slices it from the tail of its temp)
    self.gen_direct_minv_inner_function_call(
        updated_var_names={"s_Minv_name": "s_Minv", "s_temp_name": "s_fext_temp"},
        f_in_smem_expr="true")
    self.gen_add_sync()
    # densify Minv upper->full (direct_minv outputs SYMMETRIC_UPPER)
    self.gen_add_code_line("// densify Minv (direct_minv emits symmetric-upper)")
    self.gen_add_parallel_loop("ind", str(nv * nv))
    self.gen_add_code_line("int r = ind % " + str(nv) + "; int c = ind / " + str(nv) + ";")
    self.gen_add_code_line("if (c < r) { s_Minv[r + " + str(nv) + "*c] = s_Minv[c + " + str(nv) + "*r]; }")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    # dqdd = -Minv @ dtau  (both nv x 6NB; Minv is nv x nv symmetric)
    self.gen_add_code_line("// dqdd/dfext = M^{-1} J^T = -Minv @ (dtau/dfext)")
    self.gen_add_parallel_loop("ind", str(nv * out6))
    self.gen_add_code_line("int row = ind % " + str(nv) + "; int col = ind / " + str(nv) + ";")
    # (Minv @ dtau)[row,col] = sum_k Minv[row,k] dtau[k,col]
    self.gen_add_code_line("T acc = static_cast<T>(0);")
    self.gen_add_code_line("for (int k = 0; k < " + str(nv) + "; ++k) {", True)
    self.gen_add_code_line("acc += s_Minv[row + " + str(nv) + "*k] * s_dtau_dfext[k + " + str(nv) + "*col];")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("s_dqdd_dfext[ind] = -acc;")
    self.gen_add_end_control_flow()
    self.gen_add_sync()
    self.gen_add_end_function()


def gen_f_ext_gradient_kernel(self, single_call_timing=False):
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    n_pos = self.robot.get_num_pos()
    out_each = nv * 6 * NB

    func_params = [
        "d_dtau_dfext / d_dqdd_dfext are the two outputs (each NV*6*NB per timestep)",
        "d_q is the joint positions, stride_q the per-timestep stride",
        "d_robotModel is the initialized model helpers on the GPU",
        "NUM_TIMESTEPS is the trajectory length (or timing reps)",
    ]
    func_def_start = ("void f_ext_gradient_kernel(T *d_dtau_dfext, T *d_dqdd_dfext, "
                      "const T *d_q, const int stride_q, ")
    func_def_end = "const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS) {"
    func_def = func_def_start + func_def_end
    if single_call_timing:
        func_def = func_def.replace("(", "_single_timing(")
    self.gen_add_func_doc("Compute the f_ext gradient (batched kernel)",
                          [], func_params, None)
    self.gen_add_code_line("template <typename T, int RESOURCE_TIER = GRID_DEFAULT_RESOURCE_TIER>")
    self.gen_add_code_line("__global__")
    self.gen_add_code_line("__launch_bounds__(tier_max_threads<RESOURCE_TIER>())")
    self.gen_add_code_line(func_def, True)
    jt_temp = self.gen_f_ext_gradient_inner_temp_mem_size()
    minv_temp = self.gen_direct_minv_inner_temp_mem_size()
    shared_extra = nv * nv + max(jt_temp, minv_temp)
    self.gen_XImats_helpers_temp_shared_memory_code(
        shared_extra, extra_t_buffers=[("s_q", n_pos), ("s_dtau_dfext", out_each),
                                       ("s_dqdd_dfext", out_each)],
        include_linalg_scratch=True)

    def _body():
        self.gen_add_code_line("T *s_Minv = s_temp;")
        self.gen_add_code_line("T *s_fext_temp = &s_temp[" + str(nv * nv) + "];")
        self.gen_load_update_XImats_helpers_function_call()
        self.gen_f_ext_gradient_inner_function_call(
            updated_var_names={"s_temp_name": "s_fext_temp"})
        self.gen_add_sync()
        self.gen_direct_minv_inner_function_call(
            updated_var_names={"s_Minv_name": "s_Minv", "s_temp_name": "s_fext_temp"},
            f_in_smem_expr="true")
        self.gen_add_sync()
        self.gen_add_parallel_loop("ind", str(nv * nv))
        self.gen_add_code_line("int r = ind % " + str(nv) + "; int c = ind / " + str(nv) + ";")
        self.gen_add_code_line("if (c < r) { s_Minv[r + " + str(nv) + "*c] = s_Minv[c + " + str(nv) + "*r]; }")
        self.gen_add_end_control_flow()
        self.gen_add_sync()
        self.gen_add_parallel_loop("ind", str(out_each))
        self.gen_add_code_line("int row = ind % " + str(nv) + "; int col = ind / " + str(nv) + ";")
        self.gen_add_code_line("T acc = static_cast<T>(0);")
        self.gen_add_code_line("for (int k = 0; k < " + str(nv) + "; ++k) { acc += s_Minv[row + " + str(nv) + "*k] * s_dtau_dfext[k + " + str(nv) + "*col]; }")
        self.gen_add_code_line("s_dqdd_dfext[ind] = -acc;")
        self.gen_add_end_control_flow()
        self.gen_add_sync()

    if not single_call_timing:
        self.gen_add_parallel_loop("k", "NUM_TIMESTEPS", block_level=True)
        self.gen_kernel_load_inputs("q", str(n_pos), stride="stride_q")
        self.gen_add_code_line("// compute")
        _body()
        self.gen_kernel_save_result("dtau_dfext", str(out_each), stride=str(out_each))
        self.gen_kernel_save_result("dqdd_dfext", str(out_each), stride=str(out_each))
        self.gen_add_end_control_flow()
    else:
        self.gen_kernel_load_inputs("q", str(n_pos))
        self.gen_add_code_line("// compute with NUM_TIMESTEPS as NUM_REPS for timing")
        self.gen_add_code_line("for (int rep = 0; rep < NUM_TIMESTEPS; rep++){", True)
        self.gen_anti_licm_input_reload("q", str(n_pos), feedback_from="dtau_dfext")
        _body()
        self.gen_anti_licm_output_write("dtau_dfext")
        self.gen_add_end_control_flow()
        self.gen_kernel_save_result("dtau_dfext", str(out_each))
        self.gen_kernel_save_result("dqdd_dfext", str(out_each))
    self.gen_add_end_function()


def gen_f_ext_gradient_host(self, mode=0):
    single_call_timing = (mode == 1)
    compute_only = (mode == 2)
    func_params = [
        "hd_data is the packaged input and output pointers",
        "d_robotModel is the initialized model helpers on the GPU",
        "num_timesteps is the trajectory length (or timing reps)",
        "streams are CUDA streams for async transfers",
    ]
    func_def_start = ("void f_ext_gradient(gridData<T, KIND> *hd_data, "
                      "const robotModel<T> *d_robotModel, const int num_timesteps,")
    func_def_end = "                      const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {"
    if single_call_timing:
        func_def_start = func_def_start.replace("(", "_single_timing(")
        func_def_end = "              " + func_def_end
    if compute_only:
        func_def_start = func_def_start.replace("(", "_compute_only(")
        func_def_end = "             " + func_def_end.replace(", cudaStream_t *streams", "")
    self.gen_add_func_doc("Compute the f_ext gradient (host wrapper)", [], func_params, None)
    self.gen_add_code_line("template <typename T, bool USE_COMPRESSED_MEM = false, gridDataKind KIND = GRID_DATA_ALL>")
    self.gen_add_code_line("__host__")
    self.gen_add_code_line(func_def_start)
    self.gen_add_code_line(func_def_end, True)
    self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"f_ext_gradient requires all-data or dynamics gridData\");")
    NB = self.robot.get_num_bodies()
    nv = self.robot.get_num_vel()
    out_each = "NUM_VEL*6*NUM_BODIES"
    func_call_start = ("f_ext_gradient_kernel<T><<<block_dimms,thread_dimms,F_EXT_GRAD_DYNAMIC_SHARED_MEM_BYTES<T>()>>>("
                       "hd_data->d_dtau_dfext,hd_data->d_dqdd_dfext,hd_data->d_q,stride_q,")
    func_call_end = "d_robotModel,num_timesteps);"
    if single_call_timing:
        func_call_start = func_call_start.replace("kernel<T>", "kernel_single_timing<T>")
    if not compute_only:
        self.gen_add_code_lines([
            "// start code with memory transfer",
            "int stride_q;",
            "if (USE_COMPRESSED_MEM) {stride_q = NUM_JOINTS; gpuErrchk(cudaMemcpyAsync(hd_data->d_q,hd_data->h_q,stride_q*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}",
            "else {stride_q = 3*NUM_JOINTS; gpuErrchk(cudaMemcpyAsync(hd_data->d_q_qd_u,hd_data->h_q_qd_u,stride_q*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyHostToDevice,streams[0]));}",
            "gpuErrchkKernel();"])
    else:
        self.gen_add_code_line("int stride_q = USE_COMPRESSED_MEM ? NUM_JOINTS: 3*NUM_JOINTS;")
    self.gen_add_code_line("// then call the kernel")
    func_call = func_call_start + func_call_end
    func_call_mem_adjust = "if (USE_COMPRESSED_MEM) {" + func_call + "}"
    func_call_mem_adjust2 = "else                    {" + func_call.replace("hd_data->d_q", "hd_data->d_q_qd_u") + "}"
    func_call_code = [func_call_mem_adjust, func_call_mem_adjust2, "gpuErrchkKernel();"]
    if single_call_timing:
        func_call_code.insert(0, "struct timespec start, end; clock_gettime(CLOCK_MONOTONIC,&start);")
        func_call_code.append("clock_gettime(CLOCK_MONOTONIC,&end);")
    self.gen_add_code_line("gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"f_ext_gradient\", F_EXT_GRAD_DYNAMIC_SHARED_MEM_BYTES<T>()));")
    self.gen_add_code_lines(func_call_code)
    if not compute_only:
        self.gen_add_code_lines([
            "// finally transfer the result back",
            "gpuErrchk(cudaMemcpy(hd_data->h_dtau_dfext,hd_data->d_dtau_dfext," + out_each + "*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchk(cudaMemcpy(hd_data->h_dqdd_dfext,hd_data->d_dqdd_dfext," + out_each + "*" + ("num_timesteps*" if not single_call_timing else "") + "sizeof(T),cudaMemcpyDeviceToHost));",
            "gpuErrchkKernel();"])
    if single_call_timing:
        from ..algo_registry import single_call_printf_line
        self.gen_add_code_line(single_call_printf_line("f_ext_gradient"))
    self.gen_add_end_function()


def gen_f_ext_gradient(self):
    """Emit the full f_ext-gradient family: J^T inner, device, kernels, hosts.

    A.1 (-J^T) and A.2 (M^-1 J^T) are emitted for ALL base modes. A.3 (-dJ^T/dq,
    the mixed second-order block) is emitted as a composable device function for
    FIXED-BASE robots only (the FD-on-Jacobian needs the SE(3) Lie integrator for
    floating-base tangent perturbations; deferred to backlog). The numpy +
    pinocchio oracle ships A.3 for BOTH base modes."""
    self.gen_f_ext_gradient_jacobianT_inner()
    self.gen_f_ext_gradient_device()
    # A.3 (-dJ^T/dq) GPU device emit: the mixed second-order block. Emitted as a
    # composable device function for FIXED-BASE robots only (the FD-on-Jacobian
    # needs the SE(3) Lie integrator for floating-base tangent perturbations;
    # deferred to backlog). The kernel/host wire it as the third output
    # (s_did_du_dfext, size nv*6NB*nv) ONLY when emitted (fixed base); on a
    # floating base the third output is absent and the kernel keeps the two
    # first-order outputs. The numpy + pinocchio oracle ships A.3 for BOTH modes.
    self.gen_f_ext_gradient_kernel(single_call_timing=False)
    self.gen_f_ext_gradient_kernel(single_call_timing=True)
    self.gen_f_ext_gradient_host(mode=0)
    self.gen_f_ext_gradient_host(mode=1)
    self.gen_f_ext_gradient_host(mode=2)
    # A.3 (-dJ^T/dq): own kernel + host (fixed base only); separate output buffer
    # d_did_du_dfext so the first-order kernel/host stay byte-identical. The
    # _f_ext_grad_dq_emitted gate keys the KERNEL_ATTR_MANIFEST registration: True
    # only when the kernel/macro are actually emitted (fixed base).
    self._f_ext_grad_dq_emitted = not self.robot.floating_base
    if not self.robot.floating_base:
        self.gen_f_ext_gradient_dq_kernel(single_call_timing=False)
        self.gen_f_ext_gradient_dq_kernel(single_call_timing=True)
        self.gen_f_ext_gradient_dq_host(mode=0)
        self.gen_f_ext_gradient_dq_host(mode=1)
        self.gen_f_ext_gradient_dq_host(mode=2)
