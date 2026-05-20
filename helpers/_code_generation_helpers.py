import os

# GRID_NO_LICM_BARRIER=1 suppresses the anti-LICM machinery (volatile input
# reload, rep-stomp, output->input feedback, syncthreads-protected output write)
# emitted by gen_anti_licm_input_reload + gen_anti_licm_output_write. Use this
# when ptxas hangs on floating-base kernels on older toolkits/GPUs (we've seen
# this on sm_86 with stock CUDA). Batch timings (N=16..256) are unaffected;
# single-call timings on _single_timing kernels may LICM-elide and report
# sub-microsecond nonsense values for ABA/FD/MINV/CRBA/ID/EE_POSE.
# Read at call time, not import time, so CLI flags in grid/run.py can set it
# after this module is already imported.
def _no_licm_barrier() -> bool:
    return os.environ.get("GRID_NO_LICM_BARRIER", "0") == "1"


def gen_add_code_line(self, new_code_line, add_indent_after = False):
    self.code_str += self.indent_level * "    " + new_code_line + "\n"
    if add_indent_after:
        self.indent_level += 1

def gen_add_code_lines(self, new_code_lines, add_indent_after = False):
    """Emit a list of code lines.

    Items are emitted in order via `gen_add_code_line`. A literal `True`
    immediately following a line is consumed and treated as that line's
    `add_indent_after` flag (i.e., opens a new indented block). This lets
    callers write `[\"if (cond) {\", True, ...]` inline rather than splitting
    out a separate `gen_add_code_line(line, True)` call.

    The trailing `add_indent_after` parameter still controls whether the
    overall block indents one more level after all lines are emitted.
    """
    i = 0
    while i < len(new_code_lines):
        line = new_code_lines[i]
        opens = (i + 1 < len(new_code_lines)) and (new_code_lines[i + 1] is True)
        if opens:
            self.gen_add_code_line(line, True)
            i += 2
        else:
            self.gen_add_code_line(line)
            i += 1
    if add_indent_after:
        self.indent_level += 1

def gen_add_end_control_flow(self):
    self.indent_level -= 1
    self.gen_add_code_line("}")

def gen_add_end_function(self):
    self.indent_level -= 1
    self.gen_add_code_line("}\n")

def gen_add_func_doc(self, func_desc, notes = [], params = [], return_val = None):
    self.gen_add_code_line("/**")
    self.gen_add_code_line(" * " + func_desc)
    self.gen_add_code_line(" *")
    if len(notes) > 0:
        self.gen_add_code_line(" * Notes:")
        for note in notes:
            self.gen_add_code_line(" *   " + note)
        self.gen_add_code_line(" *")
    for param in params:
        self.gen_add_code_line(" * @param " + param)
    if return_val is not None:
        self.gen_add_code_line(" * @return " + return_val)
    self.gen_add_code_line(" */")

def gen_add_serial_ops(self, use_thread_group = False):
    if use_thread_group:
        self.gen_add_code_line("if(tgrp.thread_rank() == 0){", True)
    else:
        self.gen_add_code_line("if(threadIdx.x == 0 && threadIdx.y == 0){", True)

def gen_add_parallel_loop(self, var_name, max_val, use_thread_group = False, block_level = False):
    if block_level:
        if use_thread_group:
            print("![ERROR]: BLOCK LEVEL THREAD GROUP LOOP NOT IMPLEMENTED YET")
        else:
            code = "for(int " + var_name + " = blockIdx.x + blockIdx.y*gridDim.x; " + \
                        var_name + " < " + max_val + "; " + var_name + " += gridDim.x*gridDim.y){"
    else:
        if use_thread_group:
            code = "for(int " + var_name + " = tgrp.thread_rank(); " + \
                        var_name + " < " + max_val + "; " + var_name + " += tgrp.size()){"
        else:
            code = "for(int " + var_name + " = threadIdx.x + threadIdx.y*blockDim.x; " + \
                        var_name + " < " + max_val + "; " + var_name + " += blockDim.x*blockDim.y){"
    self.gen_add_code_line(code, True)

def gen_static_array_ind_2d(self, col, row, col_stride = 6):
    return col_stride*col + row

def gen_static_array_ind_3d(self, ind, col, row, ind_stride = 36, col_stride = 6):
    return ind_stride*ind + col_stride*col + row

def gen_add_sync(self, use_thread_group = False):
    if use_thread_group:
        self.gen_add_code_line("tgrp.sync();")
    else:
        self.gen_add_code_line("__syncthreads();")

def gen_add_debug_print_code_line(self, print_code_string, use_thread_group = False):
    self.gen_add_sync(use_thread_group)
    self.gen_add_serial_ops(use_thread_group)
    self.gen_add_code_line(print_code_string)
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_add_debug_print_code_lines(self, print_code_string_arr, use_thread_group = False):
    self.gen_add_sync(use_thread_group)
    self.gen_add_serial_ops(use_thread_group)
    for code_string in print_code_string_arr:
        self.gen_add_code_line(code_string)
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_var_in_list(self, var_name, option_list):
    if len(option_list) == 1:
        return "(" + var_name + " == " + option_list[0] + ")"
    else:
        return "(" + " || ".join(["(" + var_name + " == " + option + ")" for option in option_list]) + ")"

def gen_var_not_in_list(self, var_name, option_list):
    if len(option_list) == 1:
        return "(" + var_name + " != " + option_list[0] + ")"
    else:
        return "(" + " && ".join(["(" + var_name + " != " + option + ")" for option in option_list]) + ")"

def gen_add_multi_threaded_select(self, loop_counter, comparator, counts, select_tuples, USE_NON_BRANCH_ALWAYS = False):
    # first find the resulting type and variable name
    dst_code = []
    for (dst_type, dst_var, select_list) in select_tuples:
        if dst_type is None:
            dst_code.append(dst_var)
        elif "|" in dst_type:
            dst_type_parts = dst_type.split("|")
            dst_code.append(dst_type_parts[0] + dst_var + ")" + dst_type_parts[1])
        else:
            dst_code.append(dst_type + " " + dst_var)
    # then if many things to select gen it and branch
    if len(select_tuples) > 1 and not USE_NON_BRANCH_ALWAYS:
        self.gen_add_code_line("// branch to get pointer locations")
        # init pointers outside of select
        self.gen_add_code_line("; ".join(dst_code)  + ";")
        # if / else if / else to select pointers
        n = len(counts)
        code_end = "}"
        for ind in range(n):
            if ind == 0:
                code_start = "     if (" + loop_counter + " " + comparator + " " + counts[ind] + "){ "
            elif ind < n-1:
                code_start = "else if (" + loop_counter + " " + comparator + " " + counts[ind] + "){ "
            else:
                code_start = "else              { "
            code_middle = ""
            for (dst_type, dst_var, select_list) in select_tuples:
                code_middle += dst_var + " = " + select_list[ind] + "; "
            self.gen_add_code_line(code_start + code_middle + code_end)
    # else use a non-branching selector
    else:
        self.gen_add_code_line("// non-branching pointer selector")
        # get the inverse comparator
        n = len(counts)
        inverse_comparator = comparator.replace("<",">") if "<" in comparator else comparator.replace(">","<")
        inverse_comparator = inverse_comparator + "=" if len(inverse_comparator) == 1 else (inverse_comparator[0] if inverse_comparator != "==" else inverse_comparator)
        for tuple_i in range(len(select_tuples)):
            branch_code = ""
            for ind in range(n):
                if ind == 0 or comparator == "==":
                    if comparator == "==" and ind > 0:
                        branch_code += " + "
                    branch_code += "(" + loop_counter + " " + comparator + " " + counts[ind] + ")" 
                elif ind < n-1:
                    branch_code += " + (" + loop_counter + " " + comparator + " " + counts[ind] + " && " + loop_counter + " " + inverse_comparator + " " + counts[ind-1] + ")"
                else:
                    branch_code += " + (" + loop_counter + " " + inverse_comparator + " " + counts[ind-1] + ")"
                branch_code += " * " + select_tuples[tuple_i][2][ind]
            self.gen_add_code_line(dst_code[tuple_i] + " = " + branch_code + ";")

def gen_kernel_load_inputs(self, name, stride, amount, use_thread_group = False, \
                                 name2 = None, stride2 = 1, amount2 = 1, name3 = None, stride3 = 1, amount3 = 1):
    self.gen_add_code_line("// load to shared mem")
    self.gen_add_code_line("const T *d_" + name + "_k = &d_" + name + "[k*" + stride + "];")
    self.gen_add_parallel_loop("ind",amount,use_thread_group)
    self.gen_add_code_line("s_" + name + "[ind] = d_" + name + "_k[ind];")
    self.gen_add_end_control_flow()
    if name2 is not None:
        self.gen_add_code_line("const T *d_" + name2 + "_k = &d_" + name2 + "[k*" + stride2 + "];")
        self.gen_add_parallel_loop("ind",amount2,use_thread_group)
        self.gen_add_code_line("s_" + name2 + "[ind] = d_" + name2 + "_k[ind];")
        self.gen_add_end_control_flow()
    if name3 is not None:
        self.gen_add_code_line("const T *d_" + name3 + "_k = &d_" + name3 + "[k*" + stride3 + "];")
        self.gen_add_parallel_loop("ind",amount3,use_thread_group)
        self.gen_add_code_line("s_" + name3 + "[ind] = d_" + name3 + "_k[ind];")
        self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_kernel_save_result(self, store_to_name, stride, amount, use_thread_group = False, load_from_name = None):
    if load_from_name is None:
        load_from_name = "s_" + store_to_name
    self.gen_add_code_line("// save down to global")
    self.gen_add_code_line("T *d_" + store_to_name + "_k = &d_" + store_to_name + "[k*" + stride + "];")
    self.gen_add_parallel_loop("ind",amount,use_thread_group)
    self.gen_add_code_line("d_" + store_to_name + "_k[ind] = " + load_from_name + "[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_kernel_load_inputs_single_timing(self, name, amount, use_thread_group = False, \
                                               name2 = None, amount2 = 1, name3 = None, amount3 = 1):
    self.gen_add_code_line("// load to shared mem")
    self.gen_add_parallel_loop("ind",amount,use_thread_group)
    self.gen_add_code_line("s_" + name + "[ind] = d_" + name + "[ind];")
    self.gen_add_end_control_flow()
    if name2 is not None:
        self.gen_add_parallel_loop("ind",amount2,use_thread_group)
        self.gen_add_code_line("s_" + name2 + "[ind] = d_" + name2 + "[ind];")
        self.gen_add_end_control_flow()
    if name3 is not None:
        self.gen_add_parallel_loop("ind",amount3,use_thread_group)
        self.gen_add_code_line("s_" + name3 + "[ind] = d_" + name3 + "[ind];")
        self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def gen_anti_licm_input_reload(self, name, amount, use_thread_group = False, \
                                     name2 = None, amount2 = 1, name3 = None, amount3 = 1, \
                                     feedback_from = None):
    """Inside a `for (rep ...)` single_timing loop, reload all inputs from
    device memory via a `const volatile T *` cast, stomp one slot of each
    input with `static_cast<T>(rep)`, and (when `feedback_from` is given)
    inject the previous rep's output back into the input. This creates a
    true loop-carried data dependency that no LICM pass can hoist.

    The feedback chain is the strongest defense: input[N] := f(d_input,
    rep, d_output[(rep-k) & 0x3FF])  with output[N] written by the rep N
    body. Because input N depends on output N-1, output N depends on input
    N, and so on, ptxas would need to symbolically execute every iteration
    to find a fixed point — far beyond any LICM pass's budget.

    Earlier defenses tried in order of failure:
      1. `__noinline__ grid_licm_barrier()` — defeated when ptxas stripped
         the function as no-op self-stores.
      2. Rep-stomp alone (`s_input[rep % N] = static_cast<T>(rep)`) — kept
         in this helper as the first line of defense, but ptxas can still
         prove subsets of work loop-invariant via partial value-range analysis
         (observed for end_effector_pose_gradient on sm_86 / CUDA 12.6).

    Joint-position values are safe across all GRiD algos: sin/cos are
    well-defined for any float; the algos don't assert ranges on q/qd/u.
    `feedback_from` is the name of the output buffer
    (matching `gen_anti_licm_output_write`'s `store_to_name`), so the read
    is `d_<feedback_from>[(rep - k) & 0x3FF]`. The output buffer is sized
    NUM_TIMESTEPS * output_per_step which is comfortably > 1024.
    """
    if _no_licm_barrier():
        # Opt-out path: emit nothing. Inputs were loaded before the rep loop;
        # the rep body runs against stable data and nvcc is free to LICM-elide.
        # Trade-off documented at the top of this module.
        self.gen_add_code_line(
            "// anti-LICM suppressed (GRID_NO_LICM_BARRIER=1); single-call may elide"
        )
        return
    # Both sides are volatile: read forces re-load from global; write to shared
    # via `volatile T *` cast forces nvcc to emit each store and prevents CSE
    # across iterations. The volatile reload alone is insufficient (compiler
    # can prove d_* contents are loop-invariant when no kernel writes them) —
    # the rep-stomp below provides the actual LICM defense.
    self.gen_add_code_line("// anti-LICM: volatile reload of inputs each rep")
    self.gen_add_parallel_loop("_aopt_i", amount, use_thread_group)
    self.gen_add_code_line(
        "reinterpret_cast<volatile T *>(s_" + name + ")[_aopt_i] = "
        "reinterpret_cast<const volatile T *>(d_" + name + ")[_aopt_i];"
    )
    self.gen_add_end_control_flow()
    if name2 is not None:
        self.gen_add_parallel_loop("_aopt_i", amount2, use_thread_group)
        self.gen_add_code_line(
            "reinterpret_cast<volatile T *>(s_" + name2 + ")[_aopt_i] = "
            "reinterpret_cast<const volatile T *>(d_" + name2 + ")[_aopt_i];"
        )
        self.gen_add_end_control_flow()
    if name3 is not None:
        self.gen_add_parallel_loop("_aopt_i", amount3, use_thread_group)
        self.gen_add_code_line(
            "reinterpret_cast<volatile T *>(s_" + name3 + ")[_aopt_i] = "
            "reinterpret_cast<const volatile T *>(d_" + name3 + ")[_aopt_i];"
        )
        self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    # anti-LICM defense (1/2): stomp one input slot with `rep`. First line
    # of defense. The compiler cannot fold the loop induction variable, so
    # at minimum one input slot provably varies per rep.
    self.gen_add_code_line("// anti-LICM (1/2): stomp one input slot with `rep`")
    self.gen_add_code_line("if ((threadIdx.x | threadIdx.y | threadIdx.z) == 0) {")
    self.gen_add_code_line(
        "    reinterpret_cast<volatile T *>(s_" + name + ")[rep % (" + str(amount) + ")] = "
        "static_cast<T>(rep);"
    )
    if name2 is not None:
        self.gen_add_code_line(
            "    reinterpret_cast<volatile T *>(s_" + name2 + ")[rep % (" + str(amount2) + ")] = "
            "static_cast<T>(rep);"
        )
    if name3 is not None:
        self.gen_add_code_line(
            "    reinterpret_cast<volatile T *>(s_" + name3 + ")[rep % (" + str(amount3) + ")] = "
            "static_cast<T>(rep);"
        )
    self.gen_add_code_line("}")
    # anti-LICM defense (2/2): output→input feedback. Reads previous reps'
    # output values from d_<feedback_from> and adds them into input slots.
    # Combined with the output write at end of rep, this creates a closed
    # loop-carried dependency cycle: ptxas would need to symbolically
    # execute all NUM_TIMESTEPS iterations to find any fixed point, far
    # beyond any LICM budget. Reading 3 staggered slots ensures even
    # aggressive cycle analysis can't collapse the chain.
    if feedback_from is not None:
        self.gen_add_code_line(
            "// anti-LICM (2/2): feedback prev rep's d_" + feedback_from
            + " into s_" + name + " (true loop-carried dep)"
        )
        self.gen_add_code_line("if ((threadIdx.x | threadIdx.y | threadIdx.z) == 0) {")
        self.gen_add_code_line(
            "    T _aopt_fb1 = reinterpret_cast<const volatile T *>(d_"
            + feedback_from + ")[(rep + 0x3FF) & 0x3FF];"
        )
        self.gen_add_code_line(
            "    T _aopt_fb2 = reinterpret_cast<const volatile T *>(d_"
            + feedback_from + ")[(rep + 0x3FE) & 0x3FF];"
        )
        self.gen_add_code_line(
            "    T _aopt_fb3 = reinterpret_cast<const volatile T *>(d_"
            + feedback_from + ")[(rep + 0x3FD) & 0x3FF];"
        )
        self.gen_add_code_line(
            "    reinterpret_cast<volatile T *>(s_" + name
            + ")[(rep + 1) % (" + str(amount) + ")] += _aopt_fb1;"
        )
        self.gen_add_code_line(
            "    reinterpret_cast<volatile T *>(s_" + name
            + ")[(rep + 2) % (" + str(amount) + ")] += _aopt_fb2;"
        )
        self.gen_add_code_line(
            "    reinterpret_cast<volatile T *>(s_" + name
            + ")[(rep + 3) % (" + str(amount) + ")] += _aopt_fb3;"
        )
        self.gen_add_code_line("}")
    self.gen_add_sync(use_thread_group)

def gen_anti_licm_output_write(self, store_to_name, load_from_name = None):
    """Inside a `for (rep ...)` single_timing loop, write the per-iter output's
    first element to a varying global address. Companion to
    `gen_anti_licm_input_reload`: the input reload prevents LICM in the SIMT
    path, but cuBLASDx-templated paths can still prove invariance unless
    we also force a per-iter side effect that depends on the iter's work.

    The destination cycles through 1024 slots of d_<store_to_name>, well
    within the MAX_TIMESTEPS=256 × output_size_per_step allocation any
    kernel's output buffer gets in init_gridData.
    """
    if load_from_name is None:
        load_from_name = "s_" + store_to_name
    if _no_licm_barrier():
        self.gen_add_code_line(
            "// anti-LICM output write suppressed (GRID_NO_LICM_BARRIER=1)"
        )
        return
    # __syncthreads() so thread 0 sees other threads' writes to s_<output>.
    # Without this, thread 0 only sees its own writes (or stale init-load values),
    # and if the inner algo doesn't have thread 0 personally write s_<output>[rep & 7],
    # the value is invariant across reps and LICM elides the entire algo.
    self.gen_add_code_line("__syncthreads();")
    self.gen_add_code_line(
        "if ((threadIdx.x | threadIdx.y | threadIdx.z) == 0) { "
        "reinterpret_cast<volatile T *>(d_" + store_to_name + ")[rep & 1023] = "
        "reinterpret_cast<const volatile T *>(" + load_from_name + ")[rep & 7]; }"
    )

def gen_kernel_save_result_single_timing(self, store_to_name, amount, use_thread_group = False, load_from_name = None):
    if load_from_name is None:
        load_from_name = "s_" + store_to_name
    self.gen_add_code_line("// save down to global")
    self.gen_add_parallel_loop("ind",amount,use_thread_group)
    self.gen_add_code_line("d_" + store_to_name + "[ind] = " + load_from_name + "[ind];")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

def _any_algo_uses_workspace_spill(self):
    """Return True when any Phase 2/3 spill-aware algo's PERF pick is non-zero
    (i.e. the kernel writes/reads d_workspace at PERF tier). Used to default
    GRID_CUDA_ENABLE_L2_PERSISTING to 1 — spilled bytes are hot enough to
    benefit from L2 residency. Conservative: also flips on whenever any algo's
    LITE/MINIMAL pick spills (most h1_2 scenarios). False on small robots
    where every algo collapses to no-spill at all tiers."""
    spill_picks = []
    for attr in ("minv_spill_tier_3way", "id_du_spill_tier_3way", "fd_du_spill_tier_3way",
                 "d2ee_spill_tier_3way", "fdsva_so_spill_tier_3way", "fd_spill_tier_3way"):
        picks = getattr(self, attr, None)
        if picks is not None:
            spill_picks.extend(picks)
    # Existing single-pick spill flags (binary; carried from pre-Phase-2a):
    for attr in ("idsva_so_body_frame_use_global_output",
                 "idsva_so_body_frame_grav_full_spill",
                 "idsva_so_world_frame_use_global_output",
                 "fdsva_so_use_workspace_temp",
                 "fdsva_so_fd_grad_use_spill",
                 "d2ee_use_workspace_temp",
                 "d2ee_use_workspace_d2xhom",
                 "id_du_use_global_temp",
                 "fd_du_use_global_temp"):
        if getattr(self, attr, False):
            return True
    return any(p > 0 for p in spill_picks)


def gen_add_shared_memory_helpers(self):
    self.gen_add_code_lines([
        "__host__ __device__ constexpr size_t grid_align_up(size_t offset, size_t alignment) {",
        "    return (offset + alignment - 1) / alignment * alignment;",
        "}",
        "",
        "template <typename U>",
        "__device__ U *grid_arena_ptr(unsigned char *arena, size_t byte_offset) {",
        "    return reinterpret_cast<U *>(arena + byte_offset);",
        "}",
        "",
        "template <typename T>",
        "__host__ __device__ constexpr size_t grid_shared_arena_bytes(size_t t_count, size_t int_count = 0, size_t extra_byte_count = 0) {",
        "    size_t offset = 0;",
        "    offset = grid_align_up(offset, alignof(T));",
        "    offset += sizeof(T) * t_count;",
        "    if (int_count > 0) {",
        "        offset = grid_align_up(offset, alignof(int));",
        "        offset += sizeof(int) * int_count;",
        "    }",
        "    if (extra_byte_count > 0) {",
        "        offset = grid_align_up(offset, static_cast<size_t>(16));",
        "        offset += extra_byte_count;",
        "    }",
        "    return grid_align_up(offset, static_cast<size_t>(16));",
        "}",
        "",
        "#ifndef GRID_CUDA_TARGET_SHARED_MEM_BYTES",
        "#define GRID_CUDA_TARGET_SHARED_MEM_BYTES 98304",
        "#endif",
        "",
        "#ifndef GRID_WORKSPACE_SLOTS",
        "#define GRID_WORKSPACE_SLOTS 1",
        "#endif",
        "",
        "enum gridDataKind { GRID_DATA_ALL = 0, GRID_DATA_DYNAMICS = 1, GRID_DATA_KINEMATICS = 2 };",
        "enum gridSharedTier { GRID_SHARED_FULL = 0, GRID_SPILL_DA_DF_OUTPUT = 1, GRID_SPILL_DV_DA_DF_OUTPUT = 2 };",
        "",
        "#ifndef GRID_CUDA_ENABLE_L2_PERSISTING",
        # Phase 3a/b/c spill design: workspace bytes are HOT (recursion-internal
        # buffers like Minv-F, ABA's interleaved scratch, FDSVA_SO's df_du/Minv)
        # touched many times per kernel. L2 pinning narrows the smem→HBM penalty
        # to smem→L2 (~few-cycle hit instead of 100s of cycles). Default-ON
        # because every codegen target we emit either doesn't spill (= no-op)
        # or spills hot data (= L2 is the right cache). Niche concurrent-kernel
        # workloads can override with -DGRID_CUDA_ENABLE_L2_PERSISTING=0.
        "#define GRID_CUDA_ENABLE_L2_PERSISTING 1",
        "#endif",
        "",
        "__host__ inline cudaError_t grid_get_max_dynamic_shared_memory_bytes(size_t *bytes) {",
        "    int device = 0;",
        "    cudaError_t err = cudaGetDevice(&device);",
        "    if (err != cudaSuccess) { return err; }",
        "    int max_per_block = 0;",
        "    err = cudaDeviceGetAttribute(&max_per_block, cudaDevAttrMaxSharedMemoryPerBlock, device);",
        "    if (err != cudaSuccess) { return err; }",
        "    int max_optin = 0;",
        "#if CUDART_VERSION >= 9000",
        "    err = cudaDeviceGetAttribute(&max_optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, device);",
        "    if (err != cudaSuccess) { cudaGetLastError(); max_optin = 0; }",
        "#endif",
        "    *bytes = static_cast<size_t>(max_optin > max_per_block ? max_optin : max_per_block);",
        "    return cudaSuccess;",
        "}",
        "",
        "__host__ inline cudaError_t grid_check_dynamic_shared_memory_bytes(const char *kernel_name, size_t bytes) {",
        "    size_t max_bytes = 0;",
        "    cudaError_t err = grid_get_max_dynamic_shared_memory_bytes(&max_bytes);",
        "    if (err != cudaSuccess) { return err; }",
        "    if (bytes > max_bytes) {",
        "        fprintf(stderr, \"GRID shared-memory request for %s is %zu bytes, but this device supports %zu bytes per block\\n\",",
        "                kernel_name, bytes, max_bytes);",
        "        return cudaErrorInvalidConfiguration;",
        "    }",
        "    return cudaSuccess;",
        "}",
        "",
        "__host__ inline cudaError_t grid_begin_l2_persisting(cudaStream_t stream, void *ptr, size_t bytes) {",
        "#if GRID_CUDA_ENABLE_L2_PERSISTING && CUDART_VERSION >= 11000",
        "    if (ptr == nullptr || bytes == 0) { return cudaSuccess; }",
        "    int device = 0;",
        "    cudaError_t err = cudaGetDevice(&device);",
        "    if (err != cudaSuccess) { return err; }",
        "    int max_window = 0;",
        "    err = cudaDeviceGetAttribute(&max_window, cudaDevAttrMaxAccessPolicyWindowSize, device);",
        "    if (err != cudaSuccess || max_window <= 0) { cudaGetLastError(); return cudaSuccess; }",
        "    int max_persisting_l2 = 0;",
        "    err = cudaDeviceGetAttribute(&max_persisting_l2, cudaDevAttrMaxPersistingL2CacheSize, device);",
        "    if (err == cudaSuccess && max_persisting_l2 > 0) {",
        "        size_t l2_bytes = bytes < static_cast<size_t>(max_persisting_l2) ? bytes : static_cast<size_t>(max_persisting_l2);",
        "        cudaError_t limit_err = cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize, l2_bytes);",
        "        if (limit_err != cudaSuccess) { cudaGetLastError(); }",
        "    }",
        "    else { cudaGetLastError(); }",
        "    cudaStreamAttrValue attr;",
        "    memset(&attr, 0, sizeof(attr));",
        "    attr.accessPolicyWindow.base_ptr = ptr;",
        "    attr.accessPolicyWindow.num_bytes = bytes < static_cast<size_t>(max_window) ? bytes : static_cast<size_t>(max_window);",
        "    attr.accessPolicyWindow.hitRatio = 0.60;",
        "    attr.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;",
        "    attr.accessPolicyWindow.missProp = cudaAccessPropertyStreaming;",
        "    return cudaStreamSetAttribute(stream, cudaStreamAttributeAccessPolicyWindow, &attr);",
        "#else",
        "    (void)stream; (void)ptr; (void)bytes;",
        "    return cudaSuccess;",
        "#endif",
        "}",
        "",
        "__host__ inline cudaError_t grid_end_l2_persisting(cudaStream_t stream) {",
        "#if GRID_CUDA_ENABLE_L2_PERSISTING && CUDART_VERSION >= 11000",
        "    cudaStreamAttrValue attr;",
        "    memset(&attr, 0, sizeof(attr));",
        "    attr.accessPolicyWindow.num_bytes = 0;",
        "    return cudaStreamSetAttribute(stream, cudaStreamAttributeAccessPolicyWindow, &attr);",
        "#else",
        "    (void)stream;",
        "    return cudaSuccess;",
        "#endif",
        "}",
        ""
    ])

def gen_declare_shared_arena(self, t_buffers, temp_mem_size, include_topology_helpers = True,
                             ximat_name = "s_XImats", ximat_size = 0,
                             temp_name = "s_temp", topology_name = "s_topology_helpers",
                             extra_byte_regions = None,
                             tier_workspace_expr = None):
    """Emit the shared-memory arena layout.

    When ``tier_workspace_expr`` is non-None, the ``s_temp`` slot becomes
    tier-aware: at TIER_PERF the slot is allocated from the arena as usual;
    at TIER_LITE+/MINIMAL the slot is sourced from the supplied workspace
    pointer expression (e.g. ``"s_workspace"``) and the arena allocation
    skips the temp slot entirely, freeing that smem for the caller's outer
    kernel. The arena_offset variable accumulates conditionally so trailing
    slots (topology, linalg, etc.) shift up at LITE+/MINIMAL.

    The caller must declare ``RESOURCE_TIER`` as a template parameter and
    expose ``tier_workspace_expr`` as a function argument.
    """
    if extra_byte_regions is None:
        extra_byte_regions = []
    topology_count = self.gen_topology_helpers_size() if include_topology_helpers else 0
    fixed_t_region_count = sum(int(count) for _, count in t_buffers)
    if ximat_size:
        fixed_t_region_count += int(ximat_size)
    temp_size_int = int(temp_mem_size) if temp_mem_size is not None else 0
    t_region_count = fixed_t_region_count + temp_size_int
    extra_byte_expr = " + ".join(str(count) for _, count in extra_byte_regions) if extra_byte_regions else "0"
    self.gen_add_code_line("// GRID shared arena layout")
    for name, count in t_buffers:
        self.gen_add_code_line("//   T " + name + "[" + str(count) + "]")
    if ximat_size:
        self.gen_add_code_line("//   T " + ximat_name + "[" + str(ximat_size) + "]")
    if temp_size_int != 0:
        if tier_workspace_expr is not None:
            self.gen_add_code_line("//   T " + temp_name + "[" + str(temp_mem_size) + "] (TIER_PERF only; LITE/MINIMAL route to " + tier_workspace_expr + ")")
        else:
            self.gen_add_code_line("//   T " + temp_name + "[" + str(temp_mem_size) + "]")
    if topology_count > 0:
        self.gen_add_code_line("//   int " + topology_name + "[" + str(topology_count) + "]")
    for name, count in extra_byte_regions:
        self.gen_add_code_line("//   bytes " + name + "[" + str(count) + "]")
    self.gen_add_code_line("extern __shared__ __align__(16) unsigned char s_arena[];")
    self.gen_add_code_line("size_t s_arena_offset = 0;")
    for name, count in t_buffers:
        self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
        self.gen_add_code_line("T *" + name + " = grid_arena_ptr<T>(s_arena, s_arena_offset);")
        self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(count) + ");")
    if ximat_size:
        self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
        self.gen_add_code_line("T *" + ximat_name + " = grid_arena_ptr<T>(s_arena, s_arena_offset);")
        self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(ximat_size) + ");")
    if temp_size_int != 0:
        if tier_workspace_expr is not None:
            self.gen_add_code_line("T *" + temp_name + ";")
            self.gen_add_code_line("if constexpr (RESOURCE_TIER == TIER_PERF) {", True)
            self.gen_add_code_line("(void)" + tier_workspace_expr + ";")
            self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
            self.gen_add_code_line(temp_name + " = grid_arena_ptr<T>(s_arena, s_arena_offset);")
            self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(temp_mem_size) + ");")
            self.gen_add_end_control_flow()
            self.gen_add_code_line("else {", True)
            self.gen_add_code_line(temp_name + " = " + tier_workspace_expr + ";")
            self.gen_add_end_control_flow()
        else:
            self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(T));")
            self.gen_add_code_line("T *" + temp_name + " = grid_arena_ptr<T>(s_arena, s_arena_offset);")
            self.gen_add_code_line("s_arena_offset += sizeof(T) * static_cast<size_t>(" + str(temp_mem_size) + ");")
    else:
        self.gen_add_code_line("T *" + temp_name + " = nullptr;")
    if topology_count > 0:
        self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, alignof(int));")
        self.gen_add_code_line("int *" + topology_name + " = grid_arena_ptr<int>(s_arena, s_arena_offset);")
        self.gen_add_code_line("s_arena_offset += sizeof(int) * static_cast<size_t>(" + str(topology_count) + ");")
    for name, count in extra_byte_regions:
        self.gen_add_code_line("unsigned char *" + name + " = nullptr;")
        self.gen_add_code_line("if (static_cast<size_t>(" + str(count) + ") > 0) {", True)
        self.gen_add_code_line("s_arena_offset = grid_align_up(s_arena_offset, static_cast<size_t>(16));")
        self.gen_add_code_line(name + " = grid_arena_ptr<unsigned char>(s_arena, s_arena_offset);")
        self.gen_add_code_line("s_arena_offset += static_cast<size_t>(" + str(count) + ");")
        self.gen_add_end_control_flow()
    self.gen_add_code_line("#ifdef GRID_CUDA_DEBUG_LAYOUT")
    if tier_workspace_expr is not None and temp_size_int != 0:
        # Different t_region_count per tier: PERF includes temp; LITE+ excludes it
        self.gen_add_code_line("if constexpr (RESOURCE_TIER == TIER_PERF) {", True)
        self.gen_add_code_line("assert(s_arena_offset == grid_shared_arena_bytes<T>(" + str(t_region_count) + ", " + str(topology_count) + ", " + extra_byte_expr + "));")
        self.gen_add_end_control_flow()
        self.gen_add_code_line("else {", True)
        self.gen_add_code_line("assert(s_arena_offset == grid_shared_arena_bytes<T>(" + str(fixed_t_region_count) + ", " + str(topology_count) + ", " + extra_byte_expr + "));")
        self.gen_add_end_control_flow()
    else:
        self.gen_add_code_line("assert(s_arena_offset == grid_shared_arena_bytes<T>(" + str(t_region_count) + ", " + str(topology_count) + ", " + extra_byte_expr + "));")
    self.gen_add_code_line("#endif")
    self.gen_add_code_line("(void)s_arena_offset;")

def gen_shared_arena_t_count(self, t_buffers, temp_mem_size, helper_size):
    count = helper_size
    if temp_mem_size is not None:
        count += int(temp_mem_size)
    for _, region_count in t_buffers:
        count += int(region_count)
    return count
