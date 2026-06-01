import os
import numpy as np

class GRiDCodeGenerator:
    # first import helpers to write code generation, spatial algebra, and opology helpers (parent, child, Sind, XImats) and the robotModel object wrapepr
    from .helpers import gen_add_code_line, gen_add_code_lines, gen_add_end_control_flow, gen_add_end_function, \
                         gen_add_func_doc, gen_add_serial_ops, gen_add_parallel_loop, gen_minv_apply, gen_add_sync, gen_var_in_list, \
                         gen_var_not_in_list, gen_add_multi_threaded_select, gen_kernel_load_inputs, gen_kernel_save_result, \
                         gen_anti_licm_input_reload, gen_anti_licm_output_write, \
                         gen_static_array_ind_2d, gen_static_array_ind_3d, gen_add_debug_print_code_lines, \
                         gen_mx_func_call_for_cpp, gen_add_shared_memory_helpers, gen_declare_shared_arena, \
                         gen_shared_arena_t_count, gen_device_wrapper, gen_tier_dispatch, gen_spatial_algebra_helpers, \
                         gen_get_XI_size, gen_init_XImats, gen_load_update_XImats_helpers_temp_mem_size, gen_load_update_XImats_helpers_function_call, \
                         gen_XImats_helpers_temp_shared_memory_code, gen_load_update_XImats_helpers, gen_topology_helpers_size, \
                         gen_get_Xhom_size, gen_load_update_XmatsHom_helpers, gen_load_update_XmatsHom_helpers_function_call, gen_XmatsHom_helpers_temp_shared_memory_code, \
                         gen_topology_sparsity_helpers_python, gen_init_topology_helpers, gen_topology_helpers_pointers_for_cpp, \
                         gen_topology_S_sign_for_cpp, gen_insert_helpers_function_call, gen_insert_helpers_func_def_params, gen_init_robotModel, gen_joint_limits_size, gen_init_joint_limits, \
                         gen_grid_linalg_backend_helpers, gen_linalg_smem_setup, gen_invert_matrix, gen_matmul, gen_matmul_trans, gen_crm_mul, gen_crm, gen_outer_product, custom_is_constant, \
                         robot_has_mimic_joints, _v_slot_cpp, _alpha_for_jid, _alpha_prefix_cpp

    # then import all of the algorithms
    from .algorithms import gen_inverse_dynamics_inner_temp_mem_size, gen_inverse_dynamics_inner_function_call, \
                            gen_inverse_dynamics_device_temp_mem_size, gen_inverse_dynamics_inner, gen_inverse_dynamics_device, \
                            gen_inverse_dynamics_kernel, gen_inverse_dynamics_host, gen_inverse_dynamics, \
                            gen_inverse_dynamics_regressor_inner_temp_mem_size, gen_inverse_dynamics_regressor_inner_function_call, \
                            gen_inverse_dynamics_regressor_inner, gen_inverse_dynamics_regressor_device_temp_mem_size, \
                            gen_inverse_dynamics_regressor_device, gen_inverse_dynamics_regressor_kernel, \
                            gen_inverse_dynamics_regressor_host, gen_inverse_dynamics_regressor, \
                            gen_fd_parameter_gradient_inner_temp_mem_size, gen_fd_parameter_gradient_inner_function_call, \
                            gen_fd_parameter_gradient_inner, gen_fd_parameter_gradient_device_temp_mem_size, \
                            gen_fd_parameter_gradient_device, gen_fd_parameter_gradient_kernel, \
                            gen_fd_parameter_gradient_host, gen_fd_parameter_gradient, \
                            gen_direct_minv_inner_temp_mem_size, gen_direct_minv_inner_F_size, gen_direct_minv_inner_no_F_size, gen_direct_minv_inner_function_call, gen_direct_minv_inner, \
                            gen_direct_minv_device, gen_direct_minv_kernel, gen_direct_minv_host, gen_direct_minv, \
                            gen_forward_dynamics_inner_temp_mem_size, gen_forward_dynamics_inner_F_size, gen_forward_dynamics_finish_function_call, gen_forward_dynamics_finish, \
                            gen_forward_dynamics_inner_function_call, gen_forward_dynamics_inner, gen_forward_dynamics_device, \
                            gen_forward_dynamics_kernel, gen_forward_dynamics_host, gen_forward_dynamics, \
                            gen_inverse_dynamics_gradient_inner_temp_mem_size, gen_inverse_dynamics_gradient_temp_layout, _emit_fb_bfs_level_indexing, \
                            gen_inverse_dynamics_gradient_kernel_max_temp_mem_size, \
                            gen_inverse_dynamics_gradient_inner_function_call, gen_inverse_dynamics_gradient_inner, \
                            gen_inverse_dynamics_gradient_device, gen_inverse_dynamics_gradient_device_function_call, \
                            gen_inverse_dynamics_gradient_kernel, gen_inverse_dynamics_gradient_host, gen_inverse_dynamics_gradient, \
                            gen_forward_dynamics_gradient_inner_temp_mem_size, gen_forward_dynamics_gradient_kernel_max_temp_mem_size, \
                            gen_forward_dynamics_gradient_inner_python, gen_forward_dynamics_gradient_kernel, \
                            gen_forward_dynamics_gradient_device, gen_forward_dynamics_gradient_device_function_call, \
                            gen_forward_dynamics_gradient_host, gen_forward_dynamics_gradient, \
                            gen_f_ext_gradient_inner_temp_mem_size, gen_f_ext_gradient_inner_function_call, \
                            gen_f_ext_gradient_jacobianT_inner, gen_f_ext_gradient_output_size, gen_f_ext_gradient_device, \
                            gen_f_ext_gradient_dq_kernel, gen_f_ext_gradient_dq_host, \
                            gen_f_ext_gradient_kernel, gen_f_ext_gradient_host, gen_f_ext_gradient, \
                            gen_end_effector_pose_inner_temp_mem_size, gen_end_effector_pose_inner_function_call, gen_end_effector_pose_inner, \
                            gen_end_effector_pose_device_temp_mem_size, gen_end_effector_pose_device, gen_end_effector_pose_kernel, \
                            gen_end_effector_pose_host, gen_end_effector_pose_gradient_inner_temp_mem_size, gen_end_effector_pose_gradient_inner_function_call, \
                            gen_end_effector_pose_gradient_inner, gen_end_effector_pose_gradient_device, gen_end_effector_pose_gradient_kernel, \
                            gen_end_effector_pose_gradient_host, gen_end_effector_pose_gradient_hessian_output_count, gen_end_effector_pose_gradient_hessian_inner_temp_mem_size, gen_end_effector_pose_gradient_hessian_inner_function_call, \
                            gen_end_effector_pose_gradient_hessian_inner, gen_end_effector_pose_gradient_hessian_device, gen_end_effector_pose_gradient_hessian_kernel, gen_ee_pose_inner_thread, gen_ee_pose_inner_warp, \
                            gen_ee_pose_inner_xform_from_q_lines, gen_ee_pose_inner_parent_lookup, \
                            gen_ee_pose_fk_batched_kernel, gen_ee_pose_fk_batched_host, \
                            gen_end_effector_pose_gradient_hessian_host, gen_eepose_and_derivatives, \
                            gen_aba, gen_aba_inner, gen_aba_host, \
                            gen_aba_inner_function_call, gen_aba_kernel, gen_aba_device, gen_aba_inner_temp_mem_size, gen_aba_inner_cold_mem_size, \
                            gen_crba, gen_crba_inner_temp_mem_size, gen_crba_inner_function_call, gen_crba_inner, gen_crba_device_temp_mem_size, \
                            gen_crba_device, gen_crba_kernel, gen_crba_host, \
                            gen_idsva_so_xdown_plucker_inverse, gen_idsva_so_reference_order_rt_rp_assembly, \
                            gen_idsva_so_body_frame_inner_temp_mem_size, gen_idsva_so_body_frame_inner_function_call, idsva_so_needs_reference_order_output_repair, \
                            gen_idsva_so_body_frame_reference_order_output_repair, gen_idsva_so_body_frame_floating_reference_inner, gen_idsva_so_body_frame_public_dvdq_layout_repair, gen_idsva_so_body_frame_inner, \
                            gen_idsva_so_body_frame_kernel, gen_idsva_so_body_frame_host, gen_idsva_so_body_frame, \
                            gen_idsva_so_world_frame_temp_mem_size, gen_idsva_so_world_frame_inner, \
                            gen_idsva_so_world_frame_inner_function_call, gen_idsva_so_world_frame_kernel, \
                            gen_idsva_so_world_frame_host, gen_idsva_so_world_frame, \
                            gen_idsva_so_device, gen_idsva_so_dispatcher_host, gen_idsva_so_dispatcher, \
                            gen_floating_gravity_d2tau_dq_temp_mem_size, gen_floating_gravity_d2tau_dq_shared_count, \
                            gen_floating_gravity_d2tau_dq_spill_count, gen_floating_gravity_d2tau_dq_lie_inline, \
                            gen_fdsva_so, gen_fdsva_so_contract_temp_mem_size, gen_fdsva_so_fd_gradient_inline_temp_mem_size, gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled, gen_fdsva_so_fd_gradient_inline, gen_fdsva_so_contract_function_call, gen_fdsva_so_contract, \
                            gen_fdsva_so_device, gen_fdsva_so_device_function_call, gen_fdsva_so_kernel, gen_fdsva_so_host, \
                            gen_integrator_inner_temp_mem_size, gen_integrator_finish_function_call, gen_integrator_finish, \
                            gen_integrator_inner_function_call, gen_integrator_inner, gen_integrator_device, \
                            gen_integrator_kernel, gen_integrator_host, gen_integrator, gen_lie_group_helpers, \
                            gen_integrator_gradient_inner_temp_mem_size, gen_integrator_gradient_dAB_assembly, \
                            gen_integrator_gradient_inner_python, gen_integrator_gradient_multistage, \
                            gen_integrator_gradient_device, gen_integrator_gradient_device_function_call, \
                            gen_integrator_gradient_kernel, gen_integrator_gradient_host, gen_integrator_gradient, \
                            gen_plant_step, gen_plant_step_gradient, gen_quadratic_state_cost, gen_quadratic_input_cost, \
                            gen_ee_pos_cost, gen_plant_barriers, gen_grid_plant, \
                            gen_plant_step_kernel, gen_quadratic_cost_kernel, gen_ee_pos_cost_kernel, gen_plant_kernels, \
                            gen_id_bias_device, gen_id_bias_kernel, gen_id_bias_host, gen_id_bias, \
                            gen_centroidal_inner, gen_com_device, gen_ccrba_device, gen_energy_device, \
                            _gen_kin_centroidal_kernel, _gen_kin_centroidal_host, gen_com, gen_ccrba, gen_energy, \
                            gen_frame_jacobian_inner, gen_frame_jacobian_device, gen_frame_jacobian, \
                            gen_frame_jacobian_dot_device, gen_frame_jacobian_dot, \
                            gen_osc_inertia_device, gen_osc_inertia

    # finally import the test code
    from ._test import test_rnea_fpass, test_rnea_bpass, test_rnea, test_minv_bpass, test_minv_fpass, test_densify_Minv, test_minv, test_rnea_grad_inner, \
                      test_rnea_grad, test_fd_grad, mx0, mx1, mx2, mx3, mx4, mx5, mx, mxS, mxv, fx, fxS, fxv

    # initialize the object
    def __init__(self, robotObj, DEBUG_MODE = False, NEED_PRINT_MAT = False, USE_DYNAMIC_SHARED_MEM = True, FILE_NAMESPACE = "grid"):
        self.robot = robotObj
        self.code_str = ""
        self.indent_level = 0
        self.DEBUG_MODE = DEBUG_MODE
        self.gen_print_mat = DEBUG_MODE or NEED_PRINT_MAT
        # even if dynamic shared mem is not requested for large robots we need to use it
        self.use_dynamic_shared_mem_flag = USE_DYNAMIC_SHARED_MEM or (self.robot.get_num_pos() > 12)
        # check for the file/namespace name
        self.file_namespace = FILE_NAMESPACE
        self.cuda_target_shared_mem_bytes = int(os.environ.get("GRID_CUDA_TARGET_SHARED_MEM_BYTES", "98304"))
        # LITE tier: pick the lowest-spill level whose arena ≤ this target.
        # 48 KB is roughly half the sm_120 ~100 KB per-block cap, so inline-CUDA
        # callers retain ~48 KB for their own outer-kernel scratch.
        self.cuda_target_lite_shared_mem_bytes = int(os.environ.get("GRID_CUDA_TARGET_LITE_SHARED_MEM_BYTES", "49152"))
        self.cuda_shared_mem_type_size_bytes = int(os.environ.get("GRID_CUDA_SHARED_MEM_TYPE_SIZE_BYTES", "4"))

    def _normalize_codegen_algorithms(self, codegen_profile = "all", algorithm_list = None):
        all_algorithms = {
            "id", "minv", "fd", "id_du", "fd_du", "aba", "crba",
            "idsva_so_body_frame", "fdsva_so", "ee_pose", "ee_pose_gradient", "ee_pose_hessian",
            "integrator", "integrator_gradient", "integrator_with_gradient",
            "f_ext_grad", "regressor", "fd_parameter_gradient",
        }
        # E2 (additive, opt-in only): frame_jacobian is NOT part of the default
        # `all` profile so the default-profile header stays byte-identical. It is
        # a recognized key for explicit algorithm_list requests and has its own
        # `frame-jacobian` profile. Requires `ee_pose` (world-transform machinery).
        # E2 CUDA parity (additive, opt-in only): frame_jacobian_dot (Jdot) and
        # osc_inertia (Lambda) join frame_jacobian as recognized-but-not-default
        # keys so the default `all` header stays byte-identical.
        opt_in_algorithms = {"frame_jacobian", "frame_jacobian_dot", "osc_inertia"}
        profile_algorithms = {
            "all": all_algorithms,
            "frame-jacobian": {"ee_pose", "minv", "frame_jacobian",
                               "frame_jacobian_dot", "osc_inertia"},
            "dynamics": {"id", "minv", "fd", "id_du", "fd_du", "aba", "crba", "idsva_so_body_frame", "fdsva_so",
                         "integrator", "integrator_gradient", "integrator_with_gradient"},
            "dynamics-core": {"id", "minv", "fd"},
            "dynamics-gradients": {"id", "minv", "fd", "id_du", "fd_du", "f_ext_grad"},
            "regressor": {"id", "regressor"},
            "fd-param-gradient": {"id", "minv", "fd", "regressor", "fd_parameter_gradient"},
            "f-ext-gradient": {"id", "minv", "f_ext_grad"},
            "kinematics": {"ee_pose"},
            "kinematics-derivatives": {"ee_pose", "ee_pose_gradient", "ee_pose_hessian"},
            "second-order": {"id", "minv", "fd", "id_du", "fd_du", "idsva_so_body_frame", "fdsva_so"},
            "integrators": {"id", "minv", "fd", "id_du", "fd_du", "integrator", "integrator_gradient",
                            "integrator_with_gradient"},
        }
        aliases = {
            "all-dynamics": "dynamics",
            "dynamics-only": "dynamics",
            "kinematics-only": "kinematics",
            "inverse-dynamics": "id",
            "rnea": "id",
            "direct-minv": "minv",
            "forward-dynamics": "fd",
            "inverse-dynamics-gradient": "id_du",
            "id-gradient": "id_du",
            "forward-dynamics-gradient": "fd_du",
            "fd-gradient": "fd_du",
            "f-ext-grad": "f_ext_grad",
            "fext-grad": "f_ext_grad",
            "f-ext-gradient-only": "f_ext_grad",
            "joint-torque-regressor": "regressor",
            "inverse-dynamics-regressor": "regressor",
            "idsva-so": "idsva_so_body_frame",
            "fdsva-so": "fdsva_so",
            "ee-pose": "ee_pose",
            "end-effector-pose": "ee_pose",
            "ee-pose-gradient": "ee_pose_gradient",
            "end-effector-pose-gradient": "ee_pose_gradient",
            "ee-pose-hessian": "ee_pose_hessian",
            "end-effector-pose-hessian": "ee_pose_hessian",
            "euler": "integrator",
            "integrator-euler": "integrator",
            "integrator-gradient": "integrator_gradient",
            "integrator-with-gradient": "integrator_with_gradient",
        }

        def canonicalize(name):
            return aliases.get(str(name).strip().lower().replace("_", "-"), str(name).strip().lower().replace("-", "_"))

        if algorithm_list is None:
            profile_key = str(codegen_profile or "all").strip().lower().replace("_", "-")
            profile_key = aliases.get(profile_key, profile_key)
            if profile_key not in profile_algorithms:
                raise ValueError("Unknown GRiD codegen profile: " + str(codegen_profile))
            algorithms = set(profile_algorithms[profile_key])
        else:
            if isinstance(algorithm_list, str):
                requested = [item for item in algorithm_list.replace(";", ",").split(",") if item.strip()]
            else:
                requested = list(algorithm_list)
            algorithms = set()
            for item in requested:
                key = canonicalize(item)
                if key in profile_algorithms:
                    algorithms.update(profile_algorithms[key])
                elif key in all_algorithms or key in opt_in_algorithms:
                    algorithms.add(key)
                else:
                    raise ValueError("Unknown GRiD algorithm selection: " + str(item))

        # E2: frame_jacobian (+ its Jdot/Lambda siblings) need the world-transform
        # machinery (ee_pose) and, for OSC composition, minv. Pull them in plus the
        # base frame_jacobian emit (the siblings reuse frame_jacobian_inner).
        if "frame_jacobian_dot" in algorithms or "osc_inertia" in algorithms:
            algorithms.add("frame_jacobian")
        if "frame_jacobian" in algorithms:
            algorithms.update({"ee_pose", "minv"})
        if "fd_du" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du"})
        if "id_du" in algorithms:
            algorithms.add("id")
        # f_ext gradient: dtau/dfext reuses the RNEA spatial-transform load (id),
        # dqdd/dfext reuses direct_minv's inner (minv).
        if "f_ext_grad" in algorithms:
            algorithms.update({"id", "minv"})
        if "aba" in algorithms and self.robot.floating_base:
            algorithms.update({"id", "minv", "fd"})
        if "fdsva_so" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du", "fd_du", "idsva_so_body_frame"})
        # Mimic Minv routes through crba_inner (see _direct_minv.py: the reduced-space
        # M is built via CRBA then inverted), so any mimic robot emitting `minv` has a
        # hidden dependency on `crba` for the crba_inner definition. Declare it so the
        # forward-decl'd crba_inner is actually emitted (else nvlink: unresolved extern
        # crba_inner). Non-mimic minv doesn't touch crba, so this is additive — FLAG for
        # main reconcile (shared GCG.py edit; unblocks floating-mimic fdsva_so).
        if "minv" in algorithms and self.robot_has_mimic_joints():
            algorithms.add("crba")
        if "idsva_so_body_frame" in algorithms:
            algorithms.add("id")
            if self.robot.floating_base:
                algorithms.add("id_du")
        # integrator value needs forward dynamics; gradient needs FD + FD-gradient.
        # Floating-base integrator gradients are emitted for all five types
        # (Euler / SI-Euler / Midpoint / RK3 / RK4) and validated against
        # RBDReference. The earlier "forward_dynamics_gradient drops floating
        # linear<->angular velocity coupling" story was a MISDIAGNOSIS: the real
        # defect was a CUDA thread-count race in inverse_dynamics_gradient
        # (missing __syncthreads + a 6-way root accumulation), correct at 32
        # threads and racing above one warp. Fixed; see HANDOFF.md §3.
        if "integrator" in algorithms:
            algorithms.update({"id", "minv", "fd"})
        if "integrator_gradient" in algorithms or "integrator_with_gradient" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du", "fd_du"})
        return algorithms
    
    # add generic code needs and helpers (includes, memory initialization, constants, kernel settings etc.)
    def gen_add_includes(self):
        # first all of the includes
        self.gen_add_code_line("")
        self.gen_add_code_line("#include <assert.h>")
        self.gen_add_code_line("#include <cstddef>")
        self.gen_add_code_line("#include <stdint.h>")
        self.gen_add_code_line("#include <stddef.h>")
        self.gen_add_code_line("#include <math.h>")
        self.gen_add_code_line("#include <stdio.h>")
        self.gen_add_code_line("#include <stdlib.h>")
        self.gen_add_code_line("#include <string.h>")
        self.gen_add_code_line("#include <time.h>")
        self.gen_add_code_line("#include <cuda_runtime.h>")
        self.gen_add_code_lines([
            "",
            "// SIMT GLASS is the only linalg backend (cuBLASDx was removed in v2.0;",
            "// see docs/source/user_guide/concepts/cublasdx_removal_design.rst).",
            "#if defined(__has_include)",
            "#if __has_include(<cub/cub.cuh>)",
            "#define GRID_CUB_HEADER_AVAILABLE 1",
            "#else",
            "#define GRID_CUB_HEADER_AVAILABLE 0",
            "#endif",
            "#else",
            "#define GRID_CUB_HEADER_AVAILABLE 0",
            "#endif",
        ])
        # then any namespaces
        # then any #defines
        self.gen_add_code_lines(["// single kernel timing helper code", \
            "#define time_delta_us_timespec(start,end) (1e6*static_cast<double>(end.tv_sec - start.tv_sec)+1e-3*static_cast<double>(end.tv_nsec - start.tv_nsec))"])
        self.gen_add_code_line("")
        self.gen_add_code_line("#define XIMAT_SIZE 36")

    def gen_add_constants_helpers(self, include_base_inertia = False, include_homogenous_transforms = False):
        # first add constants
        n = self.robot.get_num_pos()
        nv = self.robot.get_num_vel()
        NJ = self.robot.get_num_joints()
        # Dynamics kernels only need the spatial X/I storage. Homogeneous transforms
        # are accounted separately for kinematics kernels.
        XI_size = self.gen_get_XI_size(include_base_inertia,include_homogenous_transforms=False)
        XHom_size, dXhom_size, d2Xhom_size = self.gen_get_Xhom_size()
        dva_cols_per_partial = self.robot.get_total_ancestor_count() + self.robot.get_num_joints()
        max_threads_in_comp_loop = 6*2*dva_cols_per_partial
        max_perf_level_threads = 32 * int(np.ceil(max_threads_in_comp_loop/32.0))
        # cap to 512 mirrors the constant we emit further down; expose on self so
        # _lin_alg_helpers can pin cuBLASDx's BlockDim<TC,1,1> to the same value.
        self.max_perf_level_threads = min(max_perf_level_threads, 512)
        topology_count = self.gen_topology_helpers_size()
        def py_align_up(offset, alignment):
            return ((offset + alignment - 1) // alignment) * alignment

        def py_arena_bytes(t_count, int_count = topology_count):
            offset = 0
            t_align = max(1, min(self.cuda_shared_mem_type_size_bytes, 8))
            offset = py_align_up(offset, t_align)
            offset += self.cuda_shared_mem_type_size_bytes * int(t_count)
            if int_count > 0:
                offset = py_align_up(offset, 4)
                offset += 4 * int(int_count)
            return py_align_up(offset, 16)

        def select_shared_tier_3way(*t_counts):
            """Pick the (perf, lite, minimal) spill-level indices for one algo.
            `t_counts` is the ordered list of arena t_counts at each spill level,
            least-spill first. PERF picks the lowest index whose arena fits
            cuda_target_shared_mem_bytes; LITE picks the lowest index whose
            arena fits cuda_target_lite_shared_mem_bytes (clamped to be ≥ PERF
            pick — LITE can't be less spill than PERF); MINIMAL is always the
            last (most-spill) index."""
            last = len(t_counts) - 1
            perf = next((i for i, t in enumerate(t_counts)
                         if py_arena_bytes(t) <= self.cuda_target_shared_mem_bytes), last)
            lite = next((i for i, t in enumerate(t_counts)
                         if py_arena_bytes(t) <= self.cuda_target_lite_shared_mem_bytes), last)
            lite = max(perf, lite)
            return (perf, lite, last)

        # The ID kernel's s_vaf band is body-indexed (18*NJ); for mimic robots
        # (NJ > n) size it 18*NJ so the inner's body f-writes don't overflow into
        # the XImats region. Non-mimic keeps the legacy 18*n byte-identical.
        _id_vaf = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
        id_t_count = 2*n + n + _id_vaf + n + self.gen_inverse_dynamics_inner_temp_mem_size() + XI_size
        # joint-torque regressor (E1): kernel smem = XI + s_q_qd_qdd(NUM_POS+2nv)
        # + s_Y (nv x 10*NUM_BODIES) + s_vaf(18*NUM_POS) + RNEA forward scratch.
        # n == get_num_pos() here. Additive.
        regressor_t_count = (n + 2*nv) + nv*10*self.robot.get_num_bodies() + 18*n \
            + self.gen_inverse_dynamics_regressor_inner_temp_mem_size() + XI_size
        self.regressor_t_count = regressor_t_count
        # FD param gradient: kernel smem = XI + s_q_qd_u(NUM_POS+2nv) + s_dqdd_dpi
        # + s_Minv(nv*nv) + s_Y(nv x 10*NB) + s_qdd(nv) + s_vaf(18*NUM_POS) + s_c(nv)
        # + the (max) inner forward scratch. n == get_num_pos() here. Additive.
        fd_param_grad_t_count = (n + 2*nv) + nv*10*self.robot.get_num_bodies() \
            + nv*nv + nv*10*self.robot.get_num_bodies() + nv + 18*n + nv \
            + self.gen_fd_parameter_gradient_inner_temp_mem_size() + XI_size
        self.fd_param_grad_t_count = fd_param_grad_t_count
        # FD-param-gradient g1-spill: 2-level surgical ladder. Level 0 keeps every
        # buffer in smem (current behavior on robots that fit). Level 1 spills the
        # s_Y regressor scratch (nv*10*NB, write-once / consumed-once in the final
        # -Minv.Y GEMM) to the L2-pinned d_workspace SO section -- the hot Minv +
        # vaf + inner-RNEA path stays in smem. On g1-floating this drops the arena
        # from ~135 KB to ~94 KB, under the sm_120 ~99 KB cap. The picker selects
        # level 1 for any tier whose level-0 arena overflows the smem target.
        _fpg_Y_count = nv * 10 * self.robot.get_num_bodies()
        _fpg_t_count_full     = fd_param_grad_t_count
        _fpg_t_count_surgical = fd_param_grad_t_count - _fpg_Y_count
        self.fd_param_grad_spill_tier_3way = select_shared_tier_3way(_fpg_t_count_full, _fpg_t_count_surgical)
        self.fd_param_grad_t_count_per_tier = tuple(
            (_fpg_t_count_full, _fpg_t_count_surgical)[i] for i in self.fd_param_grad_spill_tier_3way
        )
        self.fd_param_grad_spill_Y_ws_count = _fpg_Y_count
        # f_ext gradient (section A): kernel smem = XI + s_q + the two nv x (6*NB)
        # outputs + temp (nv*nv s_Minv + max(J^T-inner, direct_minv-inner) scratch).
        _n_pos = self.robot.get_num_pos()
        _NB = self.robot.get_num_bodies()
        _feg_out = nv * 6 * _NB
        _feg_temp = nv*nv + max(self.gen_f_ext_gradient_inner_temp_mem_size(),
                                self.gen_direct_minv_inner_temp_mem_size())
        f_ext_grad_t_count = _n_pos + 2*_feg_out + _feg_temp + XI_size
        # f_ext-gradient (first-order) g1-spill: 2-level surgical ladder. Level 0
        # keeps both outputs (s_dtau_dfext, s_dqdd_dfext) in smem. Level 1 spills
        # s_dqdd_dfext (the SECOND output, written write-once by the final
        # -Minv@s_dtau GEMM) to the L2-pinned d_workspace SO section; s_dtau_dfext
        # (read by that GEMM) + s_Minv + the inner stay in smem. On g1-floating the
        # level-0 arena is ~99.3 KB -- 272 bytes over the sm_120 ~99 KB cap -- so
        # level 1 (~74.6 KB) is what lets it run. Small robots keep level 0.
        _feg_t_count_full     = f_ext_grad_t_count
        _feg_t_count_surgical = f_ext_grad_t_count - _feg_out
        self.f_ext_grad_spill_tier_3way = select_shared_tier_3way(_feg_t_count_full, _feg_t_count_surgical)
        self.f_ext_grad_t_count_per_tier = tuple(
            (_feg_t_count_full, _feg_t_count_surgical)[i] for i in self.f_ext_grad_spill_tier_3way
        )
        self.f_ext_grad_spill_out_ws_count = _feg_out
        # A.3 (-dJ^T/dq) FD kernel (both base modes): arena = XI + s_q + the FD
        # scratch (s_qpert | [floating: s_dv(nv)] | 2x J^T buffers | J^T-inner temp
        # | XImats reload). Floating base adds an nv-sized velocity-perturbation
        # buffer for the SE(3) Lie-group root retract (grid_integrate_floating_q).
        _feg_dq_dv = nv if self.robot.floating_base else 0
        _feg_dq_extra = (_n_pos + _feg_dq_dv + 2*nv*6*_NB
                         + self.gen_f_ext_gradient_inner_temp_mem_size()
                         + self.gen_load_update_XImats_helpers_temp_mem_size())
        f_ext_grad_dq_t_count = _n_pos + _feg_dq_extra + XI_size
        # Minv Phase 3a: per-tier spill picks. Level 0 = F in smem (6*NV*NV
        # bytes); Level 1 = surgical F to L2-pinned workspace.
        _minv_F_count = self.gen_direct_minv_inner_F_size()
        _minv_no_F_count = self.gen_direct_minv_inner_no_F_size()
        _minv_t_count_full     = n + n*n + _minv_F_count + _minv_no_F_count + XI_size
        _minv_t_count_surgical = n + n*n                 + _minv_no_F_count + XI_size
        self.minv_spill_tier_3way = select_shared_tier_3way(_minv_t_count_full, _minv_t_count_surgical)
        self.minv_use_workspace_F = self.minv_spill_tier_3way[0] == 1
        minv_t_count = _minv_t_count_full if not self.minv_use_workspace_F else _minv_t_count_surgical
        self.minv_t_count_per_tier = tuple(
            (_minv_t_count_full, _minv_t_count_surgical)[i] for i in self.minv_spill_tier_3way
        )
        # FD inner-controlled placement: forward_dynamics_inner slices its own
        # Minv-F. The inner-temp size now bundles F (or not) per MINV_F_IN_SMEM,
        # so the full vs surgical kernel arenas come straight from the sized
        # helper (no separate F term — avoids double-counting). Level 0 = F in
        # smem; Level 1 = F in L2-pinned workspace.
        _fd_base = 3*nv + int(self.robot.floating_base) + nv + XI_size
        _fd_t_count_full      = _fd_base + self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=True)
        _fd_t_count_surgical  = _fd_base + self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=False)
        self.fd_spill_tier_3way = select_shared_tier_3way(_fd_t_count_full, _fd_t_count_surgical)
        self.fd_use_workspace_F = self.fd_spill_tier_3way[0] == 1
        fd_t_count = _fd_t_count_full if not self.fd_use_workspace_F else _fd_t_count_surgical
        self.fd_t_count_per_tier = tuple(
            (_fd_t_count_full, _fd_t_count_surgical)[i] for i in self.fd_spill_tier_3way
        )
        # Integrator: kernel-shared t-count layout is
        #   s_q_qd_u (3nv+fb) + s_qdd (nv) + s_stage_qdd ((max_stages-1)*nv)
        #   + s_stage_point ((max_stages-1)*(2nv+fb)) + s_x_kp1 (2nv+fb)
        #   + s_temp (= FD inner)
        # The "+fb" terms account for the floating-base quaternion (q has 1
        # more element than v).  max_stages = 4 (RK4) — see
        # _integrator._max_stages_in_use().
        _max_stages = 4
        _fb = int(self.robot.floating_base)
        _integrator_base = ((3*nv + _fb) + nv
                            + (_max_stages - 1) * nv
                            + (_max_stages - 1) * (2 * nv + _fb)
                            + (2 * nv + _fb) + XI_size)
        integrator_t_count = _integrator_base + self.gen_forward_dynamics_inner_temp_mem_size()
        # Integrator VALUE surgical spill. The dominant inner buffer is the FD
        # inner's Minv F-region (6*NV*NV). Level 0 keeps it in smem; level 1
        # spills ONLY F to d_workspace (the hot FD path stays in smem), mirroring
        # the standalone forward_dynamics kernel's MINV_F_IN_SMEM lever. For
        # h1_2 the value arena overflows by only a few KB, so the surgical F
        # spill is enough — no whole-arena dump.
        _integrator_t_count_full   = _integrator_base + self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=True)
        _integrator_t_count_Fspill = _integrator_base + self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=False)
        self.integrator_spill_tier_3way = select_shared_tier_3way(_integrator_t_count_full, _integrator_t_count_Fspill)
        self.integrator_t_count_per_tier = tuple(
            (_integrator_t_count_full, _integrator_t_count_Fspill)[i] for i in self.integrator_spill_tier_3way)
        # F float-count the value path spills (for grad-section sizing); 0 if no tier spills.
        self.integrator_minv_F_workspace_count = (self.gen_direct_minv_inner_F_size()
                                                  if any(p == 1 for p in self.integrator_spill_tier_3way) else 0)
        # Integrator gradient: kernel-shared t-count layout is
        #   s_q_qd_u (3nv+fb) + s_dAB (2nv*3nv) + s_df_du (nv*2nv) + s_dc_du (nv*2nv) +
        #   s_vaf (18nv) + s_Minv (nv*nv) + s_qdd (nv)
        #   + multi-stage scratch: s_q_orig (nv+fb) + s_qd_orig (nv)
        #     + s_stage_grad_qdd (max_stages*nv) + s_D_qdd_stage (max_stages*nv*3nv)
        #   + s_temp (= FD-grad inner)
        # s_q_orig holds the FULL nq pose (floating-base adds the quaternion slot),
        # so it is nv+fb — must match _emit_body's ("s_q_orig", n+fb) exactly, else
        # the launched dynamic-smem (this t_count) is fb floats short of the arena
        # the kernel slices and the tail buffer overruns shared memory (floating only).
        # max_stages = 4 (RK4) — see _integrator._max_stages_in_use().
        # The multi-stage scratch is always allocated even for single-stage IT;
        # cost is small relative to total (~12*nv² for iiwa14 ≈ 588 floats).
        _max_stages = 4
        # s_vaf is body-indexed (stride 6 over NB bodies). For a MIMIC robot
        # (fixed base) NB > nv, so the composed FD-grad inner writes 18*NB — size
        # the arena's s_vaf term 18*NB to match _emit_body's ("s_vaf", 18*NB)
        # exactly (else the launched dynamic-smem is short and the kernel overruns
        # smem). Non-mimic keeps 18*nv (byte-identical; floating nv > NB).
        _vaf_count = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else nv)
        # +72 for the two 6x6 SE(3) dIntegrate blocks (floating-base gradient;
        # allocated for fixed-base too but unused there).
        integrator_du_t_count = ((3*nv + int(self.robot.floating_base)) + 2*nv*3*nv + 2*(nv*2*nv)
                                 + _vaf_count + nv*nv + nv
                                 + (2*nv + int(self.robot.floating_base)) + _max_stages * nv + _max_stages * nv * 3*nv
                                 + 72
                                 + self.gen_forward_dynamics_gradient_inner_temp_mem_size() + XI_size)
        # The "with x_kp1" variant adds s_x_kp1 (nq+nv = 2nv+fb) on top.
        integrator_du_with_x_kp1_t_count = integrator_du_t_count + 2*nv + int(self.robot.floating_base)
        # Integrator-gradient surgical spill ladder (4 rungs, least-spill first).
        # Each rung spills only cold / output / coalesced matrices to d_workspace,
        # keeping the hot path (s_vaf + the FD-grad scaffold) in smem as long as
        # it fits. The integrator-gradient kernel never runs concurrently with
        # id_du/fd_du/fdsva_so, so its spilled buffers safely reuse those sections.
        #   rung 0: everything in smem.
        #   rung 1: s_D_qdd_stage (max_stages*nv*3nv) -> d_workspace. (g1_fixed)
        #   rung 2: + s_dAB output (2nv*3nv) -> d_workspace, + id_du da_df band
        #           SELECTIVE spill (the FD-grad inner shrinks to the selective
        #           shared count; only the da_df band leaves smem). (g1_floating)
        #   rung 3: + the WHOLE FD-grad inner s_temp -> d_workspace (id_du
        #           global_temp; the inner can't fit a 100 KB box on h1_2). The
        #           gradient scaffold (s_dc_du / s_vaf / s_Minv) stays in smem.
        # PERF picks the lowest fitting rung; MINIMAL is the last (always fits).
        _integrator_du_D_qdd_count = _max_stages * nv * 3 * nv
        _integrator_du_dAB_count = 2 * nv * 3 * nv
        _integrator_du_inner_full = self.gen_forward_dynamics_gradient_inner_temp_mem_size()
        _integrator_du_inner_selective = max(self.gen_direct_minv_inner_temp_mem_size(),
                                             self.gen_inverse_dynamics_gradient_temp_layout()["selective_shared_count"])
        _integrator_du_full = max(integrator_du_t_count, integrator_du_with_x_kp1_t_count)
        _integrator_du_arenas = (
            _integrator_du_full,                                                                                  # 0 full
            _integrator_du_full - _integrator_du_D_qdd_count,                                                     # 1 +Dqdd
            _integrator_du_full - _integrator_du_D_qdd_count - _integrator_du_dAB_count                           # 2 +dAB+selective
                - (_integrator_du_inner_full - _integrator_du_inner_selective),
            _integrator_du_full - _integrator_du_D_qdd_count - _integrator_du_dAB_count - _integrator_du_inner_full,  # 3 +whole inner
        )
        self.integrator_du_spill_tier_3way = select_shared_tier_3way(*_integrator_du_arenas)
        self.integrator_du_t_count_per_tier = tuple(_integrator_du_arenas[i] for i in self.integrator_du_spill_tier_3way)
        _picks = self.integrator_du_spill_tier_3way
        # Placement booleans per rung index: which buffers leave smem.
        self.integrator_du_dqdd_in_smem_per_tier  = tuple(p < 1 for p in _picks)  # spilled at rungs >=1
        self.integrator_du_dab_in_smem_per_tier    = tuple(p < 2 for p in _picks)  # spilled at rungs >=2
        # FD-grad inner level per tier: 0 full smem, 1 selective (da_df band), 2 global_temp (whole inner).
        self.integrator_du_inner_level_per_tier = tuple((0 if p < 2 else (1 if p == 2 else 2)) for p in _picks)
        # d_workspace floats the gradient needs when ANY tier spills: Dqdd + dAB +
        # the whole inner (rung-3 worst case; the rungs reuse the same regions).
        self.integrator_du_workspace_count = (
            (_integrator_du_D_qdd_count + _integrator_du_dAB_count + _integrator_du_inner_full)
            if any(p >= 1 for p in _picks) else 0)
        self.integrator_du_uses_da_df_spill = any(p == 2 for p in _picks)
        self._integrator_du_dqdd_count = _integrator_du_D_qdd_count
        self._integrator_du_dAB_count = _integrator_du_dAB_count
        id_du_temp_layout = self.gen_inverse_dynamics_gradient_temp_layout()
        # Mimic robots emit a DENSE serial id_du inner (no sparse-band spill), so
        # its inner scratch is the dense count; the selective-spill tier doesn't
        # apply (its band offsets are meaningless for the dense layout). Use the
        # dense count for full AND selective so the kernel always sizes smem for
        # the dense buffers and never selects a too-small selective arena.
        _id_du_has_mimic = self.robot_has_mimic_joints()
        id_du_temp_count = self.gen_inverse_dynamics_gradient_inner_temp_mem_size()
        id_du_selective_temp_count = (
            id_du_temp_count if _id_du_has_mimic
            else id_du_temp_layout["selective_shared_count"]
        )
        fd_du_temp_count = self.gen_forward_dynamics_gradient_inner_temp_mem_size()
        fd_du_selective_temp_count = max(self.gen_direct_minv_inner_temp_mem_size(), id_du_selective_temp_count)
        # s_vaf is body-indexed (NB bodies). For a MIMIC robot NB > nv so size
        # 18*NB; non-mimic keeps 18*nv/18*n (byte-identical; floating non-mimic has
        # nv > NB so 18*nv already covers the body writes). The id_device path uses
        # the get_num_pos() (==n) flavour for non-mimic to stay byte-identical with
        # the legacy 18*n; mimic robots route through 18*NB so the inner's
        # body-indexed f writes never overflow s_vaf into the XImats region.
        _vaf_cnt = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else nv)
        _vaf_cnt_id = 18 * (self.robot.get_num_joints() if self.robot_has_mimic_joints() else n)
        id_device_t_count = _vaf_cnt_id + self.gen_inverse_dynamics_inner_temp_mem_size() + XI_size
        minv_device_t_count = self.gen_direct_minv_inner_temp_mem_size() + XI_size
        fd_device_t_count = self.gen_forward_dynamics_inner_temp_mem_size() + XI_size
        id_du_device_t_count = _vaf_cnt + id_du_temp_count + XI_size
        fd_du_device_t_count = (2*nv*nv) + (_vaf_cnt) + nv + (nv*nv) + fd_du_temp_count + XI_size
        id_du_t_count_full = (nv + n) + (2*nv*nv) + (_vaf_cnt) + nv + id_du_temp_count + XI_size
        fd_du_t_count_full = (3*nv + int(self.robot.floating_base)) + (2*nv*nv) + (_vaf_cnt) + nv + (nv*nv) + fd_du_temp_count + XI_size
        id_du_t_count_selective = id_du_t_count_full - id_du_temp_count + id_du_selective_temp_count
        fd_du_t_count_selective = fd_du_t_count_full - fd_du_temp_count + fd_du_selective_temp_count
        id_du_t_count_emergency = id_du_t_count_full - id_du_temp_count
        fd_du_t_count_emergency = fd_du_t_count_full - fd_du_temp_count
        # Per-tier picks (perf, lite, minimal). The existing single-pick flags
        # (id_du_spill_tier etc.) are kept = perf pick so today's emit paths
        # are byte-for-byte unchanged; the lite/minimal indices are exposed
        # only as metadata until the per-tier emit work lands.
        _id_du_arenas = (id_du_t_count_full, id_du_t_count_selective, id_du_t_count_emergency)
        _fd_du_arenas = (fd_du_t_count_full, fd_du_t_count_selective, fd_du_t_count_emergency)
        self.id_du_spill_tier_3way = select_shared_tier_3way(*_id_du_arenas)
        self.fd_du_spill_tier_3way = select_shared_tier_3way(*_fd_du_arenas)
        self.id_du_spill_tier = self.id_du_spill_tier_3way[0]
        self.fd_du_spill_tier = self.fd_du_spill_tier_3way[0]
        self.id_du_use_selective_spill = self.id_du_spill_tier == 1
        self.fd_du_use_selective_spill = self.fd_du_spill_tier == 1
        self.id_du_use_global_temp = self.id_du_spill_tier == 2
        self.fd_du_use_global_temp = self.fd_du_spill_tier == 2
        id_du_t_count = _id_du_arenas[self.id_du_spill_tier]
        fd_du_t_count = _fd_du_arenas[self.fd_du_spill_tier]
        # Per-tier t_counts exposed for tier-aware constexpr metadata.
        self.id_du_t_count_per_tier = tuple(_id_du_arenas[i] for i in self.id_du_spill_tier_3way)
        self.fd_du_t_count_per_tier = tuple(_fd_du_arenas[i] for i in self.fd_du_spill_tier_3way)
        aba_input_t_count = n + 2*nv
        crba_input_t_count = n + nv
        # ABA surgical-spill ladder, 3 rungs. The 140*NJ+138 inner scratch band
        # keeps its hot recursion in smem and spills only the cold sub-band when
        # possible:
        #   level 0 (full)     : whole inner arena in smem (PERF, byte-identical).
        #   level 1 (surgical) : hot band in smem, cold sub-band -> d_cold. The
        #                        smem arena shrinks to the hot region only.
        #   level 2 (workspace): whole inner arena -> L2-pinned workspace (blunt
        #                        MINIMAL fallback).
        _aba_inner_temp_count = self.gen_aba_inner_temp_mem_size()
        _aba_inner_cold_count = self.gen_aba_inner_cold_mem_size()
        # Hot smem arena at the surgical rung: FIXED reclaims the whole 42*n cold
        # tail (hot ends at 98*n); FLOATING reclaims only the 138-float fb* tail
        # above tempVec (the interior vcross slot still relocates to d_cold but
        # cannot be byte-identically compacted out of smem).
        _aba_surgical_inner_count = (_aba_inner_temp_count - 138) if self.robot.floating_base else (98 * NJ)
        _aba_base_count = nv + aba_input_t_count + 12*NJ + XI_size
        _aba_t_count_full      = _aba_base_count + _aba_inner_temp_count
        _aba_t_count_surgical  = _aba_base_count + _aba_surgical_inner_count
        _aba_t_count_workspace = _aba_base_count
        self.aba_spill_tier_3way = select_shared_tier_3way(_aba_t_count_full, _aba_t_count_surgical, _aba_t_count_workspace)
        self.aba_use_workspace_temp = self.aba_spill_tier_3way[0] == 2
        _aba_arenas = (_aba_t_count_full, _aba_t_count_surgical, _aba_t_count_workspace)
        aba_t_count = _aba_arenas[self.aba_spill_tier_3way[0]]
        self.aba_t_count_per_tier = tuple(_aba_arenas[i] for i in self.aba_spill_tier_3way)
        self._aba_inner_cold_count = _aba_inner_cold_count
        # CRBA: the inner scratch band is spilled as one band to L2-pinned
        # workspace at LITE/MINIMAL. Level 0 = scratch in smem (current);
        # Level 1 = scratch redirected to workspace.
        # An intermediate surgical-spill rung (keep hot band in smem, spill only a
        # cold sub-band) was investigated for K-crbarung and DEFERRED: post-I-crba
        # the whole inner band (42*NJ / 36*NJ+slab) is hot with no cold sub-band,
        # and the full arena (<=~19 KB) already fits smem at every default tier.
        # See _crba.py header (gen_crba_inner_temp_mem_size) for the full rationale.
        # Rungs stay at 2 (full | inner-band-to-workspace).
        _crba_base_count = nv*nv + crba_input_t_count + XI_size
        _crba_t_count_full      = _crba_base_count + self.gen_crba_inner_temp_mem_size()
        _crba_t_count_workspace = _crba_base_count
        self.crba_spill_tier_3way = select_shared_tier_3way(_crba_t_count_full, _crba_t_count_workspace)
        self.crba_use_workspace_temp = self.crba_spill_tier_3way[0] == 1
        crba_t_count = _crba_t_count_full if not self.crba_use_workspace_temp else _crba_t_count_workspace
        self.crba_t_count_per_tier = tuple(
            (_crba_t_count_full, _crba_t_count_workspace)[i] for i in self.crba_spill_tier_3way
        )
        ee_t_count = n + 6*self.robot.get_total_leaf_nodes() + self.gen_end_effector_pose_inner_temp_mem_size() + XHom_size
        # Phase 3d (EE_POSE_GRAD): three-tier spill, mirrors D2EE.
        # Level 0 = full smem (inner_temp + s_deePos + dXmatsHom). Level 1 =
        # inner_temp + s_deePos -> workspace/global. Level 2 = also
        # dXmatsHom -> workspace. inner_temp is recursion-hot but L2-pinned at
        # the host wrapper for spill tiers; s_deePos is write-once output.
        _ee_grad_num_ees = self.robot.get_total_leaf_nodes()
        _ee_grad_inner_temp_count = self.gen_end_effector_pose_gradient_inner_temp_mem_size()
        _ee_grad_full_t_count        = n + 6*n*_ee_grad_num_ees + _ee_grad_inner_temp_count + XHom_size + dXhom_size
        _ee_grad_spill_temp_t_count  = n                                                    + XHom_size + dXhom_size
        _ee_grad_spill_dxhom_t_count = n                                                    + XHom_size
        _ee_grad_arenas = (_ee_grad_full_t_count, _ee_grad_spill_temp_t_count, _ee_grad_spill_dxhom_t_count)
        self.ee_grad_spill_tier_3way = select_shared_tier_3way(*_ee_grad_arenas)
        self.ee_grad_spill_tier = self.ee_grad_spill_tier_3way[0]
        self.ee_grad_use_workspace_temp = self.ee_grad_spill_tier >= 1
        self.ee_grad_use_workspace_dxhom = self.ee_grad_spill_tier >= 2
        dee_t_count = _ee_grad_arenas[self.ee_grad_spill_tier]
        self.ee_grad_t_count_per_tier = tuple(_ee_grad_arenas[i] for i in self.ee_grad_spill_tier_3way)
        # D2EE (FD-on-d/dv-Jacobian): two spill levels (the nv^2 output is the only
        # large buffer that can move out of smem). dXhom/d2Xhom are no longer used
        # by the geometric-Jacobian gradient inner the d2ee inner runs internally.
        _d2ee_num_ees = self.robot.get_total_leaf_nodes()
        d2ee_inner_temp_count = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size()
        d2ee_output_count = self.gen_end_effector_pose_gradient_hessian_output_count()
        d2ee_grad_count = 6 * nv * _d2ee_num_ees
        # full smem: q + grad + d2ee_output + inner_temp + Xhom
        d2ee_full_t_count   = n + d2ee_grad_count + d2ee_output_count + d2ee_inner_temp_count + XHom_size
        # output spilled: drop d2ee_output from smem (still need q + grad + inner_temp + Xhom)
        d2ee_spill_t_count  = n + d2ee_grad_count                     + d2ee_inner_temp_count + XHom_size
        _d2ee_arenas = (d2ee_full_t_count, d2ee_spill_t_count, d2ee_spill_t_count)
        if "ee_pose_hessian" in getattr(self, "generated_algorithms", set()):
            self.d2ee_spill_tier_3way = select_shared_tier_3way(*_d2ee_arenas)
        else:
            self.d2ee_spill_tier_3way = (0, 0, 0)
        self.d2ee_spill_tier = self.d2ee_spill_tier_3way[0]
        self.d2ee_use_workspace_output = self.d2ee_spill_tier >= 1
        # Legacy aliases (older surfaces / tests still read these names; both now
        # mean "output spilled to workspace"). d2xhom flag is permanently false:
        # the FD inner never uses d2Xhom.
        self.d2ee_use_workspace_temp = self.d2ee_use_workspace_output
        self.d2ee_use_workspace_d2xhom = False
        d2ee_t_count = _d2ee_arenas[self.d2ee_spill_tier]
        self.d2ee_t_count_per_tier = tuple(_d2ee_arenas[i] for i in self.d2ee_spill_tier_3way)
        # G2 centroidal quick-wins smem t-counts (no tier spill — new, low perf
        # priority families use the full smem arena).
        NB = self.robot.get_num_bodies()
        # generalized_gravity (the larger of the two ID-bias kernels: + s_qd0):
        #   s_q_qd(2n) + s_out(nv) + s_vaf(18n) + s_qd0(nv) + inner_temp(6n) + XI
        self.id_bias_t_count = 2*n + nv + 18*n + nv + 6*n + XI_size
        # com/ccrba/energy share one arena sizing (use the largest input/output):
        #   s_in(<=2n) + s_out(<=6nv+6) + s_A(6nv) + s_com(3) + s_extra(4)
        #   + centroidal_inner_temp + XHom_size
        _centroidal_inner_temp = 16*self.robot.get_num_joints() + 6*nv*NB + 36*NB + 6*nv + 36
        _centroidal_base = 6*nv + 3 + 4 + _centroidal_inner_temp + XHom_size
        self.com_t_count    = n + (3 + 3*nv) + _centroidal_base
        self.ccrba_t_count  = 2*n + (6*nv + 6) + _centroidal_base
        self.energy_t_count = 2*n + 3 + _centroidal_base
        # Size-triggered gravity-shim full-spill. Default OFF; if shim total shared
        # would exceed the target, set self.idsva_so_body_frame_grav_full_spill and
        # let gen_idsva_so_body_frame_inner_temp_mem_size() return the smaller value
        # (the gravity shim's dX/a/da/f/df spill into d_workspace alongside
        # d2X/d2a/d2f). Saves 50-60 KB on g1-class robots without changing
        # iiwa14/go2 behavior.
        def _compute_idsva_body_t_count():
            inner = self.gen_idsva_so_body_frame_inner_temp_mem_size()
            base = (2*nv + n) + inner + XI_size
            full = base + 4*nv**3
            use_global_output = py_arena_bytes(full) > self.cuda_target_shared_mem_bytes
            return base if use_global_output else full, use_global_output

        self.idsva_so_body_frame_grav_full_spill = False
        idsva_so_body_frame_t_count, self.idsva_so_body_frame_use_global_output = _compute_idsva_body_t_count()
        if self.robot.floating_base and py_arena_bytes(idsva_so_body_frame_t_count) > self.cuda_target_shared_mem_bytes:
            self.idsva_so_body_frame_grav_full_spill = True
            idsva_so_body_frame_t_count, self.idsva_so_body_frame_use_global_output = _compute_idsva_body_t_count()
        # Capture the inner temp count for use downstream (world-frame fallback for
        # fixed-base + FDSVA-SO inner sizing) and for the per-tier spill ladder below.
        idsva_so_body_frame_inner_temp_count = self.gen_idsva_so_body_frame_inner_temp_mem_size()

        # ----- idsva_so BODY-frame per-tier spill ladder -----
        # Rungs least->most spill. Flags = (use_global_output, s_temp_in_global, bc_in_global, tp_in_global).
        #   rung0 full:          output + s_temp + BC all in smem
        #   rung1 global_output: 4*NV^3 output tensor -> d_workspace (cheap; coalesced one-shot)
        #   rung2 output_bc:     + BC (36*NB cold buffer, dead before hot loops) -> d_workspace (surgical)
        #   rung3 output_tp:     + ancestor-pair scratch t/p1..p6 (36*len(jids_a), 30-45% of the
        #                          body arena; DEAD through the whole recursion-hot forward sweep,
        #                          live only in the final block-parallel output assembly) ->
        #                          d_workspace. BC stays in smem (slides down to fill the vacated
        #                          t/p span). This surgical rung keeps the entire recursion-hot
        #                          chain in smem and is the highest-payoff cold sub-band on
        #                          humanoid-scale robots.
        #   rung4 output_temp:   + whole s_temp inner arena -> d_workspace (guaranteed-fit fallback)
        # rungs 2 (BC) and 3 (t/p) are mutually exclusive surgical levers (inner enforces it).
        # Fixed-base is the production overflow case. Floating-base BODY is diagnostic
        # (the dispatcher routes floating to the WORLD frame) so it keeps the legacy
        # single-body emit with the gravity shim; its picks are (0,0,0) and unused.
        _idsva_bf_BC = 36 * self.robot.get_num_bodies()
        _idsva_bf_jids_a = len(self.robot.get_jid_ancestor_ids(include_joint=True)[0])
        _idsva_bf_TP = 36 * _idsva_bf_jids_a
        _idsva_bf_base_smem = (2*nv + n) + XI_size                                  # whole s_temp -> global
        _idsva_bf_full     = (2*nv + n) + idsva_so_body_frame_inner_temp_count + XI_size + 4*nv**3
        _idsva_bf_out      = (2*nv + n) + idsva_so_body_frame_inner_temp_count + XI_size
        _idsva_so_body_tiers = [
            ("full",          _idsva_bf_full,                False, False, False, False),
            ("global_output", _idsva_bf_out,                 True,  False, False, False),
            ("output_bc",     _idsva_bf_out - _idsva_bf_BC,  True,  False, True,  False),
            ("output_tp",     _idsva_bf_out - _idsva_bf_TP,  True,  False, False, True),
            ("output_temp",   _idsva_bf_base_smem,           True,  True,  False, False),
        ]
        self._idsva_so_body_tier_table = _idsva_so_body_tiers
        if self.robot.floating_base:
            # Diagnostic path: keep legacy single-body emit + grav shim (no ladder).
            self.idsva_so_body_frame_spill_tier_3way = (0, 0, 0)
            self.idsva_so_body_frame_t_count_per_tier = (idsva_so_body_frame_t_count,) * 3
            self.idsva_so_body_frame_use_ladder = False
        else:
            _idsva_so_body_arenas = tuple(t[1] for t in _idsva_so_body_tiers)
            self.idsva_so_body_frame_spill_tier_3way = select_shared_tier_3way(*_idsva_so_body_arenas)
            self.idsva_so_body_frame_t_count_per_tier = tuple(_idsva_so_body_arenas[i] for i in self.idsva_so_body_frame_spill_tier_3way)
            self.idsva_so_body_frame_use_ladder = True
            # Keep the const flag accurate to the PERF-tier pick.
            self.idsva_so_body_frame_use_global_output = _idsva_so_body_tiers[self.idsva_so_body_frame_spill_tier_3way[0]][2]

        # world-frame path has its own (smaller) scratch — no gravity-shim shared, no
        # main-sweep extras. Sized via gen_idsva_so_world_frame_temp_mem_size.
        idsva_so_world_frame_inner_temp_count = self.gen_idsva_so_world_frame_temp_mem_size() if self.robot.floating_base else idsva_so_body_frame_inner_temp_count
        idsva_so_world_frame_base_t_count = (2*nv + n) + idsva_so_world_frame_inner_temp_count + XI_size
        idsva_so_world_frame_full_t_count = idsva_so_world_frame_base_t_count + 4*nv**3
        # ----- idsva_so WORLD-frame per-tier spill ladder -----
        # Flags = (use_global_output, s_temp_in_global, cold_in_global). The world inner
        # is UN-aliased, so a surgical rung is now landed: the cold trio Xdown (36*NB,
        # dead after Step 3) + v_w/a_w (6*NB each, dead after Step 4's f_w build) can move
        # to d_workspace while the hot arena stays in smem (inner COLD_IN_SMEM=false).
        # Rungs least->most spill:
        #   full:               output + whole s_temp arena in smem
        #   global_output:      4*NV^3 output tensor -> d_idsva_so global (coalesced one-shot)
        #   output_cold:        + surgical cold trio (36*NB + 12*NB) -> d_workspace
        #   output_temp:        + whole s_temp inner arena -> d_workspace (guaranteed-fit fallback)
        _idsva_wf_cold = 36 * self.robot.get_num_bodies() + 12 * self.robot.get_num_bodies()
        _idsva_wf_base_smem = (2*nv + n) + XI_size
        _idsva_so_world_tiers = [
            ("full",          idsva_so_world_frame_full_t_count,                  False, False, False),
            ("global_output", idsva_so_world_frame_base_t_count,                  True,  False, False),
            ("output_cold",   idsva_so_world_frame_base_t_count - _idsva_wf_cold, True,  False, True),
            ("output_temp",   _idsva_wf_base_smem,                                True,  True,  False),
        ]
        self._idsva_so_world_tier_table = _idsva_so_world_tiers
        _idsva_so_world_arenas = tuple(t[1] for t in _idsva_so_world_tiers)
        self.idsva_so_world_frame_spill_tier_3way = select_shared_tier_3way(*_idsva_so_world_arenas)
        self.idsva_so_world_frame_t_count_per_tier = tuple(_idsva_so_world_arenas[i] for i in self.idsva_so_world_frame_spill_tier_3way)
        self.idsva_so_world_frame_use_global_output = _idsva_so_world_tiers[self.idsva_so_world_frame_spill_tier_3way[0]][2]
        idsva_so_world_frame_t_count = self.idsva_so_world_frame_t_count_per_tier[0]
        # d_workspace floats needed per timestep by the idsva_so spill rungs (for so_workspace sizing).
        # Body rungs: 4=output_temp (whole inner arena), 3=output_tp (36*len(jids_a) ancestor-pair
        # scratch), 2=output_bc (36*NB cold slab); 0/1 spill nothing into d_workspace.
        def _idsva_body_ws_floats(pick):
            if pick == 4:
                return idsva_so_body_frame_inner_temp_count
            if pick == 3:
                return _idsva_bf_TP
            if pick == 2:
                return _idsva_bf_BC
            return 0
        def _idsva_world_ws_floats(pick):
            # pick 3 (output_temp) spills the whole inner arena; pick 2 (output_cold)
            # spills just the surgical cold trio (Xdown 36*NB + v_w/a_w 12*NB).
            if pick == 3:
                return idsva_so_world_frame_inner_temp_count
            if pick == 2:
                return 36 * self.robot.get_num_bodies() + 12 * self.robot.get_num_bodies()
            return 0
        idsva_so_spill_ws_t_count = max(
            [_idsva_body_ws_floats(p) for p in self.idsva_so_body_frame_spill_tier_3way] +
            [_idsva_world_ws_floats(p) for p in self.idsva_so_world_frame_spill_tier_3way] + [0])

        # ----- FDSVA_SO shared-mem tier selection -----
        # Four nested tiers, ordered from least-spill to most-spill. Pick the
        # lowest-spill tier whose shared-arena bytes fit cuda_target_shared_mem.
        # Each tier sets three orthogonal state flags read by gen_fdsva_so_*:
        #   - use_global_tensors:  s_idsva_so + s_df2 (8*nv³ outputs) -> d_workspace
        #   - use_workspace_temp:  s_fdsva_temp (4*nv³ inner) -> d_workspace
        #   - fd_grad_use_spill:   fd_grad_inline's da_dq..fxvi band -> d_workspace grad section
        fdsva_so_base_t_count = 4*nv + nv*nv + nv + 2*nv*nv + XI_size
        fdsva_so_contract_temp_count = 4*nv**3
        fdsva_so_fd_gradient_inline_temp_count = self.gen_fdsva_so_fd_gradient_inline_temp_mem_size()
        fdsva_so_fd_gradient_inline_spilled_count = self.gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled()
        # fdsva_so dispatches to world_frame_inner for floating-base (smaller
        # footprint + no grav-shim spill) and body_frame_inner for fixed-base.
        fdsva_so_inner_idsva_so_temp_count = (
            idsva_so_world_frame_inner_temp_count if self.robot.floating_base
            else idsva_so_body_frame_inner_temp_count
        )
        _temp_full     = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_contract_temp_count, fdsva_so_fd_gradient_inline_temp_count)
        _temp_no_contract = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_fd_gradient_inline_temp_count)
        _temp_spilled  = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_fd_gradient_inline_spilled_count)
        # Phase 3e: extend to 6 levels. Each level pushes an additional buffer
        # to L2-pinned workspace. Tuple is
        # (name, shared_count, use_global_tensors, use_workspace_temp,
        #  fd_grad_use_spill, use_workspace_df_du, use_workspace_Minv).
        # Levels 0-3 unchanged from pre-Phase-3e. Level 4 pushes s_df_du
        # (2*NV²); Level 5 also pushes s_Minv (NV²).
        fdsva_so_base_no_df_du = fdsva_so_base_t_count - 2*nv*nv
        fdsva_so_base_no_df_du_no_Minv = fdsva_so_base_no_df_du - nv*nv
        # Level 6: pool -> global. fdsva_so_device runs with SCRATCH_IN_SMEM=false,
        # routing the WHOLE shared s_temp pool (helper sincos + minv + fd + fd_grad +
        # idsva) to d_workspace (reusing the non-concurrent contraction SO-temp region).
        # Smem then holds only the base: inputs + s_qdd + s_Minv + s_df_du + XI
        # (~46-54 KB on h1_2 -> fits the ~99 KB cap). Outputs->device arrays and
        # contraction->global as in levels >=2. Works for BOTH bases because the full
        # inner repoints s_temp and hands the placed pool to the idsva inner (body or
        # world) — the sub-inner just uses the pointer it is given (inner-owns-placement).
        # 8-tuple: (..., use_workspace_idsva_temp == pool->global). Levels 0-5 keep pool in smem.
        _fdsva_so_tiers = [
            ("full",                 fdsva_so_base_t_count + 8*nv**3 + _temp_full,  False, False, False, False, False, False),
            ("global_tensors",       fdsva_so_base_t_count + _temp_full,            True,  False, False, False, False, False),
            ("workspace_temp",       fdsva_so_base_t_count + _temp_no_contract,        True,  True,  False, False, False, False),
            ("workspace_temp_spill", fdsva_so_base_t_count + _temp_spilled,         True,  True,  True,  False, False, False),
            ("spill_df_du",          fdsva_so_base_no_df_du + _temp_spilled,        True,  True,  True,  True,  False, False),
            ("spill_Minv",           fdsva_so_base_no_df_du_no_Minv + _temp_spilled,True,  True,  True,  True,  True,  False),
            # pool->global: smem = base (inputs + qdd + Minv + df_du + XI), no pool/outputs/contraction.
            ("pool_global",          fdsva_so_base_t_count,                         True,  True,  False, False, False, True),
        ]
        _fdsva_so_arenas = tuple(t[1] for t in _fdsva_so_tiers)
        self.fdsva_so_spill_tier_3way = select_shared_tier_3way(*_fdsva_so_arenas)
        _chosen = _fdsva_so_tiers[self.fdsva_so_spill_tier_3way[0]]
        (_, fdsva_so_t_count, self.fdsva_so_use_global_tensors,
         self.fdsva_so_use_workspace_temp, self.fdsva_so_fd_grad_use_spill,
         self.fdsva_so_use_workspace_df_du, self.fdsva_so_use_workspace_Minv,
         self.fdsva_so_use_workspace_idsva_temp) = _chosen
        self.fdsva_so_t_count_per_tier = tuple(_fdsva_so_arenas[i] for i in self.fdsva_so_spill_tier_3way)
        # Phase 3a: include Minv-F count if Minv is spilling (collisions are OK
        # because Minv runs before id_du_grad / fd_grad in any kernel that
        # composes both — they sequentially reuse the same workspace bytes).
        _minv_F_workspace_count = self.gen_direct_minv_inner_F_size() if any(p == 1 for p in self.minv_spill_tier_3way) else 0
        # CRBA whole-arena spill: when crba_inner's scratch band is redirected to
        # d_workspace (LITE/MINIMAL, or a forced deep-spill tier), the per-timestep
        # workspace must be able to back the full 140*NJ-class band. Include it in
        # the grad-section max so the allocation always covers it regardless of the
        # tier the kernel template is instantiated with.
        _crba_inner_temp_count = self.gen_crba_inner_temp_mem_size()
        grad_spill_workspace_t_count = max(id_du_temp_layout["spill_count"],
                                           id_du_temp_count,
                                           fd_du_temp_count,
                                           2*nv*nv,
                                           _crba_inner_temp_count,
                                           _minv_F_workspace_count,
                                           self.integrator_minv_F_workspace_count,
                                           self.integrator_du_workspace_count)
        # D2EE no longer needs d_workspace: under the FD-on-Jacobian inner the
        # ONLY large buffer is the nv^2 output, and when it's spilled the inner
        # writes directly into d_d2eePos (the persistent output buffer) -- not
        # into a per-timestep workspace slice. The new spill thus costs no extra
        # workspace allocation.
        d2ee_workspace_t_count = 0
        # Phase 3d: max workspace required by EE_POSE_GRAD across any tier (PERF
        # may pick 0, but the workspace allocation must cover what LITE/MINIMAL
        # need at runtime when the user switches tier via the kernel template).
        ee_grad_workspace_t_count = 0
        if any(p >= 1 for p in self.ee_grad_spill_tier_3way):
            ee_grad_workspace_t_count += _ee_grad_inner_temp_count
        if any(p >= 2 for p in self.ee_grad_spill_tier_3way):
            ee_grad_workspace_t_count += dXhom_size
        # Include the floating-base gravity-shim spill (Phase D): the d2X/d2a/d2f
        # tensors live in d_workspace instead of shared memory for larger robots.
        idsva_so_body_frame_grav_spill_t_count = self.gen_floating_gravity_d2tau_dq_spill_count() if self.robot.floating_base else 0
        # g1-spill: fd_parameter_gradient (s_Y) and f_ext_gradient (s_dqdd_dfext)
        # surgically spill into this same SO workspace section when their tier picks
        # level >= 1. They never run concurrently with the SO/d2ee/ee_grad kernels,
        # so reuse is safe and costs no new allocation. Fold their spill counts into
        # the max so the per-timestep workspace always covers them.
        _fpg_spill_ws = self.fd_param_grad_spill_Y_ws_count if any(p >= 1 for p in self.fd_param_grad_spill_tier_3way) else 0
        _feg_spill_ws = self.f_ext_grad_spill_out_ws_count if any(p >= 1 for p in self.f_ext_grad_spill_tier_3way) else 0
        so_workspace_t_count = max(8*max(nv**3, 1), d2ee_workspace_t_count, ee_grad_workspace_t_count, idsva_so_body_frame_grav_spill_t_count, idsva_so_spill_ws_t_count, _fpg_spill_ws, _feg_spill_ws)
        # Deprecated launch-count constants remain for external callers that still
        # pass COUNT*sizeof(T).  Make them conservative aliases for the byte arena
        # layouts so those callers do not under-allocate int topology helpers or
        # 16-byte alignment padding.
        legacy_count_pad = topology_count + 8
        legacy_arena_count = lambda t_count: int(t_count + legacy_count_pad)
        _b = lambda flag: "true" if flag else "false"
        # GRID_LINALG_NVIDIA_MAX_HELPER_BYTES is a stub returning 0 (v2.0+).
        # Forward-declared here so the *_DYNAMIC_SHARED_MEM_BYTES helpers
        # below compile before _lin_alg_helpers emits the definition.
        self.gen_add_code_line("template <typename T> __host__ __device__ constexpr size_t GRID_LINALG_NVIDIA_MAX_HELPER_BYTES();")
        self.gen_add_code_lines(["const int NUM_JOINTS = " + str(self.robot.get_num_pos()) + ";", \
                                 "const int NUM_POS = " + str(self.robot.get_num_pos()) + ";", \
                                 "const int NUM_VEL = " + str(self.robot.get_num_vel()) + ";", \
                                 "const int NUM_BODIES = " + str(self.robot.get_num_bodies()) + ";", \
                                 "const int SECOND_ORDER_COORDS = " + str(self.robot.get_num_vel()) + ";", \
                                 "const int SECOND_ORDER_TENSOR_SIZE = " + str(4 * self.robot.get_num_vel()**3) + ";", \
                                 "const int Q_QD_U_STRIDE = " + str(self.robot.get_num_pos() + 2 * self.robot.get_num_vel()) + ";", \
                                 "const int NUM_EES = " + str(self.robot.get_total_leaf_nodes()) + ";", \
                                 "const int TOPOLOGY_HELPERS_COUNT = " + str(topology_count) + ";", \
                                 "const int DYNAMICS_XI_T_COUNT = " + str(XI_size) + ";", \
                                 "const int XHOM_T_COUNT = " + str(XHom_size) + ";", \
                                 "const int DXHOM_T_COUNT = " + str(dXhom_size) + ";", \
                                 "const int D2XHOM_T_COUNT = " + str(d2Xhom_size) + ";", \
                                 "const int GRID_ID_DU_USES_GLOBAL_TEMP = " + str(int(self.id_du_use_global_temp)) + ";", \
                                 "const int GRID_ID_DU_USES_WORKSPACE_ANY_TIER = " + str(1 if any(p >= 1 for p in self.id_du_spill_tier_3way) else 0) + ";", \
                                 "const int GRID_FD_DU_USES_GLOBAL_TEMP = " + str(int(self.fd_du_use_global_temp)) + ";", \
                                 "const int GRID_FD_DU_USES_WORKSPACE_ANY_TIER = " + str(1 if any(p >= 1 for p in self.fd_du_spill_tier_3way) else 0) + ";", \
                                 "const int GRID_ID_DU_USES_DA_DF_SPILL = " + str(int(self.id_du_use_selective_spill)) + ";", \
                                 "const int GRID_FD_DU_USES_DA_DF_SPILL = " + str(int(self.fd_du_use_selective_spill)) + ";", \
                                 "const int GRID_INTEGRATOR_USES_WORKSPACE = " + str(int(any(p == 1 for p in self.integrator_spill_tier_3way))) + ";", \
                                 "const int GRID_INTEGRATOR_DU_USES_WORKSPACE = " + str(int(any(p >= 1 for p in self.integrator_du_spill_tier_3way))) + ";", \
                                 "const int GRID_INTEGRATOR_DU_USES_DA_DF_SPILL = " + str(int(self.integrator_du_uses_da_df_spill)) + ";", \
                                 "const int GRID_GENERATES_IDSVA_SO_BODY_FRAME = " + str(int(getattr(self, "generate_idsva_so_body_frame", True))) + ";", \
                                 "const int GRID_GENERATES_FDSVA_SO = " + str(int(getattr(self, "generate_fdsva_so", True))) + ";", \
                                 "const int GRID_GENERATES_D2EE = " + str(int(getattr(self, "generate_ee_pose_hessian", True))) + ";", \
                                 "const int GRID_IDSVA_SO_USES_GLOBAL_OUTPUT = " + str(int(self.idsva_so_body_frame_use_global_output)) + ";", \
                                 "const int GRID_FDSVA_SO_USES_GLOBAL_TENSORS = " + str(int(self.fdsva_so_use_global_tensors)) + ";", \
                                 # Single-bool kept for inline-CUDA back-compat (reflects PERF-pick only).
                                 "const int GRID_FDSVA_SO_USES_WORKSPACE_TEMP = " + str(int(self.fdsva_so_use_workspace_temp)) + ";", \
                                 # Per-tier-aware gate (true if ANY of PERF/LITE/MINIMAL spills any band
                                 # of the s_temp pool to d_workspace; picks >= 2 cover spill_temp,
                                 # spill_fd_grad_band, spill_df_du, spill_Minv, pool_global). The host
                                 # uses this to pin L2 persistence on d_workspace for the per-tier path.
                                 "const int GRID_FDSVA_SO_USES_WORKSPACE_ANY_TIER = " + str(1 if any(p >= 2 for p in self.fdsva_so_spill_tier_3way) else 0) + ";", \
                                 # GRID_D2EE_USES_WORKSPACE_TEMP: 1 if the PERF tier spills the d2ee output
                                 # (the only large buffer in the new FD-on-Jacobian path) to global memory.
                                 # When spilled, the inner writes directly into d_d2eePos (the persistent
                                 # output buffer) -- no extra per-timestep workspace slice is used. The old
                                 # d2xhom-spill bit is permanently 0 (the FD inner never touches d2Xhom).
                                 "const int GRID_D2EE_USES_WORKSPACE_TEMP = " + str(int(self.d2ee_use_workspace_output)) + ";", \
                                 "const int GRID_D2EE_USES_WORKSPACE_D2XHOM = 0;", \
                                 "const int GRID_D2EE_USES_WORKSPACE_TEMP_ANY = " + str(1 if any(p >= 1 for p in self.d2ee_spill_tier_3way) else 0) + ";", \
                                 "const int GRID_D2EE_SHARED_TIER_VALUE = " + str(self.d2ee_spill_tier) + ";", \
                                 "const int GRID_EE_GRAD_USES_WORKSPACE_TEMP = " + str(int(self.ee_grad_use_workspace_temp)) + ";", \
                                 # Same per-tier gate for the EE_POSE_GRAD chain workspace.
                                 "const int GRID_EE_GRAD_USES_WORKSPACE_TEMP_ANY = " + str(1 if any(p >= 1 for p in self.ee_grad_spill_tier_3way) else 0) + ";", \
                                 "const int GRID_EE_GRAD_USES_WORKSPACE_DXHOM = " + str(int(self.ee_grad_use_workspace_dxhom)) + ";", \
                                 "const int GRID_EE_GRAD_SHARED_TIER_VALUE = " + str(self.ee_grad_spill_tier) + ";", \
                                 "const int GRID_ID_DU_SHARED_TIER_VALUE = " + str(self.id_du_spill_tier) + ";", \
                                 "const int GRID_FD_DU_SHARED_TIER_VALUE = " + str(self.fd_du_spill_tier) + ";", \
                                 "const int ID_DU_TEMP_SPILL_START = " + str(id_du_temp_layout["spill_start"]) + ";", \
                                 "const int ID_DU_TEMP_SPILL_END = " + str(id_du_temp_layout["spill_end"]) + ";", \
                                 "const int ID_DU_TEMP_SPILL_COUNT = " + str(id_du_temp_layout["spill_count"]) + ";", \
                                 "const int ID_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(id_t_count)) + ";", \
                                 "const int MINV_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(minv_t_count)) + ";", \
                                 "const int FD_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(fd_t_count)) + ";", \
                                 "const int ID_DU_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(id_du_t_count)) + ";", \
                                 "const int FD_DU_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(fd_du_t_count)) + ";", \
                                 "const int INTEGRATOR_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(integrator_t_count)) + ";", \
                                 "const int INTEGRATOR_DU_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(max(integrator_du_t_count, integrator_du_with_x_kp1_t_count))) + ";", \
                                 "const int ABA_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(aba_t_count)) + ";", \
                                 "const int CRBA_SHARED_MEM_COUNT = " + str(legacy_arena_count(crba_t_count)) + ";", \
                                 "const int ID_DU_MAX_SHARED_MEM_COUNT = " + str(legacy_arena_count(id_du_t_count_full)) + ";", \
                                 "const int FD_DU_MAX_SHARED_MEM_COUNT = " + str(legacy_arena_count(fd_du_t_count_full)) + ";", \
                                 "const int EE_POS_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(ee_t_count)) + ";", \
                                 "const int DEE_POS_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(dee_t_count)) + ";", \
                                 "const int D2EE_POS_DYNAMIC_SHARED_MEM_COUNT = " + str(legacy_arena_count(d2ee_t_count)) + ";", \
                                 f"const int IDSVA_SO_DYNAMIC_SHARED_MEM_COUNT = {legacy_arena_count(idsva_so_body_frame_t_count)};", \
                                 f"const int FDSVA_SO_DYNAMIC_SHARED_MEM_COUNT = {legacy_arena_count(fdsva_so_t_count)};", \
                                 "const int MAX_PERF_LEVEL_THREADS = " + str(self.max_perf_level_threads) + ";", \
                                 "",
                                 "// Resource-tier API (v2.0): each emitted kernel/_device/_inner takes a",
                                 "// `RESOURCE_TIER` template parameter that picks the (launch_bounds, smem,",
                                 "// register-footprint) profile. TIER_SHARED is the default and is the",
                                 "// current-best perf; TIER_LITE keeps the same launch_bounds but reduces",
                                 "// smem footprint (some intermediates moved to workspace global mem);",
                                 "// TIER_MINIMAL drops launch_bounds to 1024 for maximum block-size flexibility",
                                 "// at the cost of register slack. Inline-CUDA power users with tight outer",
                                 "// kernels pick LITE/MINIMAL to fit GRiD primitives in their resource budget.",
                                 "constexpr int TIER_SHARED    = 0;",
                                 "constexpr int TIER_LITE    = 1;",
                                 "constexpr int TIER_MINIMAL = 2;",
                                 "// Deprecated alias: TIER_PERF was renamed to TIER_SHARED (the tier keeps the",
                                 "// full inner scratch in *shared* memory). Kept so external callsites and",
                                 "// sibling branches that still spell it TIER_PERF keep compiling; prefer",
                                 "// TIER_SHARED in new code.",
                                 "constexpr int TIER_PERF    = TIER_SHARED;",
                                 "",
                                 "// Compile-time override for the default RESOURCE_TIER baked into every",
                                 "// emitted kernel template. Defaults to TIER_SHARED so existing callsites",
                                 "// (kernel<T>, host wrappers that call kernel<T><<<...>>>) keep their",
                                 "// current best-perf semantics. Bench harness sets this via",
                                 "// -DGRID_DEFAULT_RESOURCE_TIER=TIER_LITE (or TIER_MINIMAL) to sweep",
                                 "// per-tier perf without modifying host-wrapper template signatures.",
                                 "#ifndef GRID_DEFAULT_RESOURCE_TIER",
                                 "#define GRID_DEFAULT_RESOURCE_TIER TIER_SHARED",
                                 "#endif",
                                 "",
                                 "// Per-tier launch_bounds upper-bound (= max threads per block nvcc must",
                                 "// budget registers for). sm_120 has 65536 regs/block; nvcc enforces",
                                 "// regs_per_thread * max_threads <= regs_per_block, so a larger max_threads",
                                 "// directly caps regs_per_thread. PERF=SUGGESTED keeps current best perf;",
                                 "// LITE=min(2*SUGGESTED, 768) gives ~85 regs/thread cap (mid-budget);",
                                 "// MINIMAL=1024 gives ~64 regs/thread cap (maximum block-size flexibility).",
                                 "template <int TIER> __host__ __device__ constexpr int tier_max_threads() {",
                                 "    return (TIER == TIER_MINIMAL) ? 1024",
                                 "         : (TIER == TIER_LITE)    ? ((MAX_PERF_LEVEL_THREADS * 2 < 768) ? MAX_PERF_LEVEL_THREADS * 2 : 768)",
                                 "         :                          MAX_PERF_LEVEL_THREADS;",
                                 "}"])
        self.gen_add_code_lines([
                                 "#define GRID_GENERATED_NUM_JOINTS " + str(n),
                                 "#define GRID_GENERATED_NUM_EES " + str(self.robot.get_total_leaf_nodes()),
                                 ""])
        self.gen_add_code_lines([
                                 "template <typename T> __host__ __device__ inline size_t ID_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t INVERSE_DYNAMICS_REGRESSOR_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(regressor_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 # g1-spill: tier-aware. At a spilled tier the s_Y regressor
                                 # scratch moves to d_workspace, shrinking the smem arena. Default
                                 # TIER = TIER_SHARED keeps every existing single-arg call site working.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t FD_PARAMETER_GRADIENT_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.fd_param_grad_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.fd_param_grad_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.fd_param_grad_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 # g1-spill: per-tier placement of s_Y -- true => smem, false => d_workspace.
                                 "template <int TIER> __host__ __device__ constexpr bool FD_PARAMETER_GRADIENT_Y_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.fd_param_grad_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.fd_param_grad_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.fd_param_grad_spill_tier_3way[2] == 0 else "false") + "; }",
                                 # g1-spill: tier-aware. At a spilled tier s_dqdd_dfext (2nd output) moves to d_workspace.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t F_EXT_GRAD_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.f_ext_grad_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.f_ext_grad_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.f_ext_grad_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 # g1-spill: per-tier placement of s_dqdd_dfext -- true => smem, false => d_workspace.
                                 "template <int TIER> __host__ __device__ constexpr bool F_EXT_GRAD_DQDD_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.f_ext_grad_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.f_ext_grad_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.f_ext_grad_spill_tier_3way[2] == 0 else "false") + "; }",
                                 "template <typename T> __host__ __device__ inline size_t F_EXT_GRAD_DQ_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(f_ext_grad_dq_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }"] + [
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t MINV_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.minv_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.minv_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.minv_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t FD_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.fd_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.fd_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.fd_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t ID_DU_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.id_du_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.id_du_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.id_du_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t FD_DU_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.fd_du_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.fd_du_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.fd_du_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 # Tier-aware: at LITE/MINIMAL the FD inner's Minv F-region (6*nv*nv)
                                 # spills to d_workspace, so the smem arena shrinks. Default TIER keeps
                                 # the existing single-arg call sites working.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.integrator_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.integrator_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.integrator_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 # Per-robot tier->placement map for the integrator VALUE path's Minv F-region:
                                 # in smem at spill level 0, in d_workspace (grad section) at level 1.
                                 "template <int TIER> __host__ __device__ constexpr bool INTEGRATOR_MINV_F_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.integrator_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.integrator_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.integrator_spill_tier_3way[2] == 0 else "false") + "; }",
                                 # Tier-aware: at LITE/MINIMAL the s_D_qdd_stage buffer (max_stages*nv*3nv)
                                 # spills to d_workspace, so the smem arena shrinks. Default TIER keeps the
                                 # existing single-arg call sites working.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.integrator_du_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.integrator_du_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.integrator_du_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 # Per-robot tier->placement map for the integrator gradient's s_D_qdd_stage
                                 # buffer: in smem at spill level 0, in d_workspace (grad section) at level 1.
                                 # Per-robot tier->placement maps for the integrator gradient's surgical
                                 # spill ladder. Each buffer's IN_SMEM bool is keyed on RESOURCE_TIER;
                                 # INNER_LEVEL gives the FD-grad inner spill (0 full smem, 1 da_df-band
                                 # selective, 2 whole inner -> d_workspace).
                                 "template <int TIER> __host__ __device__ constexpr bool INTEGRATOR_DU_D_QDD_IN_SMEM() { return (TIER == TIER_SHARED) ? " + _b(self.integrator_du_dqdd_in_smem_per_tier[0]) + " : (TIER == TIER_LITE) ? " + _b(self.integrator_du_dqdd_in_smem_per_tier[1]) + " : " + _b(self.integrator_du_dqdd_in_smem_per_tier[2]) + "; }",
                                 "template <int TIER> __host__ __device__ constexpr bool INTEGRATOR_DU_DAB_IN_SMEM() { return (TIER == TIER_SHARED) ? " + _b(self.integrator_du_dab_in_smem_per_tier[0]) + " : (TIER == TIER_LITE) ? " + _b(self.integrator_du_dab_in_smem_per_tier[1]) + " : " + _b(self.integrator_du_dab_in_smem_per_tier[2]) + "; }",
                                 "template <int TIER> __host__ __device__ constexpr int INTEGRATOR_DU_INNER_LEVEL() { return (TIER == TIER_SHARED) ? " + str(self.integrator_du_inner_level_per_tier[0]) + " : (TIER == TIER_LITE) ? " + str(self.integrator_du_inner_level_per_tier[1]) + " : " + str(self.integrator_du_inner_level_per_tier[2]) + "; }",
                                 # d_workspace sub-offsets (within the per-timestep slot): Dqdd at 0, then dAB, then the inner-spill region.
                                 "template <typename T> __host__ __device__ inline size_t GRID_INTEGRATOR_DU_DAB_OFFSET_BYTES() { return sizeof(T) * static_cast<size_t>(" + str(self._integrator_du_dqdd_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_INTEGRATOR_DU_INNER_OFFSET_BYTES() { return sizeof(T) * static_cast<size_t>(" + str(self._integrator_du_dqdd_count + self._integrator_du_dAB_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t ID_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t MINV_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(minv_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t FD_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fd_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 # Per-tier sizes for forward_dynamics_device (inline-CUDA users only).
                                 # At TIER_SHARED the FD inner s_temp lives in the smem arena; at
                                 # TIER_LITE/MINIMAL the whole arena moves to d_workspace (this is the
                                 # device-path analog of the FD kernel's MINV_F_IN_SMEM lever, which
                                 # surgically spills only the F tail; the device path takes the
                                 # whole-arena route to keep the inline call's smem footprint minimal).
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t FD_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_SHARED)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(fd_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(fd_device_t_count - self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=True)) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t FD_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_SHARED) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=True)) + "); }",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t ABA_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.aba_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.aba_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.aba_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t CRBA_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.crba_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.crba_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.crba_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); "
                                 "}",
                                 "template <typename T> __host__ __device__ constexpr size_t GRID_EE_LINALG_SHARED_BYTES() { return static_cast<size_t>(0); }",
                                 "template <typename T> __host__ __device__ inline size_t EE_POS_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(ee_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 # Phase 3d: tier-aware. PERF/LITE/MINIMAL each report the smem
                                 # bytes their picked spill level needs. Collapsed picks (small
                                 # robots) return identical values across branches.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t DEE_POS_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.ee_grad_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.ee_grad_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.ee_grad_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "}",
                                 # Tier-aware: TIER_SHARED/LITE/MINIMAL each report the smem bytes their
                                 # picked spill level needs. When the picks collapse (small robots) the
                                 # three branches return identical values. Default TIER = TIER_SHARED
                                 # preserves all existing single-arg call sites.
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t D2EE_POS_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.d2ee_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.d2ee_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.d2ee_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); "
                                 "}",
                                 # G2 centroidal quick-wins shared-mem macros (no tier spill).
                                 "template <typename T> __host__ __device__ inline size_t ID_BIAS_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(self.id_bias_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t COM_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(self.com_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t CCRBA_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(self.ccrba_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t ENERGY_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(self.energy_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.idsva_so_body_frame_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.idsva_so_body_frame_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.idsva_so_body_frame_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.idsva_so_world_frame_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.idsva_so_world_frame_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.idsva_so_world_frame_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ inline size_t FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES() { "
                                 "if constexpr (TIER == TIER_SHARED)    return grid_shared_arena_bytes<T>(" + str(self.fdsva_so_t_count_per_tier[0]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else if constexpr (TIER == TIER_LITE) return grid_shared_arena_bytes<T>(" + str(self.fdsva_so_t_count_per_tier[1]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "else                                 return grid_shared_arena_bytes<T>(" + str(self.fdsva_so_t_count_per_tier[2]) + ", TOPOLOGY_HELPERS_COUNT); "
                                 "}",
                                 "// Per-tier scratch sizes for fdsva_so_contract (inline-CUDA users only — the host launchers always use TIER_SHARED).",
                                 "// At TIER_SHARED the 4*NV^3 inner scratch lives in s_temp; at TIER_LITE/MINIMAL it moves to d_workspace, freeing shared memory for the caller's outer kernel.",
                                 "// fdsva_so_contract scratch sizing, keyed on the INNER's placement choice",
                                 "// (SCRATCH_IN_SMEM) rather than a tier — the inner decides placement, the",
                                 "// caller sizes both arenas from these. FDSVA_SO_SCRATCH_IN_SMEM<TIER>()",
                                 "// gives the placement codegen assigned to each tier for THIS robot.",
                                 "template <typename T, bool SCRATCH_IN_SMEM = true> __host__ __device__ constexpr size_t FDSVA_SO_INNER_SMEM_BYTES() { return SCRATCH_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(4*nv**3) + ") : static_cast<size_t>(0); }",
                                 "template <typename T, bool SCRATCH_IN_SMEM = true> __host__ __device__ constexpr size_t FDSVA_SO_INNER_WORKSPACE_BYTES() { return SCRATCH_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(4*nv**3) + "); }",
                                 "// Per-robot tier->placement map: scratch stays in smem at spill levels < 2 (the use_workspace_temp threshold).",
                                 "template <int TIER> __host__ __device__ constexpr bool FDSVA_SO_SCRATCH_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.fdsva_so_spill_tier_3way[0] < 2 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.fdsva_so_spill_tier_3way[1] < 2 else "false") + " : " + ("true" if self.fdsva_so_spill_tier_3way[2] < 2 else "false") + "; }",
                                 "// Inner-controlled placement API (design rollout): each inline inner is keyed on a",
                                 "// placement bool and decides arena pointers itself. *_INNER_{SMEM,WORKSPACE}_BYTES<T, IN_SMEM>",
                                 "// give the two arena sizes; *_<...>_IN_SMEM<TIER>() give the per-robot tier->placement",
                                 "// map codegen assigned (multiple tiers may share a placement on small robots).",
                                 "// --- direct_minv_inner (F-region) ---",
                                 "template <typename T, bool F_IN_SMEM = true> __host__ __device__ constexpr size_t MINV_INNER_SMEM_BYTES() { return sizeof(T) * static_cast<size_t>(" + str(self.gen_direct_minv_inner_no_F_size()) + (" + " + str(6*nv*nv) + " * (F_IN_SMEM ? 1 : 0)") + "); }",
                                 "template <typename T, bool F_IN_SMEM = true> __host__ __device__ constexpr size_t MINV_INNER_WORKSPACE_BYTES() { return F_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(6*nv*nv) + "); }",
                                 "template <int TIER> __host__ __device__ constexpr bool MINV_F_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.minv_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.minv_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.minv_spill_tier_3way[2] == 0 else "false") + "; }",
                                 "// --- forward_dynamics_inner (internal Minv F-region) ---",
                                 "template <typename T, bool MINV_F_IN_SMEM = true> __host__ __device__ constexpr size_t FD_INNER_SMEM_BYTES() { return MINV_F_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=True)) + ") : sizeof(T) * static_cast<size_t>(" + str(self.gen_forward_dynamics_inner_temp_mem_size(minv_f_in_smem=False)) + "); }",
                                 "template <typename T, bool MINV_F_IN_SMEM = true> __host__ __device__ constexpr size_t FD_INNER_WORKSPACE_BYTES() { return MINV_F_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(6*nv*nv) + "); }",
                                 "template <int TIER> __host__ __device__ constexpr bool FD_MINV_F_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.fd_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.fd_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.fd_spill_tier_3way[2] == 0 else "false") + "; }",
                                 # The integrator value path's only inner scratch is the FD inner itself,
                                 # so its arena sizes mirror FD_INNER_* exactly: when MINV_F_IN_SMEM the
                                 # 6*NV*NV F-region is in s_temp, else it spills to d_workspace. The
                                 # per-robot tier->placement map is INTEGRATOR_MINV_F_IN_SMEM<TIER> (above).
                                 "// --- integrator_inner (forwards the FD inner's Minv F-region lever) ---",
                                 "template <typename T, bool MINV_F_IN_SMEM = true> __host__ __device__ constexpr size_t INTEGRATOR_INNER_SMEM_BYTES() { return MINV_F_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(self.gen_integrator_inner_temp_mem_size(minv_f_in_smem=True)) + ") : sizeof(T) * static_cast<size_t>(" + str(self.gen_integrator_inner_temp_mem_size(minv_f_in_smem=False)) + "); }",
                                 "template <typename T, bool MINV_F_IN_SMEM = true> __host__ __device__ constexpr size_t INTEGRATOR_INNER_WORKSPACE_BYTES() { return MINV_F_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(6*nv*nv) + "); }",
                                 "// --- aba_inner (scratch band, surgical-spill ladder) ---",
                                 "// Levels: 0=full (smem), 1=surgical (hot smem + cold d_cold), 2=workspace (whole band global).",
                                 "// TEMP_IN_SMEM is false only at the level-2 (workspace) rung; COLD_IN_SMEM is false only at the level-1 (surgical) rung.",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t ABA_INNER_SMEM_BYTES() { return TEMP_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(self.gen_aba_inner_temp_mem_size()) + ") : static_cast<size_t>(0); }",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t ABA_INNER_WORKSPACE_BYTES() { return TEMP_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(self.gen_aba_inner_temp_mem_size()) + "); }",
                                 "template <typename T> __host__ __device__ constexpr size_t ABA_INNER_COLD_BYTES() { return sizeof(T) * static_cast<size_t>(" + str(self._aba_inner_cold_count) + "); }",
                                 "template <int TIER> __host__ __device__ constexpr bool ABA_TEMP_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("false" if self.aba_spill_tier_3way[0] == 2 else "true") + " : (TIER == TIER_LITE) ? " + ("false" if self.aba_spill_tier_3way[1] == 2 else "true") + " : " + ("false" if self.aba_spill_tier_3way[2] == 2 else "true") + "; }",
                                 "template <int TIER> __host__ __device__ constexpr bool ABA_COLD_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("false" if self.aba_spill_tier_3way[0] == 1 else "true") + " : (TIER == TIER_LITE) ? " + ("false" if self.aba_spill_tier_3way[1] == 1 else "true") + " : " + ("false" if self.aba_spill_tier_3way[2] == 1 else "true") + "; }",
                                 "// --- crba_inner (scratch band) ---",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t CRBA_INNER_SMEM_BYTES() { return TEMP_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(self.gen_crba_inner_temp_mem_size()) + ") : static_cast<size_t>(0); }",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t CRBA_INNER_WORKSPACE_BYTES() { return TEMP_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(self.gen_crba_inner_temp_mem_size()) + "); }",
                                 "template <int TIER> __host__ __device__ constexpr bool CRBA_TEMP_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.crba_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.crba_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.crba_spill_tier_3way[2] == 0 else "false") + "; }",
                                 "// --- end_effector_pose_gradient_inner (chain workspace) ---",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t EE_GRAD_INNER_SMEM_BYTES() { return TEMP_IN_SMEM ? sizeof(T) * static_cast<size_t>(" + str(self.gen_end_effector_pose_gradient_inner_temp_mem_size()) + ") : static_cast<size_t>(0); }",
                                 "template <typename T, bool TEMP_IN_SMEM = true> __host__ __device__ constexpr size_t EE_GRAD_INNER_WORKSPACE_BYTES() { return TEMP_IN_SMEM ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(self.gen_end_effector_pose_gradient_inner_temp_mem_size()) + "); }",
                                 "template <int TIER> __host__ __device__ constexpr bool EE_GRAD_TEMP_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.ee_grad_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.ee_grad_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.ee_grad_spill_tier_3way[2] == 0 else "false") + "; }",
                                 "// --- end_effector_pose_gradient_hessian_inner (large nv^2 d2eePos output) ---",
                                 "// Per-tier placement of the d2ee inner's OUTPUT s_d2eePos: true => smem, false => d_workspace (which the kernel sets to d_d2eePos directly).",
                                 "template <int TIER> __host__ __device__ constexpr bool D2EE_OUT_IN_SMEM() { return (TIER == TIER_SHARED) ? " + ("true" if self.d2ee_spill_tier_3way[0] == 0 else "false") + " : (TIER == TIER_LITE) ? " + ("true" if self.d2ee_spill_tier_3way[1] == 0 else "false") + " : " + ("true" if self.d2ee_spill_tier_3way[2] == 0 else "false") + "; }",
                                 "// Per-tier sizes for forward_dynamics_gradient_device (inline-CUDA users only). At TIER_SHARED the temp scratch arena lives in s_temp; at TIER_LITE/MINIMAL it moves to d_workspace, freeing roughly " + str(fd_du_temp_count) + "*sizeof(T) bytes of smem.",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t FD_DU_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_SHARED)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(fd_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(fd_du_device_t_count - fd_du_temp_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t FD_DU_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_SHARED) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(fd_du_temp_count) + "); }",
                                 "// Per-tier sizes for end_effector_pose_gradient_hessian_device (inline-CUDA users only). At TIER_SHARED the smem arena keeps only the FD scratch + s_Xhom; at TIER_LITE/MINIMAL the device contract is unchanged (smem arena is the same -- the caller-provided s_d2eePos is what shifts), and the inner writes its " + str(d2ee_output_count) + "*sizeof(T) output bytes to d_workspace instead.",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t D2EE_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return grid_shared_arena_bytes<T>(" + str(d2ee_inner_temp_count + XHom_size) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t D2EE_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_SHARED) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(d2ee_output_count) + "); }",
                                 "// Per-tier sizes for inverse_dynamics_gradient_device (inline-CUDA users only). At TIER_SHARED temp lives in s_temp; at TIER_LITE/MINIMAL it moves to d_workspace, freeing " + str(id_du_temp_count) + "*sizeof(T) bytes of smem.",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t ID_DU_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_SHARED)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(id_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(id_du_device_t_count - id_du_temp_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t ID_DU_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_SHARED) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(id_du_temp_count) + "); }",
                                 "// Per-tier sizes for idsva_so_device (inline-CUDA users only). At TIER_SHARED temp lives in s_temp; at TIER_LITE/MINIMAL it moves to d_workspace, freeing " + str(idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + "*sizeof(T) bytes of smem. Frame picked at codegen time: " + ("world_frame" if self.robot.floating_base else "body_frame") + ".",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t IDSVA_SO_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_SHARED)",
                                 "        ? grid_shared_arena_bytes<T>(" + str((idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + XI_size) + ", TOPOLOGY_HELPERS_COUNT)",
                                 "        : grid_shared_arena_bytes<T>(" + str(XI_size) + ", TOPOLOGY_HELPERS_COUNT);",
                                 "}",
                                 "template <typename T, int TIER = GRID_DEFAULT_RESOURCE_TIER> __host__ __device__ constexpr size_t IDSVA_SO_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_SHARED) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP() { return sizeof(T) * static_cast<size_t>(" + str(grad_spill_workspace_t_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_SO_WORKSPACE_BYTES_PER_TIMESTEP() { return sizeof(T) * static_cast<size_t>(" + str(so_workspace_t_count) + "); }",
                                 # Phase 3e: sized for MINIMAL tier's spill (max across PERF/LITE/MINIMAL).
                                 # Even if PERF doesn't spill df_du/Minv, MINIMAL might — the workspace
                                 # allocation has to cover MINIMAL's needs at all times.
                                 "template <typename T> __host__ __device__ inline size_t GRID_FDSVA_SO_SPILL_BYTES_PER_TIMESTEP() { return sizeof(T) * static_cast<size_t>(" + str(3*nv*nv if any(p >= 4 for p in getattr(self, 'fdsva_so_spill_tier_3way', (0, 0, 0))) else 0) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_FDSVA_SO_SPILL_OFFSET_BYTES() { return GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_BYTES_PER_TIMESTEP<T>(); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_WORKSPACE_BYTES_PER_TIMESTEP() { return GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_FDSVA_SO_SPILL_BYTES_PER_TIMESTEP<T>(); }",
                                 "template <typename T> __host__ __device__ inline gridSharedTier GRID_ID_DU_SHARED_TIER() { return static_cast<gridSharedTier>(GRID_ID_DU_SHARED_TIER_VALUE); }",
                                 "template <typename T> __host__ __device__ inline gridSharedTier GRID_FD_DU_SHARED_TIER() { return static_cast<gridSharedTier>(GRID_FD_DU_SHARED_TIER_VALUE); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES() { return GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP<T>(); }",
                                 # Phase 3a: Minv-F lives at offset 0 of the grad section when spilled.
                                 # Safe to overlap with id_du_spill region because Minv finishes before
                                 # id_du_grad starts in any kernel that composes both.
                                 "template <typename T> __host__ __device__ inline size_t GRID_MINV_F_WORKSPACE_OFFSET_BYTES() { return static_cast<size_t>(0); }",
                                 # ABA surgical cold sub-buffer reuses the SO/grad workspace band base (ABA
                                 # never runs concurrently with SO/grad). The cold band (ABA_INNER_COLD_BYTES)
                                 # is far smaller than GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP, so it fits at
                                 # offset 0 without growing GRID_WORKSPACE_BYTES_PER_TIMESTEP.
                                 "template <typename T> __host__ __device__ inline size_t GRID_ABA_COLD_OFFSET_BYTES() { return static_cast<size_t>(0); }",
                                 # D2EE no longer uses a per-timestep d_workspace slice (the spilled
                                 # s_d2eePos is written directly into d_d2eePos); these offset macros are
                                 # retained as 0 for backward compatibility with any inline-CUDA caller
                                 # pattern that still references them. New code should not use them.
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_TEMP_OFFSET_BYTES() { return static_cast<size_t>(0); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_D2XHOM_OFFSET_BYTES() { return static_cast<size_t>(0); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_D2EETEMP_OFFSET_BYTES() { return static_cast<size_t>(0); }",
                                 # Phase 3d: EE_POSE_GRAD reuses the SO section (the kernels don't
                                 # run concurrently — d_workspace bytes are safely repurposed). When
                                 # the MINIMAL tier spills dXmatsHom, it sits before the temp arena.
                                 "template <typename T> __host__ __device__ inline size_t GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES() { return GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>(); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_EE_GRAD_WORKSPACE_TEMP_OFFSET_BYTES() { return GRID_EE_GRAD_WORKSPACE_DXHOM_OFFSET_BYTES<T>() + (GRID_EE_GRAD_USES_WORKSPACE_DXHOM ? sizeof(T) * static_cast<size_t>(DXHOM_T_COUNT) : 0); }",
                                 "template <typename T> __host__ __device__ inline bool grid_selected_shared_memory_fits() { return ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES && FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES && (!GRID_GENERATES_D2EE || D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES) && (!GRID_GENERATES_IDSVA_SO_BODY_FRAME || IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES) && (!GRID_GENERATES_FDSVA_SO || FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES); }",
                                 "// __forceinline__ used throughout the xhom helper chain so ptxas folds these into the",
                                 "// inner kernels at all opt levels. For fixed-base the body of grid_q_index_affects_joint is",
                                 "// the trivial `q_index == joint_id` check that pre-GLASS callsites used directly.",
                                 "__host__ __device__ __forceinline__ bool grid_q_index_affects_joint(const int q_index, const int joint_id) {",
                                 ("    if (joint_id == 0) { return q_index >= 0 && q_index < 7; } return q_index == joint_id + 6;" if self.robot.floating_base else "    return q_index == joint_id;"),
                                 "}",
                                 "__host__ __device__ __forceinline__ int grid_d2xhom_offset(const int q_index_i, [[maybe_unused]] const int q_index_j) {",
                                 ("    return 16 * (q_index_i * NUM_JOINTS + q_index_j);" if self.robot.floating_base else "    return 16 * q_index_i;"),
                                 "}",
                                 "template <typename T>",
                                 "__device__ __forceinline__ const T *grid_xhom_or_dxhom_ptr(const T *s_Xhom, const T *s_dXhom, const int q_index, const int joint_id) {",
                                 "    return grid_q_index_affects_joint(q_index, joint_id) ? &s_dXhom[16 * q_index] : &s_Xhom[16 * joint_id];",
                                 "}",
                                 "template <typename T>",
                                 "__device__ __forceinline__ const T *grid_xhom_or_dxhom_or_d2xhom_ptr(const T *s_Xhom, const T *s_dXhom, const T *s_d2Xhom, const int q_index_i, const int q_index_j, const int joint_id) {",
                                 "    const bool i_affects = grid_q_index_affects_joint(q_index_i, joint_id);",
                                 "    const bool j_affects = grid_q_index_affects_joint(q_index_j, joint_id);",
                                 "    if (i_affects && j_affects) { return &s_d2Xhom[grid_d2xhom_offset(q_index_i, q_index_j)]; }",
                                 "    if (i_affects) { return &s_dXhom[16 * q_index_i]; }",
                                 "    if (j_affects) { return &s_dXhom[16 * q_index_j]; }",
                                 "    return &s_Xhom[16 * joint_id];",
                                 "}",
                                 "template <typename T, bool USE_DA_DF_SPILL>",
                                 "__device__ inline T *grid_id_du_temp_ptr(T *s_temp, T *d_temp_spill, int index) {",
                                 "    if (!USE_DA_DF_SPILL) { return &s_temp[index]; }",
                                 "    if (index >= ID_DU_TEMP_SPILL_START && index < ID_DU_TEMP_SPILL_END) {",
                                 "        return &d_temp_spill[index - ID_DU_TEMP_SPILL_START];",
                                 "    }",
                                 "    if (index >= ID_DU_TEMP_SPILL_END) {",
                                 "        return &s_temp[index - ID_DU_TEMP_SPILL_COUNT];",
                                 "    }",
                                 "    return &s_temp[index];",
                                 "}",
                                 ""])
        # then the structs
        # first add the struct
        self.gen_add_code_line("// Define custom structs")
        self.gen_add_code_lines(["template <typename T>", \
                                 "struct robotModel {", \
                                 "    T *d_XImats;", \
                                 "    int *d_topology_helpers;", \
                                 "};"])
        self.gen_add_code_lines(["template <typename T, gridDataKind KIND = GRID_DATA_ALL>", \
                                 "struct gridData {", \
                                 "    // GPU INPUTS", \
                                 "    T *d_q_qd_u;", \
                                 "    T *d_q_qd;", \
                                 "    T *d_q;", \
                                 # external forces: body-major 6*NUM_BODIES local-frame, per timestep (zeroed by default)
                                 "    T *d_f_ext;", \
                                 "    // CPU INPUTS", \
                                 "    T *h_q_qd_u;", \
                                 "    T *h_q_qd;", \
                                 "    T *h_q;", \
                                 "    T *h_f_ext;", \
                                 "    // GPU OUTPUTS", \
                                 "    T *d_c;", \
                                 "    T *d_Minv;", \
                                 "    T *d_qdd;", \
                                 "    T *d_M;", \
                                 "    T *d_dc_du;", \
                                 "    T *d_df_du;", \
                                 # f_ext gradient column (section A): dtau/dfext = -J^T,
                                 # dqdd/dfext = M^-1 J^T; each nv x (6*NB), body-major.
                                 "    T *d_dtau_dfext;",
                                 "    T *d_dqdd_dfext;",
                                 "    T *d_did_du_dfext;  // -dJ^T/dq = d(id_du)/dfext, nv*6NB*nv (both base modes)"]
                                 + [
                                 "    T *d_eePos;", \
                                 "    T *d_deePos;", \
                                 "    T *d_d2eePos;", \
                                 "    unsigned char *d_workspace;", \
                                 # idsva_so - d2tau_dq2, d2tau_dqd2, d2tau_dvdq, dM_dq
                                 "    T *d_idsva_so;", \
                                 # fdsva_so - d2a_dq2, d2a_dv2, d2a_dvdq, d2a_dtdq
                                 "    T *d_df2;", \
                                 # integrator outputs
                                 "    T *d_x_kp1;", \
                                 "    T *d_dAB;", \
                                 # G2 centroidal quick-wins outputs
                                 "    T *d_com;", \
                                 "    T *d_ccrba;", \
                                 "    T *d_energy;", \
                                 "    // CPU OUTPUTS", \
                                 "    T *h_c;", \
                                 "    T *h_Minv;", \
                                 "    T *h_qdd;", \
                                 "    T *h_M;", \
                                 "    T *h_dc_du;", \
                                 "    T *h_df_du;", \
                                 "    T *h_dtau_dfext;",
                                 "    T *h_dqdd_dfext;",
                                 "    T *h_did_du_dfext;  // -dJ^T/dq, nv*6NB*nv (both base modes)"]
                                 + [
                                 "    T *h_eePos;", \
                                 "    T *h_deePos;", \
                                 "    T *h_d2eePos;", \
                                 # idsva_so - d2tau_dq2, d2tau_dqd2, d2tau_dvdq, dM_dq
                                 "    T *h_idsva_so;", \
                                 # fdsva_so - d2a_dq2, d2a_dv2, d2a_dvdq, d2a_dtdq
                                 "    T *h_df2;", \
                                 # integrator outputs
                                 "    T *h_x_kp1;", \
                                 "    T *h_dAB;", \
                                 # G2 centroidal quick-wins outputs
                                 "    T *h_com;", \
                                 "    T *h_ccrba;", \
                                 "    T *h_energy;", \
                                 "};"])

    def gen_init_gridData(self):
        code_lines = (["gridData<T, KIND> *hd_data = (gridData<T, KIND> *)calloc(1, sizeof(gridData<T, KIND>));",
                      "const bool needs_dynamics = KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS;",
                      "const bool needs_kinematics = KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS;",
                      "// input variables used by dynamics and/or kinematics",
                      "if (needs_dynamics || needs_kinematics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_q_qd_u, 3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_q, NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_q_qd_u = (T *)malloc(3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_q = (T *)malloc(NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    // external forces (body-major 6*NUM_BODIES local-frame); zeroed so the", \
                      "    // default (no-fext) path subtracts nothing. Users overwrite h_f_ext and", \
                      "    // copy to d_f_ext to apply external forces.", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_f_ext, 6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMemset(hd_data->d_f_ext, 0, 6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_f_ext = (T *)calloc(6*NUM_BODIES*NUM_TIMESTEPS, sizeof(T));", \
                      "}", \
                      "if (needs_dynamics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_q_qd, 2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_q_qd = (T *)malloc(2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "}", \
                      "// dynamics outputs and fallback workspace", \
                      "if (needs_dynamics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_c, NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_Minv, NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_qdd, NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_M, NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_dc_du, NUM_JOINTS*2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_df_du, NUM_JOINTS*2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    // f_ext gradient column (section A): dtau/dfext, dqdd/dfext are each nv x (6*NB)", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_dtau_dfext, NUM_VEL*6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_dqdd_dfext, NUM_VEL*6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_dtau_dfext = (T *)malloc(NUM_VEL*6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T));",
                      "    hd_data->h_dqdd_dfext = (T *)malloc(NUM_VEL*6*NUM_BODIES*NUM_TIMESTEPS*sizeof(T));",
                      "    // f_ext A.3: -dJ^T/dq = d(id_du)/dfext, nv*6NB*nv (both base modes)",
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_did_du_dfext, NUM_VEL*6*NUM_BODIES*NUM_VEL*NUM_TIMESTEPS*sizeof(T)));",
                      "    hd_data->h_did_du_dfext = (T *)malloc(NUM_VEL*6*NUM_BODIES*NUM_VEL*NUM_TIMESTEPS*sizeof(T));"]
                      + [
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_idsva_so, SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_df2, SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_workspace, GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*NUM_TIMESTEPS));", \
                      "    // Phase 3a/b/c/e: L2-pin d_workspace for its lifetime. Spilled buffers", \
                      "    // (Minv-F, FD's Minv-F, ABA's inner scratch, FDSVA_SO's df_du/Minv) are", \
                      "    // recursion-hot — L2 pinning narrows the smem→HBM gap to smem→L2.", \
                      "    gpuErrchk(grid_begin_l2_persisting(0, hd_data->d_workspace, GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*NUM_TIMESTEPS));", \
                      "    hd_data->h_c = (T *)malloc(NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_Minv = (T *)malloc(NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_M = (T *)malloc(NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_qdd = (T *)malloc(NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_dc_du = (T *)malloc(NUM_JOINTS*2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_df_du = (T *)malloc(NUM_JOINTS*2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_idsva_so = (T *)malloc(SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_df2 = (T *)malloc(SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_x_kp1, 2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_dAB, 2*NUM_JOINTS*3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_x_kp1 = (T *)malloc(2*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_dAB = (T *)malloc(2*NUM_JOINTS*3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "}", \
                      "// kinematics outputs", \
                      "if (needs_kinematics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_eePos, 6*NUM_EES*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_deePos, 6*NUM_EES*NUM_VEL*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_d2eePos, 6*NUM_EES*NUM_VEL*NUM_VEL*NUM_TIMESTEPS*sizeof(T)));", \
                      "    if ((GRID_D2EE_USES_WORKSPACE_TEMP || GRID_EE_GRAD_USES_WORKSPACE_TEMP) && hd_data->d_workspace == nullptr) {gpuErrchk(cudaMalloc((void**)&hd_data->d_workspace, GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*NUM_TIMESTEPS));}", \
                      "    hd_data->h_eePos = (T *)malloc(6*NUM_EES*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_deePos = (T *)malloc(6*NUM_EES*NUM_VEL*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_d2eePos = (T *)malloc(6*NUM_EES*NUM_VEL*NUM_VEL*NUM_TIMESTEPS*sizeof(T));", \
                      "}", \
                      "// G2 centroidal quick-wins outputs (com: 3+3*NV ; ccrba: 6*NV+6 ; energy: 3)", \
                      "if (needs_dynamics || needs_kinematics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_com, (3+3*NUM_VEL)*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_ccrba, (6*NUM_VEL+6)*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_energy, 3*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_com = (T *)malloc((3+3*NUM_VEL)*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_ccrba = (T *)malloc((6*NUM_VEL+6)*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_energy = (T *)malloc(3*NUM_TIMESTEPS*sizeof(T));", \
                      "}", \
                      "return hd_data;"])
        # generate as templated or not function
        self.gen_add_func_doc("Allocated device and host memory for all computations",
                              [], [], "A pointer to the gridData struct of pointers")
        self.gen_add_code_line("template <typename T, int NUM_TIMESTEPS, gridDataKind KIND = GRID_DATA_ALL>")
        self.gen_add_code_line("__host__")
        self.gen_add_code_line("gridData<T, KIND> *init_gridData(){", True)
        self.gen_add_code_lines(code_lines)
        self.gen_add_end_function()
        self.gen_add_func_doc("Allocated device and host memory for all computations",
                              [], ["Max number of timesteps in the trajectory"], "A pointer to the gridData struct of pointers")
        self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
        self.gen_add_code_line("__host__")
        self.gen_add_code_line("gridData<T, KIND> *init_gridData(int NUM_TIMESTEPS){", True)
        self.gen_add_code_lines(code_lines)
        self.gen_add_end_function()

    # Manifest of every algorithm kernel that needs cudaFuncSetAttribute
    # (MaxDynamicSharedMemorySize). Each entry:
    #   (algo_label, algo_short_name, gate_attr, bytes_macro, [(kernel_name, signature), ...])
    #
    # algo_short_name matches the keys used in self.generated_algorithms (see
    # _normalize_codegen_algorithms). gate_attr is the legacy `self.generate_*`
    # bool — when present we honor it for back-compat; otherwise we rely on
    # algo_short_name membership in self.generated_algorithms.
    #
    # Why we apply this to EVERY kernel, not just the historically-large ones:
    # without cudaFuncSetAttribute(MaxDynamicSharedMemorySize, BYTES), a kernel
    # whose dynamic shared mem at runtime exceeds the device default per-block
    # limit (48 KB on most consumer NVIDIA GPUs incl. sm_8x) launches and fails
    # silently with cudaErrorInvalidConfiguration. The launch error doesn't
    # propagate through cudaDeviceSynchronize() reliably, so timings come back
    # as bogus ~0 us values. We hit this on g1 floating where ABA / FD / MINV
    # need 52-57 KB shared mem. The call is a no-op when BYTES is already
    # under the device default.
    # Default False so the f_ext A.3 (-dJ^T/dq) kernel is registered ONLY when
    # gen_f_ext_gradient actually emitted it. gen_f_ext_gradient now emits it for
    # BOTH base modes (fixed: scalar FD; floating: SE(3) Lie-group root retract),
    # so the instance attr is set True whenever f_ext_grad is generated.
    _f_ext_grad_dq_emitted = False

    KERNEL_ATTR_MANIFEST = [
        # (algo_label, algo_short, gate_attr, bytes_macro, [(kernel_name<T>, signature), ...])
        ("inverse_dynamics", "id", None, "ID_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("inverse_dynamics_kernel<T>",
             "void (*)(T *, const T *, const int, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel<T>",
             "void (*)(T *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
        ]),
        ("direct_minv", "minv", None, "MINV_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("direct_minv_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
            ("direct_minv_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("forward_dynamics", "fd", None, "FD_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("forward_dynamics_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
        ]),
        ("aba", "aba", None, "ABA_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("aba_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
            ("aba_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
        ]),
        ("crba", "crba", None, "CRBA_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("crba_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("crba_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("end_effector_pose", "ee_pose", None, "EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("end_effector_pose_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("end_effector_pose_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("end_effector_pose_gradient", "ee_pose_gradient", None, "DEE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("end_effector_pose_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
            ("end_effector_pose_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("inverse_dynamics_gradient", "id_du", "generate_id_du", "ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("inverse_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
        ]),
        ("forward_dynamics_gradient", "fd_du", "generate_fd_du", "FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("forward_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const T *, T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, T *, const robotModel<T> *, const T, const int)"),
        ]),
        # g1-spill: f_ext_gradient_kernel gained `unsigned char *d_workspace` as its
        # 3rd arg (after the two outputs) so s_dqdd_dfext can spill there at LITE/MINIMAL.
        ("f_ext_gradient", "f_ext_grad", None, "F_EXT_GRAD_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("f_ext_gradient_kernel<T>",
             "void (*)(T *, T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
            ("f_ext_gradient_kernel_single_timing<T>",
             "void (*)(T *, T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        # A.3 (-dJ^T/dq): own kernel + smem macro, fixed-base only. Gated on the
        # instance attr _f_ext_grad_dq_emitted (set True only when the kernel is
        # actually emitted) so the floating-base header — which has neither the
        # kernel nor the F_EXT_GRAD_DQ_* macro — never references them.
        ("f_ext_gradient_dq", "f_ext_grad_dq", "_f_ext_grad_dq_emitted", "F_EXT_GRAD_DQ_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("f_ext_gradient_dq_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("f_ext_gradient_dq_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        # E1 joint-torque regressor: Y is nv x 10*NUM_BODIES, can exceed the 48 KB
        # default dynamic-smem cap on big robots (g1: ~55 KB), so it MUST opt in.
        ("regressor", "regressor", None, "INVERSE_DYNAMICS_REGRESSOR_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("inverse_dynamics_regressor_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_regressor_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        # FD param gradient dqdd/dpi = -Minv . Y: output is nv x 10*NUM_BODIES (same
        # size class as the regressor), can exceed the 48 KB default cap; opt in.
        # g1-spill: fd_parameter_gradient_kernel gained `unsigned char *d_workspace`
        # as its 2nd arg (after d_dqdd_dpi) so s_Y can spill there at LITE/MINIMAL.
        ("fd_parameter_gradient", "fd_parameter_gradient", None, "FD_PARAMETER_GRADIENT_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("fd_parameter_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("fd_parameter_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("idsva_so_body_frame", "idsva_so_body_frame", "generate_idsva_so_body_frame", "IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("idsva_so_body_frame_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("idsva_so_body_frame_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        # world-frame single-pass alternative (opt-in via enable_idsva_so_world_frame).
        # Now takes d_workspace (unified signature; cold buffers spill there at LITE/MINIMAL).
        # Uses its own shared-mem macro (~25 KB for g1 vs shim's ~162 KB).
        ("idsva_so_world_frame", "idsva_so_world_frame", "generate_idsva_so_world_frame",
         "IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("idsva_so_world_frame_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("idsva_so_world_frame_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("fdsva_so", "fdsva_so", "generate_fdsva_so", "FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("fdsva_so_kernel<T>",
             "void (*)(T *, const T *, const int, unsigned char *, T *, const robotModel<T> *, const T, const int)"),
            ("fdsva_so_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, unsigned char *, T *, const robotModel<T> *, const T, const int)"),
        ]),
        # Integrator kernels are templated on IntegratorType IT (a non-type param
        # that does not change the function signature). Each IT is a distinct
        # __global__ instantiation, so cudaFuncSetAttribute must run for ALL of
        # them — otherwise a non-Euler IT whose floating-base arena exceeds the
        # 48 KB device default launches with cudaErrorInvalidValue while Euler
        # (the only one historically registered) succeeds.
        ("integrator", "integrator", None, "INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            (f"integrator_kernel{suffix}<T, IntegratorType::{it}>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const T, const int)")
            for suffix in ("", "_single_timing")
            for it in ("EULER", "SEMI_IMPLICIT_EULER", "MIDPOINT", "RK3", "RK4")
        ]),
        ("integrator_gradient", "integrator_gradient", None, "INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            (f"integrator_gradient_kernel{suffix}<T, IntegratorType::{it}>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const T, const int)")
            for suffix in ("", "_single_timing")
            for it in ("EULER", "SEMI_IMPLICIT_EULER", "MIDPOINT", "RK3", "RK4")
        ]),
        ("integrator_gradient_with_x_kp1", "integrator_with_gradient", None, "INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            (f"integrator_gradient_with_x_kp1_kernel{suffix}<T, IntegratorType::{it}>",
             "void (*)(T *, T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const T, const int)")
            for suffix in ("", "_single_timing")
            for it in ("EULER", "SEMI_IMPLICIT_EULER", "MIDPOINT", "RK3", "RK4")
        ]),
        # ee_pose_hessian is special: only emitted when its shared-mem fits the
        # GRID_CUDA_TARGET_SHARED_MEM_BYTES budget at compile time. The runtime
        # guard wraps the cudaFuncSetAttribute call.
        ("end_effector_pose_gradient_hessian", "ee_pose_hessian", "generate_ee_pose_hessian",
         "D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("end_effector_pose_gradient_hessian_kernel<T>",
             "void (*)(T *, T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
            ("end_effector_pose_gradient_hessian_kernel_single_timing<T>",
             "void (*)(T *, T *, unsigned char *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        # G2 centroidal quick-wins. gravity/nonlinear_effects gate on `id`
        # (RNEA bias wrappers); com/ccrba/energy gate on `ee_pose` (homogeneous-
        # transform world-frame machinery). algo_short keys an entry that is in
        # generated_algorithms exactly when the dep is present.
        ("generalized_gravity", "id", None, "ID_BIAS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("generalized_gravity_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("generalized_gravity_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("nonlinear_effects", "id", None, "ID_BIAS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("nonlinear_effects_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("nonlinear_effects_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("com", "ee_pose", None, "COM_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("com_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("com_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("ccrba", "ee_pose", None, "CCRBA_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("ccrba_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("ccrba_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("energy", "ee_pose", None, "ENERGY_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("energy_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("energy_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
    ]

    def gen_init_close_grid(self):
        # set the max shared mem to account for large robots and allocate streams
        MAX_STREAMS = 3 # max needed in any of our functions
        # ----- init_grid_kernel_attrs<T>(): cudaFuncSetAttribute for every kernel --
        # Split out from init_grid so per-algo TU callers can invoke ONLY the
        # attribute-setting part (without stream allocation), needed because
        # cudaFuncSetAttribute operates on the TU-local host stub. The per-algo
        # TU split (P6-7b) calls this from a static initializer in each
        # measure_X_*_entry so its stubs get the attribute set.
        self.gen_add_func_doc("Set MaxDynamicSharedMemorySize for every algorithm kernel "
                              "(callable from any TU; idempotent). __forceinline__ is "
                              "REQUIRED so the &kernel<T> expressions resolve to the "
                              "CALLING TU's host stubs — otherwise the linker merges this "
                              "function across TUs and we set the attribute on one TU's "
                              "stubs while the launch goes through a different TU's.",
                              [], [], None)
        self.gen_add_code_line("template <typename T>")
        self.gen_add_code_line("__host__ __forceinline__")
        self.gen_add_code_line("void init_grid_kernel_attrs(){", True)
        attr_lines = ["// enable opt-in dynamic shared memory for every algorithm kernel",
                      "// Gate registration on the DEVICE opt-in max (not the codegen target):",
                      "// grid_check_dynamic_shared_memory_bytes and the bench's",
                      "// grid_kernel_fits_device both use the device cap, so registering only up",
                      "// to the smaller GRID_CUDA_TARGET_SHARED_MEM_BYTES left kernels in",
                      "// (target, device-max] checkable+launchable but UNregistered -> launching",
                      "// them failed with cudaErrorInvalidValue (e.g. the floating idsva_so",
                      "// body-frame diagnostic ~101 KB on g1). Keying on the device max keeps",
                      "// registration, the fit-check, and the launch-skip in lockstep.",
                      "size_t _grid_smem_max = 0; gpuErrchk(grid_get_max_dynamic_shared_memory_bytes(&_grid_smem_max));"]
        generated_set = getattr(self, "generated_algorithms", None)
        alias_counter = 0
        for entry in self.KERNEL_ATTR_MANIFEST:
            algo_label, algo_short, gate_attr, bytes_macro, kernels = entry
            # Honor the legacy generate_* gate when present; otherwise fall
            # back to membership in generated_algorithms; if neither is set
            # (legacy callers), assume the algo is generated.
            if gate_attr is not None and not getattr(self, gate_attr, True):
                continue
            if gate_attr is None and generated_set is not None and algo_short not in generated_set:
                continue
            # G2 centroidal kinematics-domain families are not emitted for mimic
            # robots (their per-body Jacobian fold isn't mimic-reduced yet), so
            # skip registering their (nonexistent) kernels there.
            if algo_label in ("com", "ccrba", "energy") and self.robot_has_mimic_joints():
                continue
            # Wrap EVERY kernel's attribute registration in a compile-time-
            # resolvable size guard so init_grid never hard-aborts when a kernel
            # literally can't fit a device even with cudaFuncSetAttribute (e.g.
            # the floating-base idsva_so_body_frame *diagnostic* frame on h1_2 at
            # ~168 KB, or the integrator value/gradient kernels at ~103-228 KB on
            # big floating-base robots). Such kernels simply go unregistered;
            # they aren't the dispatched production path, and any code that DOES
            # launch them still runs `grid_check_dynamic_shared_memory_bytes` at
            # the host wrapper, so the fit check + attribute setup stay in
            # lockstep at the actual use site. Small kernels are always under the
            # target, so the guard is a no-op for them.
            attr_lines.append(f"if ({bytes_macro} <= _grid_smem_max) {{")
            attr_lines.append(f"    gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"{algo_label}\", {bytes_macro}));")
            for kernel_name, signature in kernels:
                alias = f"_grid_kern_alias_{alias_counter}"
                alias_counter += 1
                attr_lines.append(f"    auto {alias} = static_cast<{signature}>(&{kernel_name});")
                attr_lines.append(f"    gpuErrchk(cudaFuncSetAttribute({alias}, cudaFuncAttributeMaxDynamicSharedMemorySize, {bytes_macro}));")
            attr_lines.append("}")
        self.gen_add_code_lines(attr_lines)
        self.gen_add_end_function()

        # ----- init_grid<T>(): full init = attrs + streams (the original API) ----
        self.gen_add_func_doc("Sets MaxDynamicSharedMemorySize for every algorithm kernel and initializes streams for host functions", \
                              [], [], "A pointer to the array of streams")
        self.gen_add_code_line("template <typename T>")
        self.gen_add_code_line("__host__")
        self.gen_add_code_line("cudaStream_t *init_grid(){", True)
        init_lines = ["init_grid_kernel_attrs<T>();",
                      "gpuErrchk(cudaDeviceSynchronize());",
                      "// allocate streams",
                      "cudaStream_t *streams = (cudaStream_t *)malloc(" + str(MAX_STREAMS) + "*sizeof(cudaStream_t));",
                      "int priority, minPriority, maxPriority;",
                      "gpuErrchk(cudaDeviceGetStreamPriorityRange(&minPriority, &maxPriority));",
                      "for(int i=0; i<" + str(MAX_STREAMS) + "; i++){",
                      "    int adjusted_max = maxPriority - i; priority = adjusted_max > minPriority ? adjusted_max : minPriority;",
                      "    gpuErrchk(cudaStreamCreateWithPriority(&(streams[i]),cudaStreamNonBlocking,priority));",
                      "}", "return streams;"]
        self.gen_add_code_lines(init_lines)
        self.gen_add_end_function()
        # free the streams and all allocated data
        self.gen_add_func_doc("Frees the memory used by grid", [], ["streams allocated by init_grid", "robotModel allocated by init_robotModel", "data allocated by init_gridData"], None)
        self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
        self.gen_add_code_line("__host__")
        self.gen_add_code_line("void close_grid(cudaStream_t *streams, robotModel<T> *d_robotModel, gridData<T, KIND> *hd_data){", True)
        self.gen_add_code_lines(["gpuErrchk(cudaFree(d_robotModel));", \
                                 "gpuErrchk(cudaFree(hd_data->d_q_qd_u)); gpuErrchk(cudaFree(hd_data->d_q_qd)); gpuErrchk(cudaFree(hd_data->d_q));", \
                                 "gpuErrchk(cudaFree(hd_data->d_f_ext)); free(hd_data->h_f_ext);", \
                                 "gpuErrchk(cudaFree(hd_data->d_c)); gpuErrchk(cudaFree(hd_data->d_Minv)); gpuErrchk(cudaFree(hd_data->d_qdd)); gpuErrchk(cudaFree(hd_data->d_M));", \
                                 "gpuErrchk(cudaFree(hd_data->d_dc_du)); gpuErrchk(cudaFree(hd_data->d_df_du));", \
                                 "gpuErrchk(cudaFree(hd_data->d_dtau_dfext)); gpuErrchk(cudaFree(hd_data->d_dqdd_dfext));",
                                 "free(hd_data->h_dtau_dfext); free(hd_data->h_dqdd_dfext);",
                                 "gpuErrchk(cudaFree(hd_data->d_did_du_dfext)); free(hd_data->h_did_du_dfext);"]
                                 + [
                                 "gpuErrchk(cudaFree(hd_data->d_eePos)); gpuErrchk(cudaFree(hd_data->d_deePos)); gpuErrchk(cudaFree(hd_data->d_d2eePos));", \
                                 # Phase 3a/b/c/e: end the L2 persisting window opened at init.
                                 "gpuErrchk(grid_end_l2_persisting(0));", \
                                 "gpuErrchk(cudaFree(hd_data->d_workspace));", \
                                # idsva_so - d2tau_dq, d2tau_dqd, d2tau_dvdq, dM_dq
                                 "gpuErrchk(cudaFree(hd_data->d_idsva_so));", \
                                 # fdsva_so - d2fd_dq2, d2fd_cross, d2fd_dqd2, d2fd_dtaudq
                                 "gpuErrchk(cudaFree(hd_data->d_df2));", \
                                 "free(hd_data->h_idsva_so); free(hd_data->h_df2);", \
                                 "free(hd_data->h_q_qd_u); free(hd_data->h_q_qd); free(hd_data->h_q);", \
                                 "free(hd_data->h_c); free(hd_data->h_Minv); free(hd_data->h_qdd); free(hd_data->h_M);", \
                                 "free(hd_data->h_dc_du); free(hd_data->h_df_du);",\
                                 "free(hd_data->h_eePos); free(hd_data->h_deePos); free(hd_data->h_d2eePos);", \
                                 "gpuErrchk(cudaFree(hd_data->d_x_kp1)); gpuErrchk(cudaFree(hd_data->d_dAB));", \
                                 "free(hd_data->h_x_kp1); free(hd_data->h_dAB);", \
                                 "for(int i=0; i<" + str(MAX_STREAMS) + "; i++){gpuErrchk(cudaStreamDestroy(streams[i]));} free(streams);"])
        self.gen_add_end_function()
        
    def gen_combination_functions(self, algorithms, fixed_target_name = ""):
        kinematics_suffix = "" if fixed_target_name == "" else "_" + fixed_target_name

        def has_all(names):
            return all(name in algorithms for name in names)

        def emit_dynamics_combo(name, description, calls):
            self.gen_add_func_doc(description, [], [], None)
            self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
            self.gen_add_code_line("__host__")
            self.gen_add_code_line("void " + name + "(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps,", False)
            self.gen_add_code_line("                   const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {", True)
            self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS, \"" + name + " requires all-data or dynamics gridData\");")
            self.gen_add_code_lines(calls)
            self.gen_add_end_function()

        def emit_kinematics_combo(name, description, calls):
            self.gen_add_func_doc(description, [], [], None)
            self.gen_add_code_line("template <typename T, gridDataKind KIND = GRID_DATA_ALL>")
            self.gen_add_code_line("__host__")
            self.gen_add_code_line("void " + name + "(gridData<T, KIND> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps,", False)
            self.gen_add_code_line("                     const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams) {", True)
            self.gen_add_code_line("static_assert(KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS, \"" + name + " requires all-data or kinematics gridData\");")
            self.gen_add_code_lines(calls)
            self.gen_add_end_function()

        id_call = "inverse_dynamics<T,false,false,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"
        minv_call = "direct_minv<T,false,KIND>(hd_data,d_robotModel,num_timesteps,block_dimms,thread_dimms,streams);"
        fd_call = "forward_dynamics<T,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"
        id_du_call = "inverse_dynamics_gradient<T,false,false,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"
        fd_du_call = "forward_dynamics_gradient<T,false,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"
        aba_call = "aba<T,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"
        crba_call = "crba<T,false,KIND>(hd_data,d_robotModel,gravity,num_timesteps,block_dimms,thread_dimms,streams);"

        if has_all(("id", "minv", "fd")):
            core_calls = [id_call, minv_call, fd_call]
            emit_dynamics_combo("dynamics_core", "Run inverse dynamics, Minv, and forward dynamics in sequence", core_calls)
            emit_dynamics_combo("id_minv_fd", "Run inverse dynamics, Minv, and forward dynamics in sequence", core_calls)

        if has_all(("id", "id_du")):
            emit_dynamics_combo("id_and_id_gradient", "Run inverse dynamics and its first derivative in sequence", [id_call, id_du_call])

        if has_all(("fd", "fd_du")):
            emit_dynamics_combo("fd_and_fd_gradient", "Run forward dynamics and its first derivative in sequence", [fd_call, fd_du_call])

        if has_all(("id_du", "fd_du")):
            emit_dynamics_combo("dynamics_gradients", "Run inverse and forward dynamics gradients in sequence", [id_du_call, fd_du_call])

        if has_all(("id", "minv", "fd", "id_du", "fd_du")):
            calls = [id_call, minv_call, fd_call, id_du_call, fd_du_call]
            if "aba" in algorithms:
                calls.append(aba_call)
            if "crba" in algorithms:
                calls.append(crba_call)
            emit_dynamics_combo("all_dynamics", "Run all generated non-second-order dynamics wrappers in sequence", calls)
            emit_dynamics_combo("dynamics_only", "Run all generated non-second-order dynamics wrappers in sequence", calls)

        kinematics_calls = []
        if "ee_pose" in algorithms:
            kinematics_calls.append("end_effector_pose" + kinematics_suffix + "<T,false,KIND>(hd_data,d_robotModel,num_timesteps,block_dimms,thread_dimms,streams);")
        if "ee_pose_gradient" in algorithms:
            kinematics_calls.append("end_effector_pose_gradient" + kinematics_suffix + "<T,false,KIND>(hd_data,d_robotModel,num_timesteps,block_dimms,thread_dimms,streams);")
        if "ee_pose_hessian" in algorithms:
            kinematics_calls.append("end_effector_pose_gradient_hessian" + kinematics_suffix + "<T,false,KIND>(hd_data,d_robotModel,num_timesteps,block_dimms,thread_dimms,streams);")
        if kinematics_calls:
            emit_kinematics_combo("kinematics_only", "Run all generated kinematics wrappers in sequence", kinematics_calls)

    def gen_add_gpu_err(self):
        # add the GPU error check code
        self.gen_add_func_doc("Check for runtime errors using the CUDA API", \
                ["Adapted from https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api"], \
                [],None)
        self.gen_add_code_line("__host__")
        # `inline` is required so the per-algo TU split (multiple .o files all
        # including grid.cuh) doesn't trip ODR multiple-definition errors at link.
        self.gen_add_code_line("inline void gpuAssert(cudaError_t code, const char *file, const int line, bool abort=true){", True)
        self.gen_add_code_line("if (code != cudaSuccess){", True)
        # note that below we need to escape the \n and "" to get it to print to a string or file correctly
        self.gen_add_code_line("fprintf(stderr,\"GPUassert: %s %s %d\\n\", cudaGetErrorString(code), file, line);")
        self.gen_add_code_line("if (abort){cudaDeviceReset(); exit(code);}")
        self.gen_add_end_control_flow()
        self.gen_add_end_control_flow() # end of function but don't want spacing
        self.gen_add_code_line("#define gpuErrchk(err) {gpuAssert(err, __FILE__, __LINE__);}")
        # gpuErrchkKernel catches BOTH (a) synchronous launch-time errors via
        # cudaPeekAtLastError — e.g. cudaErrorLaunchOutOfResources (code 701)
        # when the kernel asks for more registers than the SM can give — and
        # (b) asynchronous execution-time errors via cudaDeviceSynchronize.
        # Use this after every <<<>>> kernel launch. Plain cudaDeviceSynchronize
        # alone does NOT propagate launch-time errors: a launch can fail before
        # work is queued, leaving the stream empty, so sync returns success and
        # the next call clears the error. That's why the overnight bench was
        # silently reporting failed launches as ~2us "compute time."
        self.gen_add_code_line("#define gpuErrchkKernel() {gpuErrchk(cudaPeekAtLastError()); gpuErrchk(cudaDeviceSynchronize());}")
        self.gen_add_code_line("")

        # also add printMat for debug if requested
        if self.gen_print_mat:
            self.gen_add_code_line("template <typename T, int M, int N>")
            self.gen_add_code_line("__host__ __device__")
            self.gen_add_code_line("void printMat(T *A, int lda){", True)
            self.gen_add_code_line("for(int i=0; i<M; i++){", True)
            self.gen_add_code_line("for(int j=0; j<N; j++){printf(\"%.4f \",A[i + lda*j]);}")
            self.gen_add_code_line("printf(\"\\n\");")
            self.gen_add_end_control_flow()
            self.gen_add_end_function()
            self.gen_add_code_line("template <typename T, int M, int N>")
            self.gen_add_code_line("__host__ __device__")
            self.gen_add_code_line("void printMat(const T *A, int lda){", True)
            self.gen_add_code_line("for(int i=0; i<M; i++){", True)
            self.gen_add_code_line("for(int j=0; j<N; j++){printf(\"%.4f \",A[i + lda*j]);}")
            self.gen_add_code_line("printf(\"\\n\");")
            self.gen_add_end_control_flow()
            self.gen_add_end_function()

    def gen_centroidal_quickwins(self, algorithms):
        """Emit the G2 centroidal quick-win families (R1-R3). Each is gated on
        the grid:: deps it composes being present; a missing dep emits a comment
        instead of an undefined call (mirrors gen_grid_plant's gating). Mimic
        robots are skipped for the kinematics-domain centroidal families (the
        per-body Jacobian fold isn't mimic-reduced yet) — gravity/bias still
        emit since they reuse the mimic-aware RNEA inner."""
        # R1 generalized_gravity / nonlinear_effects: RNEA bias wrappers.
        if "id" in algorithms:
            self.gen_id_bias(gravity_only=True)
            self.gen_id_bias(gravity_only=False)
        else:
            self.gen_add_code_line("// [centroidal] generalized_gravity/nonlinear_effects skipped: require 'id' (grid::inverse_dynamics_inner).")
        # R3/R2/energy: kinematics-domain centroidal families. Need homogeneous
        # transforms (always present when ee_pose is generated).
        kin_ok = ("ee_pose" in algorithms) and not self.robot_has_mimic_joints()
        if kin_ok:
            self.gen_centroidal_inner()
            self.gen_com()
            self.gen_ccrba()
            self.gen_energy()
        elif "ee_pose" not in algorithms:
            self.gen_add_code_line("// [centroidal] com/ccrba/energy skipped: require 'ee_pose' (homogeneous-transform world-frame machinery).")
        else:
            self.gen_add_code_line("// [centroidal] com/ccrba/energy skipped: mimic robots' per-body Jacobian fold is not yet mimic-reduced.")

    # finally generate all of the code
    def gen_all_code(self, include_base_inertia = False, include_homogenous_transforms = False, fixed_target_name = "", output_path = None,
                     codegen_profile = "all", algorithm_list = None, enable_floating_second_order = True,
                     enable_idsva_so_world_frame = None):
        # Default-pick the SO variant that wins per the 2026-05 perf sweep
        # (see test/benchmarks/benchmark_multi_version_sm120_5090_full.md
        # § IDSVA_SO_BODY_FRAME vs IDSVA_SO_WORLD_FRAME):
        #
        #   Robot          body-frame µs   world-frame µs   winner   margin
        #   iiwa14_fixed       26.5             805         body      30.3×
        #   go2_fixed          36.0            1338         body      37.1×
        #   g1_fixed         1301              5804         body       4.5×
        #   iiwa14_floating  2642              1652         world      1.6×
        #   go2_floating     3951              2830         world      1.4×
        #   g1_floating     28222              8451         world      3.3×
        #
        # body-frame multi-pass amortizes well for fixed-base; world-frame's
        # single-pass + no gravity shim wins floating-base. Callers can
        # override with enable_idsva_so_world_frame=True/False to force a
        # specific variant (the bench harness exercises both for comparison).
        if enable_idsva_so_world_frame is None:
            enable_idsva_so_world_frame = self.robot.floating_base
        self.include_fixed_kinematic_targets = fixed_target_name != ""
        algorithms = self._normalize_codegen_algorithms(codegen_profile, algorithm_list)
        self.generated_algorithms = algorithms
        # G0 footgun guard: mimic robots' GRADIENT algorithms are not yet folded
        # (the reduced alpha-scaled gradient assembly is deferred to T3-finisher).
        # Until then refuse to emit a gradient algo for a mimic robot rather than
        # silently writing ZEROS. Non-gradient mimic codegen (id/fd/aba/crba/minv/
        # ee_pose/integrator value) is unaffected and still emits normally.
        if self.robot_has_mimic_joints():
            # T3-finisher: id_du + fd_du mimic gradients have landed (dense
            # serial reduced-space fold). FIXED-BASE (P3) and now FLOATING-BASE
            # (B1: the floating root's 6-DoF motion subspace is folded via a
            # per-root-DoF loop in _gen_id_du_mimic_inner, mirroring the numpy
            # reference rnea_grad_fpass_dq's `for ii in range(len(idx))`). Both
            # are removed from the refusal set. B2-ee: ee_pose_gradient +
            # ee_pose_hessian mimic folds have landed for FIXED-BASE (alpha-
            # weighted geometric-Jacobian column / world-frame generator
            # accumulate; see _eepose_gradient_hessian.py Step 3b / Step 2),
            # removed from the fixed-base refusal set. B2-ee FLOATING: the
            # floating root contributes 6 INDEPENDENT velocity slots (vi 0..5),
            # so it decomposes into 6 singleton single-column geometric-Jacobian /
            # world-generator fills — never a shared-v-slot mimic group. The mimic
            # alpha fold (Step 3b grad / Step 2 hess) operates orthogonally on the
            # 1-DoF mimic joints' shared slots, so floating + mimic compose with no
            # separate 6-DoF root fold. ee_pose_gradient/ee_pose_hessian are now
            # supported for FLOATING-base mimic robots too (removed below).
            # B2-SO: FIXED-base mimic idsva_so/fdsva_so un-refused — the body-frame
            # inner runs the per-body INTERNAL NUM_BODIES-coordinate sweep into a
            # 4*NB^3 slab and alpha-folds to the reduced 4*NV^3 public output (see
            # _idsva_so.py gen_idsva_so_body_frame_inner is_mimic path; also fixed a
            # shared matmul %NUM_JOINTS->%NUM_BODIES block-wrap bug). Floating-base
            # mimic SO stays refused (added below). f_ext_grad remains refused for
            # both bases (P4 pending).
            # B-integrator-grad (FLAG for main reconcile — additive ungate, FIXED
            # base only): the integrator gradients (integrator_gradient +
            # integrator_with_gradient) are now supported for FIXED-base mimic
            # robots, ALL integrator types (Euler / SI-Euler / Midpoint / RK3 /
            # RK4). They COMPOSE the mimic-reduced FD gradient (s_df_du / s_Minv) +
            # the SE(3) dIntegrate blocks and assemble dAB entirely in reduced
            # NV space (n = get_num_vel()), exactly as RBDReference.integrator_grad
            # composes forward_dynamics_grad / minv / dIntegrate (all mimic-aware).
            # The mimic fix was sizing s_vaf at 18*NB (NB>NV for fixed mimic) so
            # the composed FD-grad inner's body-indexed writes don't overflow into
            # s_Minv/s_qdd (mirrors fd_du's kernel sizing). fr3-fixed validates
            # bit-exact (norm_rel ~0 even on RK4 oracle magnitudes ~5e7).
            #
            # FLOATING-base mimic stays REFUSED (added below): the floating
            # SINGLE-stage (Euler / SI-Euler) gradient is correct, but the floating
            # MULTI-stage (Midpoint / RK3 / RK4) gradient is WRONG for mimic only
            # (norm_rel 0.32 / 0.50 / 3.5 on energetic samples) — the bug lives at
            # the floating + multi-stage + mimic intersection (fixed multi-stage
            # mimic is exact; floating multi-stage NON-mimic, e.g. go2, is exact;
            # floating single-stage mimic is float32-clean). Since codegen emits all
            # integrator types into one header, refusing avoids shipping a silently
            # wrong floating-mimic RK gradient. Deferred to a follow-up.
            _MIMIC_GRADIENT_ALGORITHMS = {
                "f_ext_grad",
            }
            if self.robot.floating_base:
                _MIMIC_GRADIENT_ALGORITHMS |= {
                    "integrator_gradient", "integrator_with_gradient",
                }
            # B2-SO FLOATING (FLAG for main reconcile — additive ungate): floating-base
            # mimic SECOND-ORDER is now supported via the WORLD-frame inner (the production
            # floating SO path). The world inner runs the triple ancestor walk in per-column
            # INTERNAL coordinates (n_int = total S-column count; the floating root's 6 DoF
            # get 6 distinct internal slots) into a 4*n_int^3 slab, then alpha-folds each axis
            # to the reduced 4*NV^3 public output — exactly RBDReference.idsva_so_world_frame's
            # has_mimic path. The per-root-DoF treatment is emergent from the per-column
            # internal slotting (the root's columns fold identity, alpha=1), so no separate
            # 6-DoF root fold is needed. fdsva_so composes the world inner on floating, so it
            # is supported too. (Floating mimic SO previously stayed refused here.)
            requested_gradients = sorted(algorithms & _MIMIC_GRADIENT_ALGORITHMS)
            if requested_gradients:
                raise NotImplementedError(
                    "mimic gradients not yet supported — deferred to T3-finisher. "
                    "Robot has mimic joints but the requested codegen selection includes "
                    "gradient algorithm(s) " + ", ".join(requested_gradients) + ". "
                    "These would emit silently-zeroed output (no valid mimic-reduced "
                    "gradient exists yet). Re-run with a non-gradient profile/algorithm "
                    "list (e.g. 'dynamics-core', or id/fd/aba/crba/minv/ee_pose/integrator) "
                    "to codegen this robot."
                )
        self.generate_id_du = "id_du" in algorithms
        self.generate_fd_du = "fd_du" in algorithms
        self.generate_ee_pose_hessian = "ee_pose_hessian" in algorithms
        self.enable_floating_second_order = enable_floating_second_order
        allow_second_order = (not self.robot.floating_base) or enable_floating_second_order
        self.generate_idsva_so_body_frame = ("idsva_so_body_frame" in algorithms) and allow_second_order
        self.generate_fdsva_so = ("fdsva_so" in algorithms) and allow_second_order
        self.generate_idsva_so_world_frame = bool(enable_idsva_so_world_frame) and self.generate_idsva_so_body_frame
        # fdsva_so on floating-base's body calls idsva_so_world_frame_inner<T>(...)
        # unconditionally — without world_frame emission we'd produce a header
        # that fails at link time. Catch the misconfiguration early so the
        # error points at the cause, not a downstream nvcc undefined-symbol.
        if self.generate_fdsva_so and self.robot.floating_base and not self.generate_idsva_so_world_frame:
            raise ValueError(
                "Cannot emit fdsva_so on floating-base without idsva_so_world_frame. "
                "fdsva_so's floating-base body calls idsva_so_world_frame_inner. "
                "Either keep enable_idsva_so_world_frame at its default (None → "
                "True for floating-base) or remove fdsva_so from algorithms."
            )
        include_any_kinematics = any(name in algorithms for name in ("ee_pose", "ee_pose_gradient", "ee_pose_hessian"))
        include_homogenous_transforms = include_homogenous_transforms or include_any_kinematics
        # first generate the file info
        file_notes = [ "Interface is:", \
            "    __host__   robotModel<T> *d_robotModel = init_robotModel<T>()", \
            "    __host__   cudaStream_t streams = init_grid<T>()", \
            "    __host__   gridData<T> *hd_ata = init_gridData<T,NUM_TIMESTEPS>();"
            "    __host__   close_grid<T>(cudaStream_t *streams, robotModel<T> *d_robotModel, gridData<T> *hd_data)", \
            "",\
            "    __device__ inverse_dynamics_inner<T>(T *s_c,  T *s_vaf, const T *s_q, const T *s_qd, const T *s_qdd, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ inverse_dynamics_inner<T>(T *s_c,  T *s_vaf, const T *s_q, const T *s_qd, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ inverse_dynamics_device<T>(T *s_c, const T *s_q, const T *s_qd, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __device__ inverse_dynamics_device<T>(T *s_c, const T *s_q, const T *s_qd, const T *s_qdd, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __global__ inverse_dynamics_kernel<T>(T *d_c, const T *d_q_qd, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __global__ inverse_dynamics_kernel<T>(T *d_c, const T *d_q_qd, const T *d_qdd, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   inverse_dynamics<T,USE_QDD_FLAG=false,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ inverse_dynamics_inner_vaf<T>(T *s_vaf, const T *s_q, const T *s_qd, const T *s_qdd, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ inverse_dynamics_inner_vaf<T>(T *s_vaf, const T *s_q, const T *s_qd, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ inverse_dynamics_vaf_device<T>(T *s_vaf, const T *s_q, const T *s_qd, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __device__ inverse_dynamics_vaf_device<T>(T *s_vaf, const T *s_q, const T *s_qd, const T *s_qdd, const robotModel<T> *d_robotModel, const T gravity)", \
            "",\
            "    __device__ direct_minv_inner<T>(T *s_Minv, T *s_F, const T *s_q, T *s_XImats, int *s_topology_helpers, T *s_temp)",\
            "    __device__ direct_minv_device<T>(T *s_Minv, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ direct_minv_Kernel<T>(T *d_Minv, unsigned char *d_workspace, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   direct_minv<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ forward_dynamics_inner<T>(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, T *s_minv_F, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ forward_dynamics_device<T, RESOURCE_TIER=TIER_SHARED>(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, const robotModel<T> *d_robotModel, const T gravity, T *d_workspace = nullptr)", \
            "    __global__ forward_dynamics_kernel<T>(T *d_qdd, unsigned char *d_workspace, const T *d_q_qd_u, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   forward_dynamics<T>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ inverse_dynamics_gradient_inner<T>(T *s_dc_du, const T *s_q, const T *s_qd, const T *s_vaf, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ inverse_dynamics_gradient_device<T>(T *s_dc_du, const T *s_q, const T *s_qd, const T *robotModel<T> *d_robotModel, const T gravity)", \
            "    __device__ inverse_dynamics_gradient_device<T>(T *s_dc_du, const T *s_q, const T *s_qd, const T *s_qdd, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __global__ inverse_dynamics_gradient_kernel<T>(T *d_dc_du, const T *d_q_qd, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __global__ inverse_dynamics_gradient_kernel<T>(T *d_dc_du, const T *d_q_qd, const T *d_qdd, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   inverse_dynamics_gradient<T,USE_QDD_FLAG=false,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ forward_dynamics_gradient_device<T>(T *s_df_du, const T *s_q, const T *s_qd, const T *s_u, const robotModel<T> *d_robotModel, const T gravity)",\
            "    __device__ forward_dynamics_gradient_device<T>(T *s_df_du, const T *s_q, const T *s_qd, const T *s_qdd, const T *s_Minv, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __global__ forward_dynamics_gradient_kernel<T>(T *d_df_du, const T *d_q_qd_u, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __global__ forward_dynamics_gradient_kernel<T>(T *d_df_du, const T *d_q_qd, const T *d_qdd, const T *d_Minv, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   forward_dynamics_gradient<T,USE_QDD_MINV_FLAG=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ end_effector_pose_inner<T,TEMP_IN_SMEM=true>(T *s_eePos, const T *s_q, const T *s_Xhom, int *s_topology_helpers, T *s_temp, T *d_workspace, unsigned char *s_linalg_smem)", \
            "    __device__ end_effector_pose_device<T>(T *s_eePos, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ end_effector_pose_kernel<T>(T *d_eePos, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   end_effector_pose<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ end_effector_pose_gradient_inner<T>(T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, int *s_topology_helpers, T *s_temp)", \
            "    __device__ end_effector_pose_gradient_device<T>(T *s_deePos, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ end_effector_pose_gradient_kernel<T>(T *d_deePos, unsigned char *d_workspace, const T *d_q, const int stride_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   end_effector_pose_gradient<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ end_effector_pose_gradient_hessian_inner<T>(T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, int *s_topology_helpers, T *s_temp)", \
            "    __device__ end_effector_pose_gradient_hessian_device<T>(T *s_deePos, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ end_effector_pose_gradient_hessian_kernel<T>(T *d_deePos, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   end_effector_pose_gradient_hessian<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ idsva_so_body_frame_inner(T *s_idsva_so, const T *s_q, const T *s_qd, T *s_qdd, T *s_XImats, T *s_mem, const T gravity)",\
            "    __global__ idsva_so_body_frame_kernel(T *d_idsva_so, const T *d_q_qd_u, const int stride_q_qd_u, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   idsva_so_body_frame_host<T>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ fdsva_so_contract(T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, T *s_q, T *s_qd, const T *s_qdd, const T *s_tau, T *s_XImats, T *s_temp, const T gravity)",\
            "    __device__ fdsva_so_device(T *s_df2, T *s_df_du, const T *s_q, const T *s_qd, const T *s_u, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __global__ fdsva_so_kernel(T *d_df2, const T *d_q_qd_qdd_tau, const int stride_q_qd_qdd, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
            "    __host__   fdsva_so<T>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const T gravity, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "","Suggested Type T is float",\
            "","Additional helper functions and ALGORITHM_inner functions which take in __shared__ memory temp variables exist -- see function descriptions in the file",\
            "","By default device and kernels need to be launched with dynamic shared mem of size <FUNC_CODE>_DYNAMIC_SHARED_MEM_COUNT where <FUNC_CODE> = [ID, MINV, FD, ID_DU, FD_DU]"]
        file_notes += ["", "Codegen profile: " + str(codegen_profile), "Generated algorithms: " + ", ".join(sorted(algorithms))]
        if self.include_fixed_kinematic_targets:
            file_notes += ["", "Additional EEPose Functions Included for Fixed Kinematic Target: " + fixed_target_name,""]
        self.gen_add_func_doc("This instance of grid.cuh is optimized for the urdf: " + self.robot.name,file_notes)
        # then all of the includes (and namespaces and defines)
        self.gen_add_includes()
        # then add the gpu error macro
        self.gen_add_gpu_err()
        # File-scope preprocessor mirrors of the second-order codegen gates.
        # We also emit `const int GRID_GENERATES_*_SO` inside `namespace grid`
        # below (for C++ runtime / test consumers), but `#if X` requires a
        # preprocessor macro — a namespaced const reads as the undefined-symbol
        # zero in `#if`, silently turning off any consumer that gates on it
        # (notably the bench's timeGRiD_{single,batch}.cu measure_* blocks).
        # Different names from the namespaced const avoid macro/const collision.
        self.gen_add_code_line(
            "#define GRID_HAS_IDSVA_SO_BODY_FRAME " + str(int(getattr(self, "generate_idsva_so_body_frame", True)))
        )
        self.gen_add_code_line(
            "#define GRID_HAS_FDSVA_SO " + str(int(getattr(self, "generate_fdsva_so", True)))
        )
        self.gen_add_code_line(
            "#define GRID_HAS_IDSVA_SO_WORLD_FRAME " + str(int(getattr(self, "generate_idsva_so_world_frame", False)))
        )
        # GRID_HAS_IDSVA_SO gates the dispatching `grid::idsva_so` host wrapper.
        # Emitted whenever the chosen variant for this robot's base type is
        # available: body_frame for fixed-base (always), world_frame for
        # floating-base (opt-in via enable_idsva_so_world_frame).
        has_idsva_so = getattr(self, "generate_idsva_so_body_frame", True) and (
            (not self.robot.floating_base) or getattr(self, "generate_idsva_so_world_frame", False)
        )
        self.gen_add_code_line("#define GRID_HAS_IDSVA_SO " + str(int(has_idsva_so)))
        # Integrator availability gates. Both value and gradient kernels are
        # emitted whenever requested, for fixed- and floating-base. (Floating
        # gradient is currently Euler-only and inherits the upstream
        # forward_dynamics_gradient dqdd/dqd bug; see _normalize_codegen_algorithms.)
        # Consumers must #if-guard gradient calls so a build that omits
        # integrator_gradient still compiles.
        self.gen_add_code_line(
            "#define GRID_HAS_INTEGRATOR " + str(int("integrator" in algorithms)))
        self.gen_add_code_line(
            "#define GRID_HAS_INTEGRATOR_GRADIENT " + str(int("integrator_gradient" in algorithms)))
        self.gen_add_code_line("")
        # then open our namespace
        self.gen_add_func_doc("All functions are kept in this namespace")
        self.gen_add_code_line("namespace " + self.file_namespace + " {", True)
        self.gen_add_shared_memory_helpers()
        # then generate any constants and other helpers
        self.gen_add_constants_helpers(include_base_inertia, include_homogenous_transforms)
        # then the linear algebra related helpers
        # Emit GLASS before spatial algebra because dot_prod is a compatibility
        # shim over glass::dot_strided.
        self.gen_grid_linalg_backend_helpers()
        # then the spatial algebra related helpers
        self.gen_spatial_algebra_helpers()
        self.gen_crm()
        self.gen_crm_mul()
        self.gen_invert_matrix()
        self.gen_matmul()
        self.gen_matmul_trans() 
        self.gen_outer_product()
        # then generate the robot specific transformation and inertia matricies
        self.gen_init_topology_helpers()
        self.gen_init_XImats(include_base_inertia, include_homogenous_transforms)
        self.gen_init_robotModel()
        self.gen_init_gridData()
        self.gen_joint_limits_size()
        self.gen_init_joint_limits()
        self.gen_load_update_XImats_helpers()
        if include_homogenous_transforms and include_any_kinematics:
            self.gen_load_update_XmatsHom_helpers(include_base_inertia)
            if "ee_pose_gradient" in algorithms or "ee_pose_hessian" in algorithms:
                self.gen_load_update_XmatsHom_helpers(include_base_inertia,include_gradients = True)
            if "ee_pose_hessian" in algorithms:
                self.gen_load_update_XmatsHom_helpers(include_base_inertia,include_gradients = True, include_hessians = True)
        # then generate kinematic algorithms.
        # SE(3) Lie-group helpers (grid_integrate_floating_q, grid_so3_*, grid_quat_*)
        # are needed by the FD-on-Jacobian d2ee inner on floating base. Emit them
        # here too so callers that skip the integrator codegen still get them; track
        # the emission so gen_integrator skips its own emit (avoiding redefinitions).
        self._lie_helpers_emitted = False
        if include_any_kinematics:
            if self.robot.floating_base and "ee_pose_hessian" in algorithms:
                self.gen_lie_group_helpers()
                self._lie_helpers_emitted = True
            self.gen_eepose_and_derivatives(fixed_target_name = fixed_target_name,
                                            include_pose = "ee_pose" in algorithms,
                                            include_gradient = "ee_pose_gradient" in algorithms,
                                            include_hessian = "ee_pose_hessian" in algorithms)
        if self.robot.floating_base and not enable_floating_second_order:
            print('floating-base second order dynamics are still under development')
        # then generate the dynamics algorithms
        if "id" in algorithms:
            self.gen_inverse_dynamics()
        # E1: joint-torque regressor Y (tau = Y . pi). Additive; reuses the RNEA
        # forward sweep emitted by gen_inverse_dynamics (requires "id").
        if "regressor" in algorithms:
            self.gen_inverse_dynamics_regressor()
        if "minv" in algorithms:
            self.gen_direct_minv()
        if "fd" in algorithms:
            self.gen_forward_dynamics()
        # FD parameter gradient dqdd/dpi = -Minv . Y. Additive; composes the
        # regressor (Y), direct_minv (Minv) and inverse_dynamics/forward_dynamics
        # inners, so it requires "id", "minv", "fd" and "regressor" co-emitted.
        if "fd_parameter_gradient" in algorithms:
            self.gen_fd_parameter_gradient()
        if "id_du" in algorithms:
            self.gen_inverse_dynamics_gradient()
        if "fd_du" in algorithms:
            self.gen_forward_dynamics_gradient()
        if "f_ext_grad" in algorithms:
            self.gen_f_ext_gradient()
        if "aba" in algorithms:
            self.gen_aba()
        if "crba" in algorithms:
            self.gen_crba()
        if "integrator" in algorithms:
            self.gen_integrator()
        if ("integrator_gradient" in algorithms) or ("integrator_with_gradient" in algorithms):
            self.gen_integrator_gradient()
        if not self.robot.floating_base or enable_floating_second_order:
            if "idsva_so_body_frame" in algorithms:
                self.gen_idsva_so_body_frame()
                # Optional: emit the world-frame single-pass alternative path alongside
                # the existing emission. Co-exists with `idsva_so_body_frame_kernel`/`idsva_so_body_frame_host`;
                # the new entry point is `idsva_so_world_frame_kernel`/`idsva_so_world_frame_host`.
                if enable_idsva_so_world_frame:
                    self.gen_idsva_so_world_frame()
                # Emit the dispatching `idsva_so` host wrapper. For floating-base
                # robots, requires world_frame to be enabled (it forwards there).
                # For fixed-base, forwards to body_frame.
                if (not self.robot.floating_base) or enable_idsva_so_world_frame:
                    self.gen_idsva_so_dispatcher()
            if "fdsva_so" in algorithms:
                self.gen_fdsva_so()
        # G2 centroidal quick-wins (R1-R3): additive families gated on their
        # grid:: deps. generalized_gravity / nonlinear_effects are RNEA bias
        # wrappers (need `id`); com / ccrba / energy live in the kinematics
        # (homogeneous-transform) domain and reuse the world-transform machinery
        # (need `ee_pose`). All are NEW emitters appended after the existing
        # algorithms, so existing emission is byte-identical.
        self.gen_centroidal_quickwins(algorithms)
        # E2 (additive, opt-in): general-frame geometric Jacobian. Only emitted
        # when the `frame_jacobian` key is explicitly selected, so every existing
        # profile's header is byte-identical. Needs ee_pose's world-transform
        # machinery (pulled in by _normalize_codegen_algorithms).
        if "frame_jacobian" in algorithms and "ee_pose" in algorithms and not self.robot_has_mimic_joints():
            NJ_fj = self.robot.get_num_joints()
            nv_fj = self.robot.get_num_vel()
            n_pos_fj = self.robot.get_num_pos()
            Xhom_size_fj, _, _ = self.gen_get_Xhom_size()
            # arena = s_XmatsHom(Xhom_size) + inner_temp(16*NJ); s_J is a caller param.
            fj_t_count = Xhom_size_fj + (16 * NJ_fj)
            self.gen_add_code_line(
                "template <typename T> __host__ __device__ inline size_t FRAME_JACOBIAN_DYNAMIC_SHARED_MEM_BYTES() "
                "{ return grid_shared_arena_bytes<T>(" + str(fj_t_count) +
                ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }")
            self.gen_frame_jacobian()
            # E2 CUDA parity (opt-in siblings). Jdot/Lambda reuse frame_jacobian_inner.
            if "frame_jacobian_dot" in algorithms:
                # Jdot arena = s_XmatsHom + extras(s_qpert[n_pos] + s_Jp[6nv] + s_Jm[6nv]) + inner_temp(16*NJ).
                fjd_t_count = Xhom_size_fj + n_pos_fj + (2 * 6 * nv_fj) + (16 * NJ_fj)
                self.gen_add_code_line(
                    "template <typename T> __host__ __device__ inline size_t FRAME_JACOBIAN_DOT_DYNAMIC_SHARED_MEM_BYTES() "
                    "{ return grid_shared_arena_bytes<T>(" + str(fjd_t_count) +
                    ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }")
                # Floating-base Jdot integrates q on the SE(3) group; emit the Lie
                # helpers if no other kinematics path already did.
                if self.robot.floating_base and not getattr(self, "_lie_helpers_emitted", False):
                    self.gen_lie_group_helpers()
                    self._lie_helpers_emitted = True
                self.gen_frame_jacobian_dot()
            if "osc_inertia" in algorithms:
                # Lambda is SELF-CONTAINED: it composes Minv on device via
                # direct_minv_inner, so the arena carries BOTH transform families
                # (spatial s_XImats for minv + homogeneous s_XmatsHom for J) plus
                # the minv buffers (s_Minv + the spilled F-region passed as
                # d_workspace) and the J*Minv*J^T compose scratch.
                # arena = s_XImats(XI) + extras + s_temp(max(no_F, 16*NJ)) where
                # extras = s_XmatsHom + s_Minv + s_F + s_Jfj + s_MJt + s_task + s_taskinv.
                osc_XI_size = self.gen_get_XI_size(False, False)
                osc_noF = self.gen_direct_minv_inner_no_F_size()
                osc_F = self.gen_direct_minv_inner_F_size()
                osc_temp = max(osc_noF, 16 * NJ_fj)
                osc_t_count = (osc_XI_size + Xhom_size_fj + (nv_fj * nv_fj) + osc_F
                               + (6 * nv_fj) + (nv_fj * 6) + 36 + 36 + osc_temp)
                self.gen_add_code_line(
                    "template <typename T> __host__ __device__ inline size_t OSC_INERTIA_DYNAMIC_SHARED_MEM_BYTES() "
                    "{ return grid_shared_arena_bytes<T>(" + str(osc_t_count) +
                    ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }")
                self.gen_osc_inertia()
        self.gen_combination_functions(algorithms, fixed_target_name)
        # then finally the master init and close the namespace
        self.gen_init_close_grid()
        self.gen_add_end_control_flow()
        # T6: emit the sibling `grid_plant` namespace (cost/constraint/plant-step
        # primitives composed over the grid:: surface). Additive: this runs AFTER
        # the grid namespace closes and makes ZERO edits to any grid:: emit path.
        self.gen_grid_plant(algorithms)
        # then output to a file
        if output_path is None:
            output_path = self.file_namespace + ".cuh"
        file = open(output_path, "w")
        file.write(self.code_str)
        file.close()
