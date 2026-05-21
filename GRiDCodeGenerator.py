import os
import numpy as np

class GRiDCodeGenerator:
    # first import helpers to write code generation, spatial algebra, and opology helpers (parent, child, Sind, XImats) and the robotModel object wrapepr
    from .helpers import gen_add_code_line, gen_add_code_lines, gen_add_end_control_flow, gen_add_end_function, \
                         gen_add_func_doc, gen_add_serial_ops, gen_add_parallel_loop, gen_add_sync, gen_var_in_list, \
                         gen_var_not_in_list, gen_add_multi_threaded_select, gen_kernel_load_inputs, gen_kernel_save_result, \
                         gen_kernel_load_inputs_single_timing, gen_kernel_save_result_single_timing, gen_anti_licm_input_reload, gen_anti_licm_output_write, \
                         gen_static_array_ind_2d, gen_static_array_ind_3d, gen_add_debug_print_code_lines, \
                         gen_mx_func_call_for_cpp, gen_add_shared_memory_helpers, gen_declare_shared_arena, \
                         gen_shared_arena_t_count, gen_spatial_algebra_helpers, \
                         gen_get_XI_size, gen_init_XImats, gen_load_update_XImats_helpers_temp_mem_size, gen_load_update_XImats_helpers_function_call, \
                         gen_XImats_helpers_temp_shared_memory_code, gen_load_update_XImats_helpers, gen_topology_helpers_size, \
                         gen_get_Xhom_size, gen_load_update_XmatsHom_helpers, gen_load_update_XmatsHom_helpers_function_call, gen_XmatsHom_helpers_temp_shared_memory_code, \
                         gen_topology_sparsity_helpers_python, gen_init_topology_helpers, gen_topology_helpers_pointers_for_cpp, \
                         gen_topology_S_sign_for_cpp, gen_insert_helpers_function_call, gen_insert_helpers_func_def_params, gen_init_robotModel, gen_joint_limits_size, gen_init_joint_limits, \
                         gen_grid_linalg_backend_helpers, gen_linalg_smem_setup, gen_invert_matrix, gen_matmul, gen_matmul_trans, gen_crm_mul, gen_crm, gen_outer_product, custom_is_constant

    # then import all of the algorithms
    from .algorithms import gen_inverse_dynamics_inner_temp_mem_size, gen_inverse_dynamics_inner_function_call, \
                            gen_inverse_dynamics_device_temp_mem_size, gen_inverse_dynamics_inner, gen_inverse_dynamics_device, \
                            gen_inverse_dynamics_kernel, gen_inverse_dynamics_host, gen_inverse_dynamics, \
                            gen_direct_minv_inner_temp_mem_size, gen_direct_minv_inner_function_call, gen_direct_minv_inner, \
                            gen_direct_minv_device, gen_direct_minv_kernel, gen_direct_minv_host, gen_direct_minv, \
                            gen_forward_dynamics_inner_temp_mem_size, gen_forward_dynamics_finish_function_call, gen_forward_dynamics_finish, \
                            gen_forward_dynamics_inner_function_call, gen_forward_dynamics_inner, gen_forward_dynamics_device, \
                            gen_forward_dynamics_kernel, gen_forward_dynamics_host, gen_forward_dynamics, \
                            gen_inverse_dynamics_gradient_inner_temp_mem_size, gen_inverse_dynamics_gradient_temp_layout, \
                            gen_inverse_dynamics_gradient_kernel_max_temp_mem_size, \
                            gen_inverse_dynamics_gradient_inner_function_call, gen_inverse_dynamics_gradient_inner, gen_inverse_dynamics_gradient_device, \
                            gen_inverse_dynamics_gradient_kernel, gen_inverse_dynamics_gradient_host, gen_inverse_dynamics_gradient, \
                            gen_forward_dynamics_gradient_inner_temp_mem_size, gen_forward_dynamics_gradient_kernel_max_temp_mem_size, \
                            gen_forward_dynamics_gradient_inner_python, gen_forward_dynamics_gradient_device, gen_forward_dynamics_gradient_kernel, \
                            gen_forward_dynamics_gradient_host, gen_forward_dynamics_gradient, gen_forward_dynamics_gradient_device_function_call, \
                            gen_end_effector_pose_inner_temp_mem_size, gen_end_effector_pose_inner_function_call, gen_end_effector_pose_inner, \
                            gen_end_effector_pose_device_temp_mem_size, gen_end_effector_pose_device, gen_end_effector_pose_kernel, \
                            gen_end_effector_pose_host, gen_end_effector_pose_gradient_inner_temp_mem_size, gen_end_effector_pose_gradient_inner_function_call, \
                            gen_end_effector_pose_gradient_inner, gen_end_effector_pose_gradient_device, gen_end_effector_pose_gradient_kernel, \
                            gen_end_effector_pose_gradient_host, gen_end_effector_pose_gradient_hessian_d2_temp_mem_size, gen_end_effector_pose_gradient_hessian_inner_temp_mem_size, gen_end_effector_pose_gradient_hessian_inner_function_call, \
                            gen_end_effector_pose_gradient_hessian_inner, gen_end_effector_pose_gradient_hessian_device, gen_end_effector_pose_gradient_hessian_kernel, gen_X_single_thread, gen_X_warp, \
                            gen_end_effector_pose_gradient_hessian_host, gen_eepose_and_derivatives, \
                            gen_aba, gen_aba_inner, gen_aba_host, \
                            gen_aba_inner_function_call, gen_aba_kernel, gen_aba_device, gen_aba_inner_temp_mem_size, \
                            gen_crba, gen_crba_inner_temp_mem_size, gen_crba_inner_function_call, gen_crba_inner, gen_crba_device_temp_mem_size, \
                            gen_crba_device, gen_crba_kernel, gen_crba_host, \
                            gen_idsva_so_body_frame_inner_temp_mem_size, gen_idsva_so_body_frame_inner_function_call, idsva_so_needs_reference_order_output_repair, \
                            gen_idsva_so_body_frame_reference_order_output_repair, gen_idsva_so_body_frame_floating_reference_inner, gen_idsva_so_body_frame_public_dvdq_layout_repair, gen_idsva_so_body_frame_inner, gen_idsva_so_body_frame_device_temp_mem_size, \
                            gen_idsva_so_body_frame_device, gen_idsva_so_body_frame_kernel, gen_idsva_so_body_frame_host, gen_idsva_so_body_frame, \
                            gen_idsva_so_world_frame_temp_mem_size, gen_idsva_so_world_frame_inner, \
                            gen_idsva_so_world_frame_inner_function_call, gen_idsva_so_world_frame_kernel, \
                            gen_idsva_so_world_frame_host, gen_idsva_so_world_frame, \
                            gen_idsva_so_device, gen_idsva_so_dispatcher_host, gen_idsva_so_dispatcher, \
                            gen_floating_gravity_d2tau_dq_temp_mem_size, gen_floating_gravity_d2tau_dq_shared_count, \
                            gen_floating_gravity_d2tau_dq_spill_count, gen_floating_gravity_d2tau_dq_lie_inline, \
                            gen_fdsva_so, gen_fdsva_so_inner_temp_mem_size, gen_fdsva_so_fd_gradient_inline_temp_mem_size, gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled, gen_fdsva_so_fd_gradient_inline, gen_fdsva_so_inner_function_call, gen_fdsva_so_inner, gen_fdsva_so_device_temp_mem_size, \
                            gen_fdsva_so_device, gen_fdsva_so_kernel, gen_fdsva_so_host, \
                            gen_integrator_inner_temp_mem_size, gen_integrator_finish_function_call, gen_integrator_finish, \
                            gen_integrator_inner_function_call, gen_integrator_inner, gen_integrator_device, \
                            gen_integrator_kernel, gen_integrator_host, gen_integrator, gen_lie_group_helpers, \
                            gen_integrator_gradient_inner_temp_mem_size, gen_integrator_gradient_dAB_assembly, \
                            gen_integrator_gradient_inner_python, gen_integrator_gradient_multistage, gen_integrator_gradient_device, \
                            gen_integrator_gradient_kernel, gen_integrator_gradient_host, gen_integrator_gradient

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
        self.cuda_shared_mem_type_size_bytes = int(os.environ.get("GRID_CUDA_SHARED_MEM_TYPE_SIZE_BYTES", "4"))

    def _normalize_codegen_algorithms(self, codegen_profile = "all", algorithm_list = None):
        all_algorithms = {
            "id", "minv", "fd", "id_du", "fd_du", "aba", "crba",
            "idsva_so_body_frame", "fdsva_so", "ee_pose", "ee_pose_gradient", "ee_pose_hessian",
            "integrator", "integrator_gradient", "integrator_with_gradient",
        }
        profile_algorithms = {
            "all": all_algorithms,
            "dynamics": {"id", "minv", "fd", "id_du", "fd_du", "aba", "crba", "idsva_so_body_frame", "fdsva_so",
                         "integrator", "integrator_gradient", "integrator_with_gradient"},
            "dynamics-core": {"id", "minv", "fd"},
            "dynamics-gradients": {"id", "minv", "fd", "id_du", "fd_du"},
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
                elif key in all_algorithms:
                    algorithms.add(key)
                else:
                    raise ValueError("Unknown GRiD algorithm selection: " + str(item))

        if "fd_du" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du"})
        if "id_du" in algorithms:
            algorithms.add("id")
        if "aba" in algorithms and self.robot.floating_base:
            algorithms.update({"id", "minv", "fd"})
        if "fdsva_so" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du", "fd_du", "idsva_so_body_frame"})
        if "idsva_so_body_frame" in algorithms:
            algorithms.add("id")
            if self.robot.floating_base:
                algorithms.add("id_du")
        # integrator value needs forward dynamics; gradient needs FD + FD-gradient.
        # The integrator GRADIENT kernels are not yet implemented for
        # floating-base (the value kernels are), so drop them from the set on
        # floating-base robots — this keeps the kernel-attr manifest and
        # gen_all_code consistent (no dangling references to ungenerated
        # kernels). The Python RBDReference.integrator_grad still supports
        # floating-base for CPU reference use.
        if self.robot.floating_base:
            algorithms.discard("integrator_gradient")
            algorithms.discard("integrator_with_gradient")
        if "integrator" in algorithms:
            algorithms.update({"id", "minv", "fd"})
        if "integrator_gradient" in algorithms or "integrator_with_gradient" in algorithms:
            algorithms.update({"id", "minv", "fd", "id_du", "fd_du"})
        return algorithms
    
    # add generic code needs and helpers (includes, memory initialization, constants, kernel settings etc.)
    def gen_add_includes(self, use_thread_group = False):
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
        if use_thread_group:
            self.gen_add_code_line("#include <cooperative_groups.h>")
            self.gen_add_code_line("#include <cooperative_groups/memcpy_async.h>")
        # then any namespaces
        if use_thread_group:
            self.gen_add_code_line("namespace cgrps = cooperative_groups;")
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
        suggested_threads = 32 * int(np.ceil(max_threads_in_comp_loop/32.0))
        # cap to 512 mirrors the constant we emit further down; expose on self so
        # _lin_alg_helpers can pin cuBLASDx's BlockDim<TC,1,1> to the same value.
        self.suggested_threads = min(suggested_threads, 512)
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

        def select_shared_tier(full_t_count, selective_t_count):
            if py_arena_bytes(full_t_count) <= self.cuda_target_shared_mem_bytes:
                return 0
            if py_arena_bytes(selective_t_count) <= self.cuda_target_shared_mem_bytes:
                return 1
            return 2

        id_t_count = 2*n + n + 18*n + n + self.gen_inverse_dynamics_inner_temp_mem_size() + XI_size
        minv_t_count = n + n*n + self.gen_direct_minv_inner_temp_mem_size() + XI_size
        fd_t_count = 3*nv + int(self.robot.floating_base) + nv + self.gen_forward_dynamics_inner_temp_mem_size() + XI_size
        # Integrator: kernel-shared t-count layout is
        #   s_q_qd_u (3nv+fb) + s_qdd (nv) + s_stage_qdd ((max_stages-1)*nv)
        #   + s_stage_point ((max_stages-1)*(2nv+fb)) + s_x_kp1 (2nv+fb)
        #   + s_temp (= FD inner)
        # The "+fb" terms account for the floating-base quaternion (q has 1
        # more element than v).  max_stages = 4 (RK4) — see
        # _integrator._max_stages_in_use().
        _max_stages = 4
        _fb = int(self.robot.floating_base)
        integrator_t_count = ((3*nv + _fb) + nv
                              + (_max_stages - 1) * nv
                              + (_max_stages - 1) * (2 * nv + _fb)
                              + (2 * nv + _fb)
                              + self.gen_forward_dynamics_inner_temp_mem_size() + XI_size)
        # Integrator gradient: kernel-shared t-count layout is
        #   s_q_qd_u (3nv+fb) + s_dAB (2nv*3nv) + s_df_du (nv*2nv) + s_dc_du (nv*2nv) +
        #   s_vaf (18nv) + s_Minv (nv*nv) + s_qdd (nv)
        #   + multi-stage scratch: s_q_orig (nv) + s_qd_orig (nv)
        #     + s_stage_grad_qdd (max_stages*nv) + s_D_qdd_stage (max_stages*nv*3nv)
        #   + s_temp (= FD-grad inner)
        # max_stages = 4 (RK4) — see _integrator._max_stages_in_use().
        # The multi-stage scratch is always allocated even for single-stage IT;
        # cost is small relative to total (~12*nv² for iiwa14 ≈ 588 floats).
        _max_stages = 4
        # +72 for the two 6x6 SE(3) dIntegrate blocks (floating-base gradient;
        # allocated for fixed-base too but unused there).
        integrator_du_t_count = ((3*nv + int(self.robot.floating_base)) + 2*nv*3*nv + 2*(nv*2*nv)
                                 + 18*nv + nv*nv + nv
                                 + 2*nv + _max_stages * nv + _max_stages * nv * 3*nv
                                 + 72
                                 + self.gen_forward_dynamics_gradient_inner_temp_mem_size() + XI_size)
        # The "with x_kp1" variant adds s_x_kp1 (nq+nv = 2nv+fb) on top.
        integrator_du_with_x_kp1_t_count = integrator_du_t_count + 2*nv + int(self.robot.floating_base)
        id_du_temp_layout = self.gen_inverse_dynamics_gradient_temp_layout()
        id_du_temp_count = id_du_temp_layout["full_count"]
        id_du_selective_temp_count = id_du_temp_layout["selective_shared_count"]
        fd_du_temp_count = self.gen_forward_dynamics_gradient_inner_temp_mem_size()
        fd_du_selective_temp_count = max(self.gen_direct_minv_inner_temp_mem_size(), id_du_selective_temp_count)
        id_device_t_count = 18*n + self.gen_inverse_dynamics_inner_temp_mem_size() + XI_size
        minv_device_t_count = self.gen_direct_minv_inner_temp_mem_size() + XI_size
        fd_device_t_count = self.gen_forward_dynamics_inner_temp_mem_size() + XI_size
        id_du_device_t_count = 18*nv + id_du_temp_count + XI_size
        fd_du_device_t_count = (2*nv*nv) + (18*nv) + nv + (nv*nv) + fd_du_temp_count + XI_size
        id_du_t_count_full = (nv + n) + (2*nv*nv) + (18*nv) + nv + id_du_temp_count + XI_size
        fd_du_t_count_full = (3*nv + int(self.robot.floating_base)) + (2*nv*nv) + (18*nv) + nv + (nv*nv) + fd_du_temp_count + XI_size
        id_du_t_count_selective = id_du_t_count_full - id_du_temp_count + id_du_selective_temp_count
        fd_du_t_count_selective = fd_du_t_count_full - fd_du_temp_count + fd_du_selective_temp_count
        id_du_t_count_emergency = id_du_t_count_full - id_du_temp_count
        fd_du_t_count_emergency = fd_du_t_count_full - fd_du_temp_count
        self.id_du_spill_tier = select_shared_tier(id_du_t_count_full, id_du_t_count_selective)
        self.fd_du_spill_tier = select_shared_tier(fd_du_t_count_full, fd_du_t_count_selective)
        self.id_du_use_selective_spill = self.id_du_spill_tier == 1
        self.fd_du_use_selective_spill = self.fd_du_spill_tier == 1
        self.id_du_use_global_temp = self.id_du_spill_tier == 2
        self.fd_du_use_global_temp = self.fd_du_spill_tier == 2
        id_du_t_count = [id_du_t_count_full, id_du_t_count_selective, id_du_t_count_emergency][self.id_du_spill_tier]
        fd_du_t_count = [fd_du_t_count_full, fd_du_t_count_selective, fd_du_t_count_emergency][self.fd_du_spill_tier]
        aba_input_t_count = n + 2*nv
        crba_input_t_count = n + nv
        aba_t_count = nv + aba_input_t_count + 12*NJ + self.gen_aba_inner_temp_mem_size() + XI_size
        crba_t_count = nv*nv + crba_input_t_count + self.gen_crba_inner_temp_mem_size() + XI_size
        ee_t_count = n + 6*self.robot.get_total_leaf_nodes() + self.gen_end_effector_pose_inner_temp_mem_size() + XHom_size
        dee_t_count = n + 6*n*self.robot.get_total_leaf_nodes() + self.gen_end_effector_pose_gradient_inner_temp_mem_size() + XHom_size + dXhom_size
        d2ee_inner_temp_count_full = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size()
        d2ee_inner_temp_count_shared = self.gen_end_effector_pose_gradient_hessian_inner_temp_mem_size(include_d2_temp=False)
        d2ee_workspace_temp_count = self.gen_end_effector_pose_gradient_hessian_d2_temp_mem_size()
        d2ee_full_t_count = n + 6*n*n*self.robot.get_total_leaf_nodes() + 6*n*self.robot.get_total_leaf_nodes() + d2ee_inner_temp_count_full + XHom_size + dXhom_size + d2Xhom_size
        d2ee_spill_t_count = n + d2ee_inner_temp_count_shared + XHom_size + dXhom_size + d2Xhom_size
        d2ee_spill_d2xhom_t_count = n + d2ee_inner_temp_count_shared + XHom_size + dXhom_size
        self.d2ee_spill_tier = 0
        if "ee_pose_hessian" in getattr(self, "generated_algorithms", set()):
            if py_arena_bytes(d2ee_full_t_count) > self.cuda_target_shared_mem_bytes:
                self.d2ee_spill_tier = 1 if py_arena_bytes(d2ee_spill_t_count) <= self.cuda_target_shared_mem_bytes else 2
        self.d2ee_use_workspace_temp = self.d2ee_spill_tier >= 1
        self.d2ee_use_workspace_d2xhom = self.d2ee_spill_tier >= 2
        d2ee_t_count = [d2ee_full_t_count, d2ee_spill_t_count, d2ee_spill_d2xhom_t_count][self.d2ee_spill_tier]
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
        # After the spill decision is final, capture the inner temp count for use
        # downstream (world-frame fallback for fixed-base + FDSVA-SO inner sizing).
        idsva_so_body_frame_inner_temp_count = self.gen_idsva_so_body_frame_inner_temp_mem_size()

        # world-frame path has its own (smaller) scratch — no gravity-shim shared, no
        # main-sweep extras. Sized via gen_idsva_so_world_frame_temp_mem_size.
        idsva_so_world_frame_inner_temp_count = self.gen_idsva_so_world_frame_temp_mem_size() if self.robot.floating_base else idsva_so_body_frame_inner_temp_count
        idsva_so_world_frame_base_t_count = (2*nv + n) + idsva_so_world_frame_inner_temp_count + XI_size
        idsva_so_world_frame_full_t_count = idsva_so_world_frame_base_t_count + 4*nv**3
        self.idsva_so_world_frame_use_global_output = py_arena_bytes(idsva_so_world_frame_full_t_count) > self.cuda_target_shared_mem_bytes
        idsva_so_world_frame_t_count = idsva_so_world_frame_base_t_count if self.idsva_so_world_frame_use_global_output else idsva_so_world_frame_full_t_count

        # ----- FDSVA_SO shared-mem tier selection -----
        # Four nested tiers, ordered from least-spill to most-spill. Pick the
        # lowest-spill tier whose shared-arena bytes fit cuda_target_shared_mem.
        # Each tier sets three orthogonal state flags read by gen_fdsva_so_*:
        #   - use_global_tensors:  s_idsva_so + s_df2 (8*nv³ outputs) -> d_workspace
        #   - use_workspace_temp:  s_fdsva_temp (4*nv³ inner) -> d_workspace
        #   - fd_grad_use_spill:   fd_grad_inline's da_dq..fxvi band -> d_workspace grad section
        fdsva_so_base_t_count = 4*nv + nv*nv + nv + 2*nv*nv + XI_size
        fdsva_so_inner_temp_count = 4*nv**3
        fdsva_so_fd_gradient_inline_temp_count = self.gen_fdsva_so_fd_gradient_inline_temp_mem_size()
        fdsva_so_fd_gradient_inline_spilled_count = self.gen_fdsva_so_fd_gradient_inline_temp_mem_size_spilled()
        # fdsva_so dispatches to world_frame_inner for floating-base (smaller
        # footprint + no grav-shim spill) and body_frame_inner for fixed-base.
        fdsva_so_inner_idsva_so_temp_count = (
            idsva_so_world_frame_inner_temp_count if self.robot.floating_base
            else idsva_so_body_frame_inner_temp_count
        )
        _temp_full     = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_inner_temp_count, fdsva_so_fd_gradient_inline_temp_count)
        _temp_no_inner = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_fd_gradient_inline_temp_count)
        _temp_spilled  = max(fdsva_so_inner_idsva_so_temp_count, fdsva_so_fd_gradient_inline_spilled_count)
        # (name, shared_count, use_global_tensors, use_workspace_temp, fd_grad_use_spill)
        _fdsva_so_tiers = [
            ("full",                 fdsva_so_base_t_count + 8*nv**3 + _temp_full,  False, False, False),
            ("global_tensors",       fdsva_so_base_t_count + _temp_full,            True,  False, False),
            ("workspace_temp",       fdsva_so_base_t_count + _temp_no_inner,        True,  True,  False),
            ("workspace_temp_spill", fdsva_so_base_t_count + _temp_spilled,         True,  True,  True),
        ]
        _chosen = next(
            (t for t in _fdsva_so_tiers if py_arena_bytes(t[1]) <= self.cuda_target_shared_mem_bytes),
            _fdsva_so_tiers[-1],  # fallback: most-spill tier even if it still exceeds (runtime SKIP handles that)
        )
        _, fdsva_so_t_count, self.fdsva_so_use_global_tensors, self.fdsva_so_use_workspace_temp, self.fdsva_so_fd_grad_use_spill = _chosen
        grad_spill_workspace_t_count = max(id_du_temp_layout["spill_count"],
                                           id_du_temp_count,
                                           fd_du_temp_count,
                                           2*nv*nv)
        d2ee_workspace_t_count = 0
        if self.d2ee_use_workspace_temp:
            d2ee_workspace_t_count += d2ee_workspace_temp_count
        if self.d2ee_use_workspace_d2xhom:
            d2ee_workspace_t_count += d2Xhom_size
        # Include the floating-base gravity-shim spill (Phase D): the d2X/d2a/d2f
        # tensors live in d_workspace instead of shared memory for larger robots.
        idsva_so_body_frame_grav_spill_t_count = self.gen_floating_gravity_d2tau_dq_spill_count() if self.robot.floating_base else 0
        so_workspace_t_count = max(8*max(nv**3, 1), d2ee_workspace_t_count, idsva_so_body_frame_grav_spill_t_count)
        # Deprecated launch-count constants remain for external callers that still
        # pass COUNT*sizeof(T).  Make them conservative aliases for the byte arena
        # layouts so those callers do not under-allocate int topology helpers or
        # 16-byte alignment padding.
        legacy_count_pad = topology_count + 8
        legacy_arena_count = lambda t_count: int(t_count + legacy_count_pad)
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
                                 "const int GRID_FD_DU_USES_GLOBAL_TEMP = " + str(int(self.fd_du_use_global_temp)) + ";", \
                                 "const int GRID_ID_DU_USES_DA_DF_SPILL = " + str(int(self.id_du_use_selective_spill)) + ";", \
                                 "const int GRID_FD_DU_USES_DA_DF_SPILL = " + str(int(self.fd_du_use_selective_spill)) + ";", \
                                 "const int GRID_GENERATES_IDSVA_SO_BODY_FRAME = " + str(int(getattr(self, "generate_idsva_so_body_frame", True))) + ";", \
                                 "const int GRID_GENERATES_FDSVA_SO = " + str(int(getattr(self, "generate_fdsva_so", True))) + ";", \
                                 "const int GRID_GENERATES_D2EE = " + str(int(getattr(self, "generate_ee_pose_hessian", True))) + ";", \
                                 "const int GRID_IDSVA_SO_USES_GLOBAL_OUTPUT = " + str(int(self.idsva_so_body_frame_use_global_output)) + ";", \
                                 "const int GRID_FDSVA_SO_USES_GLOBAL_TENSORS = " + str(int(self.fdsva_so_use_global_tensors)) + ";", \
                                 "const int GRID_FDSVA_SO_USES_WORKSPACE_TEMP = " + str(int(self.fdsva_so_use_workspace_temp)) + ";", \
                                 "const int GRID_D2EE_USES_WORKSPACE_TEMP = " + str(int(self.d2ee_use_workspace_temp)) + ";", \
                                 "const int GRID_D2EE_USES_WORKSPACE_D2XHOM = " + str(int(self.d2ee_use_workspace_d2xhom)) + ";", \
                                 "const int GRID_D2EE_SHARED_TIER_VALUE = " + str(self.d2ee_spill_tier) + ";", \
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
                                 "const int SUGGESTED_THREADS = " + str(self.suggested_threads) + ";", \
                                 "",
                                 "// Resource-tier API (v2.0): each emitted kernel/_device/_inner takes a",
                                 "// `RESOURCE_TIER` template parameter that picks the (launch_bounds, smem,",
                                 "// register-footprint) profile. TIER_PERF is the default and is the",
                                 "// current-best perf; TIER_LITE keeps the same launch_bounds but reduces",
                                 "// smem footprint (some intermediates moved to workspace global mem);",
                                 "// TIER_MINIMAL drops launch_bounds to 1024 for maximum block-size flexibility",
                                 "// at the cost of register slack. Inline-CUDA power users with tight outer",
                                 "// kernels pick LITE/MINIMAL to fit GRiD primitives in their resource budget.",
                                 "constexpr int TIER_PERF    = 0;",
                                 "constexpr int TIER_LITE    = 1;",
                                 "constexpr int TIER_MINIMAL = 2;",
                                 "",
                                 "// Per-tier launch_bounds upper-bound (= max threads per block nvcc must",
                                 "// budget registers for). sm_120 has 65536 regs/block; nvcc enforces",
                                 "// regs_per_thread * max_threads <= regs_per_block, so a larger max_threads",
                                 "// directly caps regs_per_thread. PERF=SUGGESTED keeps current best perf;",
                                 "// LITE=min(2*SUGGESTED, 768) gives ~85 regs/thread cap (mid-budget);",
                                 "// MINIMAL=1024 gives ~64 regs/thread cap (maximum block-size flexibility).",
                                 "template <int TIER> __host__ __device__ constexpr int tier_max_threads() {",
                                 "    return (TIER == TIER_MINIMAL) ? 1024",
                                 "         : (TIER == TIER_LITE)    ? ((SUGGESTED_THREADS * 2 < 768) ? SUGGESTED_THREADS * 2 : 768)",
                                 "         :                          SUGGESTED_THREADS;",
                                 "}"])
        self.gen_add_code_lines([
                                 "#define GRID_GENERATED_NUM_JOINTS " + str(n),
                                 "#define GRID_GENERATED_NUM_EES " + str(self.robot.get_total_leaf_nodes()),
                                 ""])
        self.gen_add_code_lines([
                                 "template <typename T> __host__ __device__ inline size_t ID_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t MINV_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(minv_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t FD_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fd_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t ID_DU_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_du_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t FD_DU_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fd_du_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(integrator_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(max(integrator_du_t_count, integrator_du_with_x_kp1_t_count)) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t ID_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t MINV_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(minv_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t FD_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fd_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t ID_DU_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(id_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t FD_DU_DEVICE_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fd_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t ABA_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(aba_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t CRBA_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(crba_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ constexpr size_t GRID_EE_LINALG_SHARED_BYTES() { return static_cast<size_t>(0); }",
                                 "template <typename T> __host__ __device__ inline size_t EE_POS_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(ee_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t DEE_POS_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(dee_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t D2EE_POS_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(d2ee_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>()); }",
                                 "template <typename T> __host__ __device__ inline size_t IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(idsva_so_body_frame_t_count) + ", TOPOLOGY_HELPERS_COUNT); }",
                                 "template <typename T> __host__ __device__ inline size_t IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(idsva_so_world_frame_t_count) + ", TOPOLOGY_HELPERS_COUNT); }",
                                 "template <typename T> __host__ __device__ inline size_t FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES() { return grid_shared_arena_bytes<T>(" + str(fdsva_so_t_count) + ", TOPOLOGY_HELPERS_COUNT); }",
                                 "// Per-tier scratch sizes for fdsva_so_inner (inline-CUDA users only — the host launchers always use TIER_PERF).",
                                 "// At TIER_PERF the 4*NV^3 inner scratch lives in s_temp; at TIER_LITE/MINIMAL it moves to s_workspace, freeing shared memory for the caller's outer kernel.",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t FDSVA_SO_INNER_SMEM_BYTES() { return (TIER == TIER_PERF) ? sizeof(T) * static_cast<size_t>(" + str(4*nv**3) + ") : static_cast<size_t>(0); }",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t FDSVA_SO_INNER_WORKSPACE_BYTES() { return (TIER == TIER_PERF) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(4*nv**3) + "); }",
                                 "// Per-tier sizes for forward_dynamics_gradient_device (inline-CUDA users only). At TIER_PERF the temp scratch arena lives in s_temp; at TIER_LITE/MINIMAL it moves to s_workspace, freeing roughly " + str(fd_du_temp_count) + "*sizeof(T) bytes of smem.",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t FD_DU_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_PERF)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(fd_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(fd_du_device_t_count - fd_du_temp_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t FD_DU_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_PERF) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(fd_du_temp_count) + "); }",
                                 "// Per-tier sizes for end_effector_pose_gradient_hessian_device (inline-CUDA users only). At TIER_PERF d2eeTemp lives in the shared arena; at TIER_LITE/MINIMAL it moves to s_workspace, freeing " + str(d2ee_workspace_temp_count) + "*sizeof(T) bytes of smem.",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t D2EE_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_PERF)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(d2ee_inner_temp_count_shared + d2ee_workspace_temp_count + XHom_size + dXhom_size + d2Xhom_size) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(d2ee_inner_temp_count_shared + XHom_size + dXhom_size + d2Xhom_size) + ", TOPOLOGY_HELPERS_COUNT, GRID_EE_LINALG_SHARED_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t D2EE_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_PERF) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(d2ee_workspace_temp_count) + "); }",
                                 "// Per-tier sizes for inverse_dynamics_gradient_device (inline-CUDA users only). At TIER_PERF temp lives in s_temp; at TIER_LITE/MINIMAL it moves to s_workspace, freeing " + str(id_du_temp_count) + "*sizeof(T) bytes of smem.",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t ID_DU_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_PERF)",
                                 "        ? grid_shared_arena_bytes<T>(" + str(id_du_device_t_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>())",
                                 "        : grid_shared_arena_bytes<T>(" + str(id_du_device_t_count - id_du_temp_count) + ", TOPOLOGY_HELPERS_COUNT, GRID_LINALG_NVIDIA_MAX_HELPER_BYTES<T>());",
                                 "}",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t ID_DU_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_PERF) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(id_du_temp_count) + "); }",
                                 "// Per-tier sizes for idsva_so_device (inline-CUDA users only). At TIER_PERF temp lives in s_temp; at TIER_LITE/MINIMAL it moves to s_workspace, freeing " + str(idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + "*sizeof(T) bytes of smem. Frame picked at codegen time: " + ("world_frame" if self.robot.floating_base else "body_frame") + ".",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t IDSVA_SO_DEVICE_INLINE_SMEM_BYTES() {",
                                 "    return (TIER == TIER_PERF)",
                                 "        ? grid_shared_arena_bytes<T>(" + str((idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + XI_size) + ", TOPOLOGY_HELPERS_COUNT)",
                                 "        : grid_shared_arena_bytes<T>(" + str(XI_size) + ", TOPOLOGY_HELPERS_COUNT);",
                                 "}",
                                 "template <typename T, int TIER = TIER_PERF> __host__ __device__ constexpr size_t IDSVA_SO_DEVICE_INLINE_WORKSPACE_BYTES() { return (TIER == TIER_PERF) ? static_cast<size_t>(0) : sizeof(T) * static_cast<size_t>(" + str(idsva_so_world_frame_inner_temp_count if self.robot.floating_base else idsva_so_body_frame_inner_temp_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP() { return sizeof(T) * static_cast<size_t>(" + str(grad_spill_workspace_t_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_SO_WORKSPACE_BYTES_PER_TIMESTEP() { return sizeof(T) * static_cast<size_t>(" + str(so_workspace_t_count) + "); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_WORKSPACE_BYTES_PER_TIMESTEP() { return GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP<T>() + GRID_SO_WORKSPACE_BYTES_PER_TIMESTEP<T>(); }",
                                 "template <typename T> __host__ __device__ inline gridSharedTier GRID_ID_DU_SHARED_TIER() { return static_cast<gridSharedTier>(GRID_ID_DU_SHARED_TIER_VALUE); }",
                                 "template <typename T> __host__ __device__ inline gridSharedTier GRID_FD_DU_SHARED_TIER() { return static_cast<gridSharedTier>(GRID_FD_DU_SHARED_TIER_VALUE); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES() { return GRID_GRAD_WORKSPACE_BYTES_PER_TIMESTEP<T>(); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_TEMP_OFFSET_BYTES() { return GRID_SO_WORKSPACE_TEMP_OFFSET_BYTES<T>(); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_D2XHOM_OFFSET_BYTES() { return GRID_D2EE_WORKSPACE_TEMP_OFFSET_BYTES<T>(); }",
                                 "template <typename T> __host__ __device__ inline size_t GRID_D2EE_WORKSPACE_D2EETEMP_OFFSET_BYTES() { return GRID_D2EE_WORKSPACE_TEMP_OFFSET_BYTES<T>() + (GRID_D2EE_USES_WORKSPACE_D2XHOM ? sizeof(T) * static_cast<size_t>(D2XHOM_T_COUNT) : 0); }",
                                 "template <typename T> __host__ __device__ inline bool grid_selected_shared_memory_fits() { return ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES && FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES && (!GRID_GENERATES_D2EE || D2EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES) && (!GRID_GENERATES_IDSVA_SO_BODY_FRAME || IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES) && (!GRID_GENERATES_FDSVA_SO || FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>() <= GRID_CUDA_TARGET_SHARED_MEM_BYTES); }",
                                 "// __forceinline__ used throughout the xhom helper chain so ptxas folds these into the",
                                 "// inner kernels at all opt levels. For fixed-base the body of grid_q_index_affects_joint is",
                                 "// the trivial `q_index == joint_id` check that pre-GLASS callsites used directly.",
                                 "__host__ __device__ __forceinline__ bool grid_q_index_affects_joint(const int q_index, const int joint_id) {",
                                 ("    if (joint_id == 0) { return q_index >= 0 && q_index < 7; } return q_index == joint_id + 6;" if self.robot.floating_base else "    return q_index == joint_id;"),
                                 "}",
                                 "__host__ __device__ __forceinline__ int grid_d2xhom_offset(const int q_index_i, const int q_index_j) {",
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
                                 "__device__ inline T *grid_id_du_temp_ptr(T *s_temp, T *s_temp_spill, int index) {",
                                 "    if (!USE_DA_DF_SPILL) { return &s_temp[index]; }",
                                 "    if (index >= ID_DU_TEMP_SPILL_START && index < ID_DU_TEMP_SPILL_END) {",
                                 "        return &s_temp_spill[index - ID_DU_TEMP_SPILL_START];",
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
                                 "    // CPU INPUTS", \
                                 "    T *h_q_qd_u;", \
                                 "    T *h_q_qd;", \
                                 "    T *h_q;", \
                                 "    // GPU OUTPUTS", \
                                 "    T *d_c;", \
                                 "    T *d_Minv;", \
                                 "    T *d_qdd;", \
                                 "    T *d_M;", \
                                 "    T *d_dc_du;", \
                                 "    T *d_df_du;", \
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
                                 "    // CPU OUTPUTS", \
                                 "    T *h_c;", \
                                 "    T *h_Minv;", \
                                 "    T *h_qdd;", \
                                 "    T *h_M;", \
                                 "    T *h_dc_du;", \
                                 "    T *h_df_du;", \
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
                                 "};"])

    def gen_init_gridData(self):
        code_lines = ["gridData<T, KIND> *hd_data = (gridData<T, KIND> *)calloc(1, sizeof(gridData<T, KIND>));",
                      "const bool needs_dynamics = KIND == GRID_DATA_ALL || KIND == GRID_DATA_DYNAMICS;",
                      "const bool needs_kinematics = KIND == GRID_DATA_ALL || KIND == GRID_DATA_KINEMATICS;",
                      "// input variables used by dynamics and/or kinematics",
                      "if (needs_dynamics || needs_kinematics) {", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_q_qd_u, 3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_q, NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    hd_data->h_q_qd_u = (T *)malloc(3*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_q = (T *)malloc(NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
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
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_idsva_so, SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_df2, SECOND_ORDER_TENSOR_SIZE*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_workspace, GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*NUM_TIMESTEPS));", \
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
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_deePos, 6*NUM_EES*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    gpuErrchk(cudaMalloc((void**)&hd_data->d_d2eePos, 6*NUM_EES*NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T)));", \
                      "    if (GRID_D2EE_USES_WORKSPACE_TEMP && hd_data->d_workspace == nullptr) {gpuErrchk(cudaMalloc((void**)&hd_data->d_workspace, GRID_WORKSPACE_BYTES_PER_TIMESTEP<T>()*GRID_WORKSPACE_SLOTS*NUM_TIMESTEPS));}", \
                      "    hd_data->h_eePos = (T *)malloc(6*NUM_EES*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_deePos = (T *)malloc(6*NUM_EES*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "    hd_data->h_d2eePos = (T *)malloc(6*NUM_EES*NUM_JOINTS*NUM_JOINTS*NUM_TIMESTEPS*sizeof(T));", \
                      "}", \
                      "return hd_data;"]
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
    KERNEL_ATTR_MANIFEST = [
        # (algo_label, algo_short, gate_attr, bytes_macro, [(kernel_name<T>, signature), ...])
        ("inverse_dynamics", "id", None, "ID_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("inverse_dynamics_kernel<T>",
             "void (*)(T *, const T *, const int, const T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("direct_minv", "minv", None, "MINV_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("direct_minv_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("direct_minv_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("forward_dynamics", "fd", None, "FD_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("forward_dynamics_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("aba", "aba", None, "ABA_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("aba_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("aba_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("crba", "crba", None, "CRBA_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("crba_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("crba_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("end_effector_pose", "ee_pose", None, "EE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("end_effector_pose_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("end_effector_pose_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("end_effector_pose_gradient", "ee_pose_gradient", None, "DEE_POS_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("end_effector_pose_gradient_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
            ("end_effector_pose_gradient_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const int)"),
        ]),
        ("inverse_dynamics_gradient", "id_du", "generate_id_du", "ID_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("inverse_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const robotModel<T> *, const T, const int)"),
            ("inverse_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("forward_dynamics_gradient", "fd_du", "generate_fd_du", "FD_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("forward_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const T *, const T *, const robotModel<T> *, const T, const int)"),
            ("forward_dynamics_gradient_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("idsva_so_body_frame", "idsva_so_body_frame", "generate_idsva_so_body_frame", "IDSVA_SO_BODY_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("idsva_so_body_frame_kernel<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("idsva_so_body_frame_kernel_single_timing<T>",
             "void (*)(T *, unsigned char *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        # world-frame single-pass alternative (opt-in via enable_idsva_so_world_frame).
        # No d_workspace param — gravity is handled in the main sweep, not via shim.
        # Uses its own shared-mem macro (~25 KB for g1 vs shim's ~162 KB).
        ("idsva_so_world_frame", "idsva_so_world_frame", "generate_idsva_so_world_frame",
         "IDSVA_SO_WORLD_FRAME_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("idsva_so_world_frame_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
            ("idsva_so_world_frame_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const int)"),
        ]),
        ("fdsva_so", "fdsva_so", "generate_fdsva_so", "FDSVA_SO_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("fdsva_so_kernel<T>",
             "void (*)(T *, const T *, const int, unsigned char *, T *, const robotModel<T> *, const T, const int)"),
            ("fdsva_so_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, unsigned char *, T *, const robotModel<T> *, const T, const int)"),
        ]),
        ("integrator", "integrator", None, "INTEGRATOR_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("integrator_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
            ("integrator_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
        ]),
        ("integrator_gradient", "integrator_gradient", None, "INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("integrator_gradient_kernel<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
            ("integrator_gradient_kernel_single_timing<T>",
             "void (*)(T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
        ]),
        ("integrator_gradient_with_x_kp1", "integrator_with_gradient", None, "INTEGRATOR_DU_DYNAMIC_SHARED_MEM_BYTES<T>()", [
            ("integrator_gradient_with_x_kp1_kernel<T>",
             "void (*)(T *, T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
            ("integrator_gradient_with_x_kp1_kernel_single_timing<T>",
             "void (*)(T *, T *, const T *, const int, const robotModel<T> *, const T, const T, const int)"),
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
        attr_lines = ["// enable opt-in dynamic shared memory for every algorithm kernel"]
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
            # Kernels whose shared-mem may exceed the compile target for some
            # robots (d2ee + fdsva_so on large floating-base) get wrapped in a
            # compile-time-resolvable size guard so init_grid doesn't fail
            # registration when the kernel literally can't fit on a device
            # even with cudaFuncSetAttribute. idsva_so / world_frame are not
            # guarded — their runtime `grid_check_dynamic_shared_memory_bytes`
            # picks up the actual per-device cap (which can exceed the codegen
            # target on some GPUs), and we want them registered so that the
            # check + attribute setup happens in lockstep.
            guarded = algo_label in {
                "end_effector_pose_gradient_hessian",
                "fdsva_so",
            }
            indent = "    " if guarded else ""
            if guarded:
                attr_lines.append(f"if ({bytes_macro} <= GRID_CUDA_TARGET_SHARED_MEM_BYTES) {{")
            attr_lines.append(f"{indent}gpuErrchk(grid_check_dynamic_shared_memory_bytes(\"{algo_label}\", {bytes_macro}));")
            for kernel_name, signature in kernels:
                alias = f"_grid_kern_alias_{alias_counter}"
                alias_counter += 1
                attr_lines.append(f"{indent}auto {alias} = static_cast<{signature}>(&{kernel_name});")
                attr_lines.append(f"{indent}gpuErrchk(cudaFuncSetAttribute({alias}, cudaFuncAttributeMaxDynamicSharedMemorySize, {bytes_macro}));")
            if guarded:
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
                                 "gpuErrchk(cudaFree(hd_data->d_c)); gpuErrchk(cudaFree(hd_data->d_Minv)); gpuErrchk(cudaFree(hd_data->d_qdd)); gpuErrchk(cudaFree(hd_data->d_M));", \
                                 "gpuErrchk(cudaFree(hd_data->d_dc_du)); gpuErrchk(cudaFree(hd_data->d_df_du));", \
                                 "gpuErrchk(cudaFree(hd_data->d_eePos)); gpuErrchk(cudaFree(hd_data->d_deePos)); gpuErrchk(cudaFree(hd_data->d_d2eePos));", \
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

    # finally generate all of the code
    def gen_all_code(self, use_thread_group = False, include_base_inertia = False, include_homogenous_transforms = False, fixed_target_name = "", output_path = None,
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
            "    __device__ direct_minv_inner<T>(T *s_Minv, const T *s_q, T *s_XImats, int *s_topology_helpers, T *s_temp)",\
            "    __device__ direct_minv_device<T>(T *s_Minv, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ direct_minv_Kernel<T>(T *d_Minv, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   direct_minv<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ forward_dynamics_inner<T>(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, T *s_XImats, int *s_topology_helpers, T *s_temp, const T gravity)",\
            "    __device__ forward_dynamics_device<T>(T *s_qdd, const T *s_q, const T *s_qd, const T *s_u, const robotModel<T> *d_robotModel, const T gravity)", \
            "    __global__ forward_dynamics_kernel<T>(T *d_qdd, const T *d_q_qd_u, const robotModel<T> *d_robotModel, const T gravity, const int NUM_TIMESTEPS)", \
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
            "    __device__ end_effector_pose_inner<T>(T *s_eePos, const T *s_q, const T *s_Xhom, int *s_topology_helpers, T *s_temp)", \
            "    __device__ end_effector_pose_device<T>(T *s_eePos, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ end_effector_pose_kernel<T>(T *d_eePos, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
            "    __host__   end_effector_pose<T,USE_COMPRESSED_MEM=false>(gridData<T> *hd_data, const robotModel<T> *d_robotModel, const int num_timesteps, const dim3 block_dimms, const dim3 thread_dimms, cudaStream_t *streams)", \
            "",\
            "    __device__ end_effector_pose_gradient_inner<T>(T *s_deePos, const T *s_q, const T *s_Xhom, const T *s_dXhom, int *s_topology_helpers, T *s_temp)", \
            "    __device__ end_effector_pose_gradient_device<T>(T *s_deePos, const T *s_q, const robotModel<T> *d_robotModel)", \
            "    __global__ end_effector_pose_gradient_kernel<T>(T *d_deePos, const T *d_q, const robotModel<T> *d_robotModel, const int NUM_TIMESTEPS)", \
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
            "    __device__ fdsva_so_inner(T *s_df2, T *s_idsva_so, T *s_Minv, T *s_df_du, T *s_q, T *s_qd, const T *s_qdd, const T *s_tau, T *s_XImats, T *s_temp, const T gravity)",\
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
        self.gen_add_includes(use_thread_group)
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
        self.gen_invert_matrix(use_thread_group)
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
        self.gen_load_update_XImats_helpers(use_thread_group)
        if include_homogenous_transforms and include_any_kinematics:
            self.gen_load_update_XmatsHom_helpers(use_thread_group,include_base_inertia)
            if "ee_pose_gradient" in algorithms or "ee_pose_hessian" in algorithms:
                self.gen_load_update_XmatsHom_helpers(use_thread_group,include_base_inertia,include_gradients = True)
            if "ee_pose_hessian" in algorithms:
                self.gen_load_update_XmatsHom_helpers(use_thread_group,include_base_inertia,include_gradients = True, include_hessians = True)
        # then generate kinematic algorithms
        if include_any_kinematics:
            self.gen_eepose_and_derivatives(use_thread_group, fixed_target_name = fixed_target_name,
                                            include_pose = "ee_pose" in algorithms,
                                            include_gradient = "ee_pose_gradient" in algorithms,
                                            include_hessian = "ee_pose_hessian" in algorithms)
        if self.robot.floating_base and not enable_floating_second_order:
            print('floating-base second order dynamics are still under development')
        # then generate the dynamics algorithms
        if "id" in algorithms:
            self.gen_inverse_dynamics(use_thread_group)
        if "minv" in algorithms:
            self.gen_direct_minv(use_thread_group)
        if "fd" in algorithms:
            self.gen_forward_dynamics(use_thread_group)
        if "id_du" in algorithms:
            self.gen_inverse_dynamics_gradient(use_thread_group)
        if "fd_du" in algorithms:
            self.gen_forward_dynamics_gradient(use_thread_group)
        if "aba" in algorithms:
            self.gen_aba(use_thread_group)
        if "crba" in algorithms:
            self.gen_crba(use_thread_group)
        if "integrator" in algorithms:
            self.gen_integrator(use_thread_group)
        if ("integrator_gradient" in algorithms) or ("integrator_with_gradient" in algorithms):
            self.gen_integrator_gradient(use_thread_group)
        if not self.robot.floating_base or enable_floating_second_order:
            if "idsva_so_body_frame" in algorithms:
                self.gen_idsva_so_body_frame(use_thread_group)
                # Optional: emit the world-frame single-pass alternative path alongside
                # the existing emission. Co-exists with `idsva_so_body_frame_kernel`/`idsva_so_body_frame_host`;
                # the new entry point is `idsva_so_world_frame_kernel`/`idsva_so_world_frame_host`.
                if enable_idsva_so_world_frame:
                    self.gen_idsva_so_world_frame(use_thread_group)
                # Emit the dispatching `idsva_so` host wrapper. For floating-base
                # robots, requires world_frame to be enabled (it forwards there).
                # For fixed-base, forwards to body_frame.
                if (not self.robot.floating_base) or enable_idsva_so_world_frame:
                    self.gen_idsva_so_dispatcher()
            if "fdsva_so" in algorithms:
                self.gen_fdsva_so(use_thread_group)
        self.gen_combination_functions(algorithms, fixed_target_name)
        # then finally the master init and close the namespace
        self.gen_init_close_grid()
        self.gen_add_end_control_flow()
        # then output to a file
        if output_path is None:
            output_path = self.file_namespace + ".cuh"
        file = open(output_path, "w")
        file.write(self.code_str)
        file.close()
