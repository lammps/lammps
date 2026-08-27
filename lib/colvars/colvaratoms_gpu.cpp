#include "colvaratoms_gpu.h"
#include "colvar_gpu_support.h"
#include "colvardeps.h"
#include "colvarproxy.h"
#include "colvarmodule.h"
#include "colvarcomp.h"
#include "cuda/colvaratoms_kernel.h"

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

colvaratoms_gpu::colvaratoms_gpu(colvarmodule *cvmodule_in)
  : cvmodule(cvmodule_in) {
  std::memset(&gpu_buffers, 0, sizeof(gpu_buffers));
  std::memset(&debug_graphs, 0, sizeof(debug_graphs));
  std::memset(&calc_fit_gradients_gpu_info, 0, sizeof(calc_fit_gradients_gpu_info));
  std::memset(&calc_fit_forces_gpu_info, 0, sizeof(calc_fit_forces_gpu_info));
  h_sum_applied_colvar_force = nullptr;
  rot_deriv_gpu = nullptr;
  use_group_force = false;
  use_apply_colvar_force = false;
}

colvaratoms_gpu::~colvaratoms_gpu() {
  destroy_gpu();
}

int colvaratoms_gpu::init_gpu() {
  int error_code = COLVARS_OK;
  colvarproxy *p = cvmodule->proxy;
  // error_code |= checkGPUError(cudaStreamCreate(&stream));
  error_code |= p->reallocate_device(&gpu_buffers.d_com, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_com_tmp, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_cog, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_cog_tmp, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_cog_orig, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_ref_pos_cog, 1);
  error_code |= p->reallocate_device(&gpu_buffers.d_com_cog_tbcount, 1);
  error_code |= p->clear_device_array(gpu_buffers.d_com_cog_tbcount, 1);
  error_code |= p->clear_device_array(gpu_buffers.d_com_tmp, 1);
  error_code |= p->clear_device_array(gpu_buffers.d_cog_tmp, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_com, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_cog, 1);
  error_code |= p->reallocate_host(&gpu_buffers.h_cog_orig, 1);
  rot_deriv_gpu = nullptr;
  // error_code |= checkGPUError(cudaStreamCreate(&stream_ag_force));
  // std::memset(&calc_fit_gradients_info, 0, sizeof(calc_fit_gradients_info));
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_atom_grad, 1);
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_sum_dxdq, 4);
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_dxdC, 1);
  error_code |= p->reallocate_device(&calc_fit_gradients_gpu_info.d_tbcount, 1);
  error_code |= p->clear_device_array(calc_fit_gradients_gpu_info.d_sum_dxdq, 4);
  error_code |= p->clear_device_array(calc_fit_gradients_gpu_info.d_tbcount, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_atom_grad, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_sum_dxdq, 4);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_dxdC, 1);
  error_code |= p->reallocate_device(&calc_fit_forces_gpu_info.d_tbcount, 1);
  error_code |= p->clear_device_array(calc_fit_forces_gpu_info.d_sum_dxdq, 4);
  error_code |= p->clear_device_array(calc_fit_forces_gpu_info.d_tbcount, 1);
  error_code |= p->reallocate_host(&h_sum_applied_colvar_force, 1);
  h_sum_applied_colvar_force[0] = 0;
  use_apply_colvar_force = false;
  use_group_force = false;
  if (debug_graphs.graph_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphDestroy(
      debug_graphs.graph_calc_required_properties));
    debug_graphs.graph_calc_required_properties = nullptr;
  }
  if (debug_graphs.graph_exec_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphExecDestroy(
      debug_graphs.graph_exec_calc_required_properties));
    debug_graphs.graph_exec_calc_required_properties = nullptr;
  }
  debug_graphs.initialized = false;
  return error_code;
}

int colvaratoms_gpu::destroy_gpu() {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_index);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_pos);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_charge);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_vel);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_mass);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_grad);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_total_force);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_weight);
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_applied_force);
  error_code |= p->deallocate_device(&gpu_buffers.d_fit_gradients);
  error_code |= p->deallocate_device(&gpu_buffers.d_ref_pos);
  if (rot_deriv_gpu) {
    rot_deriv_gpu->~rotation_derivative_gpu();
    error_code |= p->deallocate_host(&rot_deriv_gpu);
  }
  error_code |= p->deallocate_device(&gpu_buffers.d_atoms_pos_unrotated);
  error_code |= p->deallocate_device(&gpu_buffers.d_com);
  error_code |= p->deallocate_device(&gpu_buffers.d_cog);
  error_code |= p->deallocate_device(&gpu_buffers.d_cog_orig);
  error_code |= p->deallocate_device(&gpu_buffers.d_ref_pos_cog);
  error_code |= p->deallocate_device(&gpu_buffers.d_com_cog_tbcount);
  error_code |= p->deallocate_device(&gpu_buffers.d_com_tmp);
  error_code |= p->deallocate_device(&gpu_buffers.d_cog_tmp);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_atom_grad);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_sum_dxdq);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_dxdC);
  error_code |= p->deallocate_device(&calc_fit_gradients_gpu_info.d_tbcount);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_atom_grad);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_sum_dxdq);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_dxdC);
  error_code |= p->deallocate_device(&calc_fit_forces_gpu_info.d_tbcount);
  error_code |= p->deallocate_host(&h_sum_applied_colvar_force);
  error_code |= p->deallocate_host(&gpu_buffers.h_com);
  error_code |= p->deallocate_host(&gpu_buffers.h_cog);
  error_code |= p->deallocate_host(&gpu_buffers.h_cog_orig);
  if (debug_graphs.graph_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphDestroy(
      debug_graphs.graph_calc_required_properties));
    debug_graphs.graph_calc_required_properties = nullptr;
  }
  if (debug_graphs.graph_exec_calc_required_properties) {
    error_code |= checkGPUError(cudaGraphExecDestroy(
      debug_graphs.graph_exec_calc_required_properties));
    debug_graphs.graph_exec_calc_required_properties = nullptr;
  }
  debug_graphs.initialized = false;
  return error_code;
}

int colvaratoms_gpu::sync_to_gpu_buffers(const cvm::atom_group* cpu_atoms) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_index, cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_charge, cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_mass, cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_weight, cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_pos, 3 * cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_vel, 3 * cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_total_force, 3 * cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_grad, 3 * cpu_atoms->num_atoms);
  error_code |= p->reallocate_device(&this->gpu_buffers.d_atoms_applied_force, 3 * cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_index.data(), this->gpu_buffers.d_atoms_index, cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_charge.data(), this->gpu_buffers.d_atoms_charge, cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_mass.data(), this->gpu_buffers.d_atoms_mass, cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_pos.data(), this->gpu_buffers.d_atoms_pos, 3 * cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_vel.data(), this->gpu_buffers.d_atoms_vel, 3 * cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_grad.data(), this->gpu_buffers.d_atoms_grad, 3 * cpu_atoms->num_atoms);
  error_code |= p->copy_HtoD(cpu_atoms->atoms_total_force.data(), this->gpu_buffers.d_atoms_total_force, 3 * cpu_atoms->num_atoms);
#elif defined(COLVARS_SYCL)
  // TODO: SYCL
#endif
  return error_code;
}

int colvaratoms_gpu::clear_gpu_buffers(const cvm::atom_group* cpu_atoms) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  size_t num_atoms = cpu_atoms->size();
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_index, num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_mass, num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_charge, num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_weight, num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_pos, 3 * num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_grad, 3 * num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_total_force, 3 * num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_vel, 3 * num_atoms);
  error_code |= p->clear_device_array(gpu_buffers.d_atoms_applied_force, 3 * num_atoms);
#elif defined(COLVARS_SYCL)
// TODO: COLVARS_SYCL
#endif
  return error_code;
}

int colvaratoms_gpu::add_reset_atoms_data_nodes(
  const cvm::atom_group* cpu_atoms,
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (!cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
#define ADD_CLEAR_FIELD_NODE(fieldName, numAtoms) do {\
cudaGraphNode_t clear_ ## fieldName ## _node ;\
error_code |= colvars_gpu::add_clear_array_node( \
  gpu_buffers.d_ ## fieldName, 3 * numAtoms, \
  clear_ ## fieldName ## _node , graph, {});\
  nodes_map[COLVARS_STRINGIFY(clear_ ## fieldName)] = clear_ ## fieldName ## _node;\
} while (0);
    // The clearing of atom positions is skipped since they are overwritten by the reading
    // ADD_CLEAR_FIELD_NODE(atoms_pos, cpu_atoms->num_atoms);
    // ADD_CLEAR_FIELD_NODE(atoms_vel, cpu_atoms->num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_grad, cpu_atoms->num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_total_force, cpu_atoms->num_atoms);
    ADD_CLEAR_FIELD_NODE(atoms_applied_force, cpu_atoms->num_atoms);
#undef ADD_CLEAR_FIELD_NODE
  }
  return error_code;
}

int colvaratoms_gpu::add_read_positions_nodes(
  const cvm::atom_group* cpu_atoms,
  cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (cpu_atoms->b_dummy) return error_code;
  colvarproxy *p = cvmodule->proxy;
  if (!cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
    cudaGraphNode_t read_positions_main;
    std::vector<cudaGraphNode_t> dependencies;
    error_code |= colvars_gpu::atoms_pos_from_proxy(
      gpu_buffers.d_atoms_index, p->proxy_atoms_positions_gpu(),
      gpu_buffers.d_atoms_pos, cpu_atoms->num_atoms, p->get_atom_ids()->size(),
      read_positions_main,
      graph, dependencies);
    nodes_map["read_positions_main"] = read_positions_main;
  }
  if (cpu_atoms->fitting_group) {
    cudaGraphNode_t read_positions_fitting;
    error_code |= colvars_gpu::atoms_pos_from_proxy(
      cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_index,
      p->proxy_atoms_positions_gpu(),
      cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
      cpu_atoms->fitting_group->num_atoms,
      p->get_atom_ids()->size(),
      read_positions_fitting,
      graph, {});
    nodes_map["read_positions_fitting"] =
      read_positions_fitting;
  }
  return error_code;
}

int colvaratoms_gpu::add_calc_required_properties_nodes(
  const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  const std::vector<cudaGraphNode_t>& extra_initial_dependencies) {
  int error_code = COLVARS_OK;
  // COM and COG of the main group
  if (cpu_atoms->b_dummy) {
    return error_code;
  } else if (cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
    error_code |= cvmodule->error("BUG: GPU buffers are not implemented with scalable atom group.\n");
    return error_code;
  } else {
    // Add kernel node for COG and COM
    cudaGraphNode_t calc_com_cog_main;
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    error_code |= colvars_gpu::prepare_dependencies(
      {{"read_positions_main", true}},
      dependencies, nodes_map, "calc_cog_com_main");
    error_code |= colvars_gpu::atoms_calc_cog_com(
      gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_mass, cpu_atoms->num_atoms,
      gpu_buffers.d_cog_tmp, gpu_buffers.d_cog, gpu_buffers.d_cog_orig,
      gpu_buffers.d_com_tmp, gpu_buffers.d_com,
      gpu_buffers.h_cog, gpu_buffers.h_com, cpu_atoms->total_mass,
      gpu_buffers.d_com_cog_tbcount,
      calc_com_cog_main, graph, dependencies);
    nodes_map["calc_com_cog_main"] = calc_com_cog_main;
  }

  // Fitting group cog
  if ((!cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) &&
      (cpu_atoms->is_enabled(colvardeps::f_ag_center) ||
       cpu_atoms->is_enabled(colvardeps::f_ag_rotate))) {
    if (cpu_atoms->fitting_group) {
      if (cpu_atoms->fitting_group->b_dummy) {
        // fitting_group->cog = fitting_group->dummy_atom_pos;
      } else {
        // Add COG kernel node
        std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
        error_code |= colvars_gpu::prepare_dependencies(
          {{"read_positions_fitting", true}},
          dependencies, nodes_map, "calc_cog_fitting");
        cudaGraphNode_t calc_cog_fitting;
        error_code |= colvars_gpu::atoms_calc_cog(
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
          cpu_atoms->fitting_group->num_atoms,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog_tmp,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog_orig,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.h_cog,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_com_cog_tbcount,
          calc_cog_fitting,
          graph, dependencies);
        nodes_map["calc_cog_fitting"] = calc_cog_fitting;
      }
    }

    // center on the origin first
    if (cpu_atoms->is_enabled(colvardeps::f_ag_center)) {
      const cvm::rvector* d_rpg_cog =
        cpu_atoms->fitting_group ?
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog :
        this->gpu_buffers.d_cog;
      // Apply translation
      std::vector<cudaGraphNode_t> dependencies;
      // d_rpg_cog could be the COG of the fitting group, so we need to wait for it
      error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_com_cog_main", false},
         {"read_positions_main", true},
         {"calc_cog_fitting", true}},
        dependencies, nodes_map, "apply_translation_main");
      cudaGraphNode_t move_to_origin_main;
      error_code |= colvars_gpu::apply_translation(
        gpu_buffers.d_atoms_pos, -1.0, d_rpg_cog,
        cpu_atoms->num_atoms, move_to_origin_main,
        graph, dependencies
      );
      nodes_map["move_to_origin_main"] = move_to_origin_main;
      if (cpu_atoms->fitting_group) {
        std::vector<cudaGraphNode_t> dependencies_fitting_group_translate;
        error_code |= colvars_gpu::prepare_dependencies(
          {{"calc_cog_fitting", false},
           {"read_positions_fitting", true}},
           dependencies_fitting_group_translate, nodes_map, "apply_translation_fitting");
        cudaGraphNode_t move_to_origin_fitting;
        error_code |= colvars_gpu::apply_translation(
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos, -1.0,
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog,
          cpu_atoms->fitting_group->num_atoms,
          move_to_origin_fitting,
          graph,
          dependencies_fitting_group_translate);
        nodes_map["move_to_origin_fitting"] = move_to_origin_fitting;
      }
    }

    // Save the unrotated frame for fit gradients
    if (cpu_atoms->is_enabled(colvardeps::f_ag_fit_gradients) && !cpu_atoms->b_dummy) {
      std::vector<cudaGraphNode_t> dependencies;
      error_code |= colvars_gpu::prepare_dependencies(
        {{"read_positions_main", true}, {"move_to_origin_main", true}},
        dependencies, nodes_map, "copy_unrotated_positions");
      cudaGraphNode_t copy_unrotated_positions;
      error_code |= colvars_gpu::add_copy_node(
        gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_pos_unrotated,
        3 * cpu_atoms->num_atoms, cudaMemcpyDeviceToDevice,
        copy_unrotated_positions,
        graph, dependencies);
      nodes_map["copy_unrotated_positions"] = copy_unrotated_positions;
    }

    // rotate the group (around the center of geometry if f_ag_center is
    // enabled, around the origin otherwise)
    if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
      auto* group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
      error_code |= rot_gpu.add_optimal_rotation_nodes(
        group_for_fit->gpu_atom_group->gpu_buffers.d_atoms_pos,
        gpu_buffers.d_ref_pos, group_for_fit->size(),
        cpu_atoms->num_ref_pos, graph, nodes_map);
      // Rotate atoms
      std::vector<cudaGraphNode_t> dependencies_rotate;
      // If we have f_ag_fit_gradients, then we need to wait for the copying of unrotated positions.
      // The atom group may be moved to origin (but not always).
      // The rotation should also wait after the COG and COM done
      error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_optimal_rotation", false},
         {"copy_unrotated_positions", true},
         {"move_to_origin_main", true},
         {"calc_com_cog_main", true}},
        dependencies_rotate, nodes_map, "rotate_main");
      cudaGraphNode_t rotate_main;
      error_code |= colvars_gpu::rotate_with_quaternion(
        gpu_buffers.d_atoms_pos, rot_gpu.get_q(), cpu_atoms->num_atoms,
        rotate_main, graph,
        dependencies_rotate);
      nodes_map["rotate_main"] = rotate_main;
      if (cpu_atoms->fitting_group) {
        std::vector<cudaGraphNode_t> dependencies_rotate_fitting_group;
        error_code |= colvars_gpu::prepare_dependencies(
          {{"calc_optimal_rotation", false},
           {"move_to_origin_fitting", true},
           {"calc_cog_fitting", true}},
          dependencies_rotate_fitting_group, nodes_map, "rotate_fitting");
        cudaGraphNode_t rotate_fitting;
        error_code |= colvars_gpu::rotate_with_quaternion(
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
          rot_gpu.get_q(),
          cpu_atoms->fitting_group->num_atoms,
          rotate_fitting, graph,
          dependencies_rotate_fitting_group);
        nodes_map["rotate_fitting"] = rotate_fitting;
      }
    }

    // align with the center of geometry of ref_pos
    if (cpu_atoms->is_enabled(colvardeps::f_ag_center) && !cpu_atoms->is_enabled(colvardeps::f_ag_center_origin)) {
      std::vector<cudaGraphNode_t> dependencies_move_to_ref_cog;
      // It looks like that the moving to COG of ref_pos can depend on
      // any of the following operations, so I add them all.
      error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_com_cog_main", false},
         {"rotate_main", true},
         {"copy_unrotated_positions", true},
         {"move_to_origin_main", true}},
         dependencies_move_to_ref_cog, nodes_map, "apply_translation_main_2");
      cudaGraphNode_t move_to_ref_cog_main;
      error_code |= colvars_gpu::apply_translation(
        gpu_buffers.d_atoms_pos, 1.0, gpu_buffers.d_ref_pos_cog, cpu_atoms->num_atoms,
        move_to_ref_cog_main, graph,
        dependencies_move_to_ref_cog);
      nodes_map["move_to_ref_cog_main"] = move_to_ref_cog_main;
      if (cpu_atoms->fitting_group) {
        std::vector<cudaGraphNode_t> dependencies_move_fitting_group_to_ref_cog;
        error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_cog_fitting", true},
         {"move_to_origin_fitting", true},
         {"rotate_fitting", true}},
         dependencies_move_fitting_group_to_ref_cog,
         nodes_map, "apply_translation_fitting_2");
        cudaGraphNode_t move_to_ref_cog_fitting;
        error_code |= colvars_gpu::apply_translation(
          cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos, 1.0,
          gpu_buffers.d_ref_pos_cog, cpu_atoms->fitting_group->num_atoms,
          move_to_ref_cog_fitting,
          graph,
          dependencies_move_fitting_group_to_ref_cog);
        nodes_map["move_to_ref_cog_fitting"] =
          move_to_ref_cog_fitting;
      }
    }

    // update COM and COG after fitting
    std::vector<cudaGraphNode_t> dependencies;
    error_code |= colvars_gpu::prepare_dependencies(
      {{"calc_com_cog_main", false},
      {"move_to_origin_main", true},
      {"rotate_main", true},
      {"move_to_ref_cog_main", true}},
      dependencies, nodes_map, "calc_cog_com_main_2");
    cudaGraphNode_t calc_com_cog_main_2;
    error_code |= colvars_gpu::atoms_calc_cog_com(
      gpu_buffers.d_atoms_pos, gpu_buffers.d_atoms_mass, cpu_atoms->num_atoms,
      gpu_buffers.d_cog_tmp, gpu_buffers.d_cog, nullptr,
      gpu_buffers.d_com_tmp, gpu_buffers.d_com,
      gpu_buffers.h_cog, gpu_buffers.h_com,
      cpu_atoms->total_mass, gpu_buffers.d_com_cog_tbcount,
      calc_com_cog_main_2, graph, dependencies);
    nodes_map["calc_com_cog_main_2"] = calc_com_cog_main_2;
    if (cpu_atoms->fitting_group) {
      dependencies.clear();
      // The following operations may or may not exist,
      // but if any of them exist, we need to wait for them
      error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_cog_fitting", false},
        {"move_to_origin_fitting", true},
        {"rotate_fitting", true},
        {"move_to_ref_cog_fitting", true}},
        dependencies, nodes_map, "calc_cog_fitting_2");
      cudaGraphNode_t calc_cog_fitting_2;
      error_code |= colvars_gpu::atoms_calc_cog(
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
        cpu_atoms->fitting_group->num_atoms,
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog_tmp,
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_cog,
        nullptr,
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.h_cog,
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_com_cog_tbcount,
        calc_cog_fitting_2, graph, dependencies);
      nodes_map["calc_cog_fitting_2"] = calc_cog_fitting_2;
    }
  }

  return error_code;
}

int colvaratoms_gpu::add_update_cpu_buffers_nodes(
  cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map) {
  int error_code = COLVARS_OK;
  if (cpu_atoms->b_dummy) return error_code;
  if (!cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
    std::vector<cudaGraphNode_t> dependencies;
    error_code |= colvars_gpu::prepare_dependencies(
      {{"read_positions_main", true},
       {"move_to_origin_main", true},
       {"rotate_main", true},
       {"move_to_ref_cog_main", true}},
      dependencies, nodes_map, "copy_atoms_to_host");
    cudaGraphNode_t copy_atoms_to_host_node;
    error_code |= colvars_gpu::add_copy_node(
      gpu_buffers.d_atoms_pos, cpu_atoms->atoms_pos.data(), 3 * cpu_atoms->num_atoms,
      cudaMemcpyDeviceToHost, copy_atoms_to_host_node,
      graph, dependencies);
    nodes_map["copy_atoms_to_host"] = copy_atoms_to_host_node;
    if (cpu_atoms->fitting_group) {
      dependencies.clear();
      error_code |= colvars_gpu::prepare_dependencies(
      {{"read_positions_fitting", true},
       {"move_to_origin_fitting", true},
       {"rotate_fitting", true},
       {"move_to_ref_cog_fitting", true}},
      dependencies, nodes_map, "copy_fitting_group_atoms_to_host");
      cudaGraphNode_t copy_fitting_group_atoms_to_host_node;
      error_code |= colvars_gpu::add_copy_node(
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
        cpu_atoms->fitting_group->atoms_pos.data(),
        3 * cpu_atoms->fitting_group->num_atoms,
        cudaMemcpyDeviceToHost,
        copy_fitting_group_atoms_to_host_node,
        graph, dependencies);
      nodes_map["copy_fitting_group_atoms_to_host"] = copy_fitting_group_atoms_to_host_node;
    }
    if (cpu_atoms->is_enabled(colvardeps::f_ag_center) || cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
      if (cpu_atoms->is_enabled(colvardeps::f_ag_fit_gradients) && !cpu_atoms->b_dummy) {
        if (cpu_atoms->atoms_pos_unrotated.size() != 3 * cpu_atoms->num_atoms) {
          cpu_atoms->atoms_pos_unrotated.resize(3 * cpu_atoms->num_atoms);
        }
        dependencies.clear();
        error_code |= colvars_gpu::prepare_dependencies(
          {{"copy_unrotated_positions", true}},
          dependencies, nodes_map, "copy_unrotated_positions_to_host");
        cudaGraphNode_t copy_unrotated_positions_to_host_node;
        error_code |= colvars_gpu::add_copy_node(
          gpu_buffers.d_atoms_pos_unrotated,
          cpu_atoms->atoms_pos_unrotated.data(), 3 * cpu_atoms->num_atoms,
          cudaMemcpyDeviceToHost,
          copy_unrotated_positions_to_host_node,
          graph, dependencies);
        nodes_map["copy_unrotated_positions_to_host"] = copy_unrotated_positions_to_host_node;
      }
    }
  }
  return error_code;
}

int colvaratoms_gpu::after_read_data_sync(
  cvm::atom_group* cpu_atoms, bool copy_to_cpu, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaStreamSynchronize(stream));
  // Update the COM
  if (cpu_atoms->b_dummy) {
    cpu_atoms->com = cpu_atoms->dummy_atom_pos;
    if (cvm::debug()) {
      cvmodule->log("Dummy atom center of mass = "+cvm::to_str(cpu_atoms->com)+"\n");
    }
  } else if (cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
    cpu_atoms->com = (cvmodule->proxy)->get_atom_group_com(cpu_atoms->index);
  } else {
    cpu_atoms->com.reset();
    if (copy_to_cpu) {
      cpu_atoms->com = *(gpu_buffers.h_com);
    }
  }
  // Update the COG
  if (cpu_atoms->b_dummy) {
    cpu_atoms->cog = cpu_atoms->dummy_atom_pos;
  } else {
    cpu_atoms->cog.reset();
    if (copy_to_cpu) {
      cpu_atoms->cog = *(gpu_buffers.h_cog);
    }
  }

  if (!cpu_atoms->is_enabled(colvardeps::f_ag_scalable)) {
    if (cpu_atoms->is_enabled(colvardeps::f_ag_center) || cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
      if (cpu_atoms->fitting_group) {
        // Update the fitting group COG
        if (cpu_atoms->fitting_group->b_dummy) {
          cpu_atoms->fitting_group->cog = cpu_atoms->fitting_group->dummy_atom_pos;
        } else {
          cpu_atoms->fitting_group->cog.reset();
          if (copy_to_cpu) {
            cpu_atoms->fitting_group->cog = *(cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.h_cog);
          }
        }
      }

      if (copy_to_cpu) {
        cpu_atoms->cog_orig = *(gpu_buffers.h_cog_orig);
      }
      if (cpu_atoms->fitting_group) {
        if (copy_to_cpu) {
          cpu_atoms->fitting_group->cog_orig = *(cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.h_cog_orig);
        }
      }

      if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
        if (copy_to_cpu) rot_gpu.to_cpu(cpu_atoms->rot);
        rot_gpu.after_sync_check();
      }
    }
  }
  return error_code;
}

int colvaratoms_gpu::add_calc_fit_gradients_nodes(
  cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  bool use_cpu_buffers) {
  int error_code = COLVARS_OK;
  if (cpu_atoms->b_dummy || !cpu_atoms->is_enabled(colvardeps::f_ag_fit_gradients)) return error_code;
  cvm::atom_group *group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
  // First, clear the temporary variables
  cudaGraphNode_t clear_atoms_grad_node;
  // cudaGraphNode_t clear_sum_dxdq_node;
  // cudaGraphNode_t clear_tbcount_node;
  error_code |= colvars_gpu::add_clear_array_node(
    calc_fit_gradients_gpu_info.d_atom_grad,
    1, clear_atoms_grad_node, graph, {});
  // error_code |= colvars_gpu::add_clear_array_node(
  //   calc_fit_gradients_gpu_info.d_sum_dxdq,
  //   1, clear_sum_dxdq_node, graph, {});
  nodes_map["clear_atoms_grad"] = clear_atoms_grad_node;
  // nodes_map["clear_sum_dxdq"] = clear_sum_dxdq_node;
  // If the CVC updates the gradients on CPU, then we need to copy them to GPU
  if (use_cpu_buffers) {
    cudaGraphNode_t copy_grad_HtoD_node;
    error_code |= colvars_gpu::add_copy_node(
      cpu_atoms->atoms_grad.data(), gpu_buffers.d_atoms_grad, 3 * cpu_atoms->num_atoms,
      cudaMemcpyHostToDevice, copy_grad_HtoD_node, graph,
      {clear_atoms_grad_node});
    nodes_map["copy_grad_HtoD"] = copy_grad_HtoD_node;
  }
  // Prepare the rotation derivative
  if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate) && rot_deriv_gpu) {
    error_code |= rot_deriv_gpu->add_prepare_derivative_nodes(
      rotation_derivative_dldq::use_dq, graph, nodes_map);
  }
  // Loop over the main group atoms (does not require rot_deriv_gpu)
  if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate) || cpu_atoms->is_enabled(colvardeps::f_ag_center)) {
    std::vector<cudaGraphNode_t> dependencies_main;
    error_code |= colvars_gpu::prepare_dependencies(
    {{"clear_atoms_grad", false},
     {"prepare_rotation_derivative", true},
     {"copy_grad_HtoD", true}}, dependencies_main,
     nodes_map, "calc_fit_gradients_loop1");
    cudaGraphNode_t calc_fit_forces_loop1_node;
    error_code |= colvars_gpu::calc_fit_gradients_impl_loop1(
      gpu_buffers.d_atoms_pos_unrotated, gpu_buffers.d_atoms_grad,
      rot_deriv_gpu, rot_gpu.get_q(),
      cpu_atoms->num_atoms, group_for_fit->size(),
      calc_fit_gradients_gpu_info.d_atom_grad,
      calc_fit_gradients_gpu_info.d_sum_dxdq,
      calc_fit_gradients_gpu_info.d_dxdC,
      calc_fit_gradients_gpu_info.d_tbcount,
      cpu_atoms->is_enabled(colvardeps::f_ag_center),
      cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
      calc_fit_forces_loop1_node,
      graph, dependencies_main);
    nodes_map["calc_fit_gradients_loop1"] = calc_fit_forces_loop1_node;
    // Loop over the fitting group
    cudaGraphNode_t calc_fit_forces_loop2_node;
    std::vector<cudaGraphNode_t> dependencies_fit_gradients;
    error_code |= colvars_gpu::prepare_dependencies(
      {{"calc_fit_gradients_loop1", false}},
      dependencies_fit_gradients, nodes_map, "calc_fit_gradients_loop2");
    error_code |= colvars_gpu::calc_fit_gradients_impl_loop2(
      group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients, rot_deriv_gpu,
      calc_fit_gradients_gpu_info.d_atom_grad,
      calc_fit_gradients_gpu_info.d_dxdC,
      group_for_fit->size(),
      cpu_atoms->is_enabled(colvardeps::f_ag_center),
      cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
      calc_fit_forces_loop2_node,
      graph, dependencies_fit_gradients);
    nodes_map["calc_fit_gradients_loop2"] = calc_fit_forces_loop2_node;
    if (use_cpu_buffers) {
      if (group_for_fit->fit_gradients.empty()) {
        group_for_fit->fit_gradients.resize(3 * group_for_fit->size());
      }
      cudaGraphNode_t copy_fit_gradients_to_host;
      std::vector<cudaGraphNode_t> dependencies_copy_fit_gradients;
      error_code |= colvars_gpu::prepare_dependencies(
        {{"calc_fit_gradients_loop2", false}},
        dependencies_copy_fit_gradients, nodes_map,
        "copy_fit_gradients_to_host");
      error_code |= colvars_gpu::add_copy_node(
        group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients,
        group_for_fit->fit_gradients.data(), 3 * group_for_fit->size(),
        cudaMemcpyDeviceToHost, copy_fit_gradients_to_host, graph,
        dependencies_copy_fit_gradients);
      nodes_map["copy_fit_gradients_to_host"] = copy_fit_gradients_to_host;
    }
  }
  return error_code;
}

int colvaratoms_gpu::begin_apply_force_gpu() {
  int error_code = COLVARS_OK;
  h_sum_applied_colvar_force[0] = 0;
  use_apply_colvar_force = false;
  use_group_force = false;
  return error_code;
}

int colvaratoms_gpu::add_apply_force_nodes(
  const cvm::atom_group* cpu_atoms, cudaGraph_t& graph,
  std::unordered_map<std::string, cudaGraphNode_t>& nodes_map,
  const std::vector<cudaGraphNode_t>& extra_initial_dependencies) {
  int error_code = COLVARS_OK;
  if (cpu_atoms->b_dummy) return error_code;
  colvarproxy* p = cvmodule->proxy;
  // Check if any of the parent CVCs require CPU buffers
  // Get all parents
  const std::vector<colvardeps*> parents = cpu_atoms->get_parents();
  std::vector<colvar::cvc*> parent_components;
  // Cast the parents to colvarcomp objects if possible
  for (size_t i = 0; i < parents.size(); ++i) {
    try {
      parent_components.push_back(dynamic_cast<colvar::cvc*>(parents[i]));
    } catch (const std::bad_cast& exception) {
      // Ignore the bad cast
    }
  }
  auto check_cvc_cpu_buffers = [](const colvar::cvc* obj){return obj->is_enabled(colvardeps::f_cvc_require_cpu_buffers);};
  const bool any_require_cpu_buffers = std::any_of(parent_components.begin(), parent_components.end(), check_cvc_cpu_buffers);
  const bool all_require_cpu_buffers = std::all_of(parent_components.begin(), parent_components.end(), check_cvc_cpu_buffers);
  if (use_apply_colvar_force) {
    const cvm::quaternion* q = rot_gpu.initialized() ? rot_gpu.get_q() : nullptr;
    if (any_require_cpu_buffers) {
      cudaGraphNode_t copy_grad_to_device;
      error_code |= colvars_gpu::add_copy_node(
        cpu_atoms->atoms_grad.data(), gpu_buffers.d_atoms_grad, 3 * cpu_atoms->num_atoms,
        cudaMemcpyHostToDevice, copy_grad_to_device, graph, {});
      nodes_map["copy_grad_to_device"] = copy_grad_to_device;
      if (!all_require_cpu_buffers) {
        const std::string error = "BUG: either none of the CVCs or "
          " all of the CVCs require CPU buffers!\n";
        return cvmodule->error(error);
      }
    }
    cudaGraphNode_t apply_main_colvar_force_to_proxy_node;
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    error_code |= colvars_gpu::prepare_dependencies(
      {{"copy_grad_to_device", true}},
      dependencies, nodes_map, "apply_main_colvar_force_to_proxy");
    error_code |= colvars_gpu::apply_main_colvar_force_to_proxy(
      gpu_buffers.d_atoms_index,
      p->proxy_atoms_new_colvar_forces_gpu(),
      gpu_buffers.d_atoms_grad,
      h_sum_applied_colvar_force,
      cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
      q, cpu_atoms->num_atoms, p->get_atom_ids()->size(),
      apply_main_colvar_force_to_proxy_node, graph, dependencies);
    nodes_map["apply_main_colvar_force_to_proxy"] = apply_main_colvar_force_to_proxy_node;
    if ((cpu_atoms->is_enabled(colvardeps::f_ag_center) ||
         cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) &&
        cpu_atoms->is_enabled(colvardeps::f_ag_fit_gradients)) {
      cudaGraphNode_t apply_fitting_colvar_force_to_proxy_node;
      auto* group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
      dependencies = extra_initial_dependencies;
      error_code |= colvars_gpu::apply_fitting_colvar_force_to_proxy(
        group_for_fit->gpu_atom_group->gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(),
        group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients,
        h_sum_applied_colvar_force, group_for_fit->size(),
        p->get_atom_ids()->size(),
        apply_fitting_colvar_force_to_proxy_node, graph, dependencies);
      nodes_map["apply_fitting_colvar_force_to_proxy"] = apply_fitting_colvar_force_to_proxy_node;
    }
  }
  if (use_group_force) {
    std::vector<cudaGraphNode_t> dependencies = extra_initial_dependencies;
    if (all_require_cpu_buffers) {
      // All CVCs using this group require CPU buffers, so
      // we can directly copy the data from group_forces to
      // d_atoms_applied_force.
      cudaGraphNode_t copy_forces_to_device;
      error_code |= colvars_gpu::add_copy_node(
        cpu_atoms->group_forces.data(), gpu_buffers.d_atoms_applied_force, 3 * cpu_atoms->num_atoms,
        cudaMemcpyHostToDevice, copy_forces_to_device, graph, dependencies);
      nodes_map["copy_forces_to_device"] = copy_forces_to_device;
    } else {
      if (any_require_cpu_buffers) {
        // Only some of the CVCs require the CPU buffers, so we
        // need to add the CPU buffers to d_atoms_applied_force
        cudaGraphNode_t accumulate_cpu_force_node;
        error_code |= colvars_gpu::accumulate_cpu_force(
          cpu_atoms->group_forces.data(), gpu_buffers.d_atoms_applied_force, cpu_atoms->num_atoms,
          accumulate_cpu_force_node, graph, dependencies);
        nodes_map["accumulate_cpu_force"] = accumulate_cpu_force_node;
      }
    }
    error_code |= colvars_gpu::prepare_dependencies(
      {{"copy_forces_to_device", true}, {"accumulate_cpu_force", true}},
      dependencies, nodes_map, "apply_force_main");
    if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
      // Rotate the forces back and add them to proxy
      cudaGraphNode_t apply_force_with_inverse_rotation_node;
      error_code |= colvars_gpu::apply_force_with_inverse_rotation(
        cpu_atoms->group_forces.data(), rot_gpu.get_q(), gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(), cpu_atoms->num_atoms,
        p->get_atom_ids()->size(),
        apply_force_with_inverse_rotation_node,
        graph, dependencies);
      nodes_map["apply_force_with_inverse_rotation"] =
        apply_force_with_inverse_rotation_node;
      // Prepare the derivative with respect to fit gradients
      if (rot_deriv_gpu) {
        error_code |= rot_deriv_gpu->add_prepare_derivative_nodes(
          rotation_derivative_dldq::use_dq, graph, nodes_map);
      }
    } else {
      // Just add the forces to proxy
      cudaGraphNode_t apply_force_main_node;
      error_code |= colvars_gpu::apply_force(
        cpu_atoms->group_forces.data(), gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(), cpu_atoms->num_atoms,
        p->get_atom_ids()->size(), apply_force_main_node,
        graph, dependencies);
      nodes_map["apply_force_main"] = apply_force_main_node;
    }
    if (cpu_atoms->is_enabled(colvardeps::f_ag_fit_gradients)) {
      auto* group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
      // Compute the forces on the fitting group and add them to proxy
      // Clear the temporary variables
      cudaGraphNode_t clear_atoms_grad_node;
      error_code |= colvars_gpu::add_clear_array_node(
        calc_fit_forces_gpu_info.d_atom_grad,
        1, clear_atoms_grad_node, graph, {});
      nodes_map["clear_atoms_grad"] = clear_atoms_grad_node;
      error_code |= colvars_gpu::prepare_dependencies(
        {{"clear_atoms_grad", false},
         {"prepare_rotation_derivative", true}},
        dependencies, nodes_map, "calc_fit_forces_loop1");
      cudaGraphNode_t calc_fit_forces_loop1_node;
      error_code |= colvars_gpu::calc_fit_forces_impl_loop1(
        gpu_buffers.d_atoms_pos_unrotated,
        gpu_buffers.d_atoms_applied_force,
        rot_deriv_gpu, rot_gpu.get_q(),
        cpu_atoms->num_atoms, group_for_fit->size(),
        calc_fit_forces_gpu_info.d_atom_grad,
        calc_fit_forces_gpu_info.d_sum_dxdq,
        calc_fit_forces_gpu_info.d_dxdC,
        calc_fit_forces_gpu_info.d_tbcount,
        cpu_atoms->is_enabled(colvardeps::f_ag_center),
        cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
        calc_fit_forces_loop1_node,
        graph, dependencies);
      nodes_map["calc_fit_forces_loop1"] = calc_fit_forces_loop1_node;
      // Compute the forces on the fitting group and add them back to proxy
      cudaGraphNode_t calc_fit_forces_loop2_node;
      error_code |= colvars_gpu::calc_fit_forces_impl_loop2(
        rot_deriv_gpu,
        calc_fit_forces_gpu_info.d_atom_grad,
        calc_fit_forces_gpu_info.d_dxdC,
        group_for_fit->gpu_atom_group->gpu_buffers.d_atoms_index,
        p->proxy_atoms_new_colvar_forces_gpu(),
        group_for_fit->size(), p->get_atom_ids()->size(),
        cpu_atoms->is_enabled(colvardeps::f_ag_center),
        cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
        calc_fit_forces_loop2_node, graph,
        {calc_fit_forces_loop1_node});
      nodes_map["calc_fit_forces_loop2"] = calc_fit_forces_loop2_node;
    }
  }
  return error_code;
}

int colvaratoms_gpu::read_positions_gpu_debug(
  cvm::atom_group* cpu_atoms, bool change_fitting_group,
  size_t change_atom_i, int xyz,
  bool to_cpu, double sign, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  colvarproxy *p = cvmodule->proxy;
  error_code |= colvars_gpu::atoms_pos_from_proxy(
    gpu_buffers.d_atoms_index, p->proxy_atoms_positions_gpu(),
    gpu_buffers.d_atoms_pos, cpu_atoms->num_atoms, p->get_atom_ids()->size(),
    stream);
  if (cpu_atoms->fitting_group) {
    error_code |= colvars_gpu::atoms_pos_from_proxy(
      cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_index,
      p->proxy_atoms_positions_gpu(),
      cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
      cpu_atoms->fitting_group->num_atoms, p->get_atom_ids()->size(),
      stream);
  }
  if (!change_fitting_group) {
    error_code |= colvars_gpu::change_one_coordinate(
      gpu_buffers.d_atoms_pos, change_atom_i, xyz,
      sign * cvmodule->debug_gradients_step_size, cpu_atoms->num_atoms, stream);
  } else {
    if (cpu_atoms->fitting_group) {
      error_code |= colvars_gpu::change_one_coordinate(
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos, change_atom_i, xyz,
        sign * cvmodule->debug_gradients_step_size, cpu_atoms->fitting_group->num_atoms, stream);
    }
  }
  if (to_cpu) {
    error_code |= p->copy_DtoH_async(
      gpu_buffers.d_atoms_pos, cpu_atoms->atoms_pos.data(), 3 * cpu_atoms->num_atoms, stream);
    if (cpu_atoms->fitting_group) {
      error_code |= p->copy_DtoH_async(
        cpu_atoms->fitting_group->gpu_atom_group->gpu_buffers.d_atoms_pos,
        cpu_atoms->fitting_group->atoms_pos.data(),
        3 * cpu_atoms->fitting_group->num_atoms, stream);
    }
  }
  error_code |= checkGPUError(cudaStreamSynchronize(stream));
  return error_code;
}

int colvaratoms_gpu::calc_required_properties_gpu_debug(
  cvm::atom_group* cpu_atoms, bool to_cpu, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (!debug_graphs.initialized) {
    // Create the debug graph
    error_code |= checkGPUError(cudaGraphCreate(&debug_graphs.graph_calc_required_properties, 0));
    std::unordered_map<std::string, cudaGraphNode_t> nodes_map;
    error_code |= add_calc_required_properties_nodes(
      cpu_atoms, debug_graphs.graph_calc_required_properties, nodes_map);
    if (to_cpu) {
      error_code |= add_update_cpu_buffers_nodes(
        cpu_atoms, debug_graphs.graph_calc_required_properties, nodes_map);
    }
    cudaGraphInstantiateParams params{0};
    params.flags = cudaGraphInstantiateFlagUpload;
    params.uploadStream = stream;
    error_code |= checkGPUError(cudaGraphInstantiateWithParams(
      &debug_graphs.graph_exec_calc_required_properties, debug_graphs.graph_calc_required_properties, &params));
    if (params.result_out != cudaGraphInstantiateSuccess) {
      error_code |= cvmodule->error("Failed to instantiate CUDA graph!", COLVARS_ERROR);
    }
    debug_graphs.initialized = true;
  }
  error_code |= checkGPUError(cudaGraphLaunch(
    debug_graphs.graph_exec_calc_required_properties, stream));
  return error_code;
}

void colvaratoms_gpu::do_feature_side_effects_gpu(
  cvm::atom_group* cpu_atoms, int id) {
  if (cvm::debug()) {
    cvmodule->log("cvm::atom_group::do_feature_side_effects_gpu.\n");
  }
  switch (id) {
    case colvardeps::f_ag_fit_gradients: {
      colvarproxy* p = cvmodule->proxy;
      if (gpu_buffers.d_atoms_pos_unrotated == nullptr) {
        p->allocate_device(&gpu_buffers.d_atoms_pos_unrotated, 3 * cpu_atoms->num_atoms);
      }
      if (cpu_atoms->is_enabled(colvardeps::f_ag_center) ||
          cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
        cvm::atom_group *group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
        if (group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients == nullptr) {
          p->allocate_device(&group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients,
                              3 * group_for_fit->size());
          p->clear_device_array(group_for_fit->gpu_atom_group->gpu_buffers.d_fit_gradients, 3 * group_for_fit->size());
        }
      }
      break;
    }
    case colvardeps::f_ag_rotate: {
      rot_gpu.init();
      break;
    }
  }
}

int colvaratoms_gpu::setup_rotation(const cvm::atom_group* cpu_atoms) {
  int error_code = COLVARS_OK;
  colvarproxy* p = cvmodule->proxy;
  error_code |= p->reallocate_device(&gpu_buffers.d_ref_pos, cpu_atoms->ref_pos.size());
  error_code |=p->copy_HtoD(cpu_atoms->ref_pos.data(), gpu_buffers.d_ref_pos, cpu_atoms->ref_pos.size());
  error_code |=p->copy_HtoD(&cpu_atoms->ref_pos_cog, gpu_buffers.d_ref_pos_cog, 1);
  rot_gpu.init();
  return error_code;
}

int colvaratoms_gpu::setup_rotation_derivative(const cvm::atom_group* cpu_atoms) {
  int error_code = COLVARS_OK;
  auto* group_for_fit = cpu_atoms->fitting_group ? cpu_atoms->fitting_group : cpu_atoms;
  colvarproxy* p = cvmodule->proxy;
  if (rot_deriv_gpu != nullptr) {
    rot_deriv_gpu->~rotation_derivative_gpu();
    error_code |= p->deallocate_host(&rot_deriv_gpu);
  }
  error_code |= p->allocate_host(&rot_deriv_gpu, 1);
  rot_deriv_gpu = new (rot_deriv_gpu) colvars_gpu::rotation_derivative_gpu(cvmodule);
  error_code |= rot_deriv_gpu->init(
    &rot_gpu,
    group_for_fit->gpu_atom_group->gpu_buffers.d_atoms_pos,
    gpu_buffers.d_ref_pos, group_for_fit->size(), cpu_atoms->num_ref_pos);
  return error_code;
}

int colvaratoms_gpu::read_total_forces(cvm::atom_group* cpu_atoms) {
  int error_code = COLVARS_OK;
  colvarproxy *p = cvmodule->proxy;
  cvm::quaternion* q_ptr = nullptr;
  if (cpu_atoms->is_enabled(colvardeps::f_ag_rotate)) {
    q_ptr = rot_gpu.get_q();
  }
  error_code |= colvars_gpu::atoms_total_force_from_proxy(
    gpu_buffers.d_atoms_index, p->proxy_atoms_total_forces_gpu(),
    gpu_buffers.d_atoms_total_force, cpu_atoms->is_enabled(colvardeps::f_ag_rotate),
    q_ptr, cpu_atoms->num_atoms, p->get_atom_ids()->size(), p->get_default_stream());
  // TODO: How can I check if the CVC has GPU implementation?
  // If the CVC only supports CPU, copy the data to host
  error_code |= p->copy_DtoH_async(
    gpu_buffers.d_atoms_total_force, cpu_atoms->atoms_total_force.data(),
    3 * cpu_atoms->num_atoms, p->get_default_stream());
  error_code |= checkGPUError(cudaStreamSynchronize(p->get_default_stream()));
  return error_code;
}

void colvaratoms_gpu::apply_colvar_force_from_cpu(cvm::real const& cpu_force) {
  use_apply_colvar_force = true;
  h_sum_applied_colvar_force[0] += cpu_force;
}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
}
