#ifndef COLVARATOMS_KERNEL_H
#define COLVARATOMS_KERNEL_H

#include "colvar_rotation_derivative.h"
#include "colvarmodule.h"
// #include "colvar_gpu_support.h"
#include "colvartypes.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {

int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

// For debug gradients
int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaStream_t stream);

int change_one_coordinate(
  cvm::real* atoms_pos_ag,
  size_t atom_id_in_group, int xyz,
  cvm::real step_size,
  unsigned int num_atoms,
  cudaStream_t stream);

int atoms_calc_cog_com(
  const cvm::real* atoms_pos_ag,
  const cvm::real* atoms_mass,
  unsigned int num_atoms,
  cvm::rvector* d_cog_tmp,
  cvm::rvector* d_cog_out,
  cvm::rvector* d_cog_origin,
  cvm::rvector* d_com_tmp,
  cvm::rvector* d_com_out,
  cvm::rvector* h_cog_out,
  cvm::rvector* h_com_out,
  cvm::real total_mass,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int atoms_calc_cog(
  const cvm::real* atoms_pos_ag,
  unsigned int num_atoms,
  cvm::rvector* d_cog_tmp,
  cvm::rvector* d_cog_out,
  cvm::rvector* d_cog_origin,
  cvm::rvector* h_cog_out,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int atoms_total_force_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_total_force_proxy,
  cvm::real* atoms_total_force_ag,
  bool rotate,
  const cvm::quaternion* q,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaStream_t stream);

int apply_main_colvar_force_to_proxy(
  const int* atoms_proxy_index,
  cvm::real* atoms_applied_force_proxy,
  const cvm::real* atoms_grad_ag,
  cvm::real* colvar_force,
  bool rotate,
  const cvm::quaternion* q,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int apply_fitting_colvar_force_to_proxy(
  const int* atoms_proxy_index,
  cvm::real* atoms_applied_force_proxy,
  const cvm::real* atoms_grad_ag,
  cvm::real* colvar_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int accumulate_cpu_force(
  const cvm::real* h_atoms_force,
  cvm::real* d_atoms_force,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_fit_gradients_impl_loop1(
  const cvm::real* pos_unrotated,
  cvm::real* main_grad,
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const cvm::quaternion* q,
  unsigned int num_atoms_main,
  unsigned int num_atoms_fitting,
  double3* atom_grad,
  double* sum_dxdq,
  cvm::rmatrix* dxdC,
  unsigned int* tbcount,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_fit_gradients_impl_loop2(
  cvm::real* fit_grad,
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* atom_grad,
  const cvm::rmatrix* dxdC,
  unsigned int group_for_fit_size,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int apply_translation(
  cvm::real* atoms_pos_ag,
  cvm::real translation_vector_factor,
  const cvm::rvector* translation_vector,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int rotate_with_quaternion(
  cvm::real* atoms_pos_ag,
  const cvm::quaternion* q,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int apply_force_with_inverse_rotation(
  const cvm::real* atoms_force,
  const cvm::quaternion* q,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int apply_force(
  const cvm::real* atoms_force,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_fit_forces_impl_loop1(
  const cvm::real* pos_unrotated,
  cvm::real* main_force,
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const cvm::quaternion* q,
  unsigned int num_atoms_main,
  unsigned int num_atoms_fitting,
  double3* atom_grad,
  double* sum_dxdq,
  cvm::rmatrix* dxdC,
  unsigned int* tbcount,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_fit_forces_impl_loop2(
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* atom_grad,
  const cvm::rmatrix* dxdC,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  unsigned int group_for_fit_size,
  unsigned int proxy_stride,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

}

#elif defined(COLVARS_SYCL)
// TODO: SYCL
#endif
#endif // COLVARATOMS_KERNEL_H
