#ifndef COLVARCOMP_DISTANCE_KERNEL_H
#define COLVARCOMP_DISTANCE_KERNEL_H

#include "colvarmodule.h"
#include "colvar_gpu_support.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {

int calc_value_rmsd(
  const cvm::real* ref_pos,
  const cvm::real* ag_pos,
  cvm::real* d_permutation_msds,
  cvm::real* h_permutation_msds,
  unsigned int num_atoms,
  unsigned int num_permutations,
  unsigned int num_ref_pos,
  unsigned int* tbcounts,
  std::vector<cudaGraphNode_t>& nodes,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_gradients_rmsd(
  const cvm::real* h_rmsd,
  const size_t* h_best_perm_index,
  const cvm::real* ref_pos,
  const cvm::real* ag_pos,
  cvm::real* grad,
  unsigned int num_atoms,
  unsigned int num_ref_pos,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int calc_force_invgrads_rmsd(
  const bool rotate,
  const int* atoms_proxy_index,
  const cvm::real* atoms_total_force_proxy,
  const cvm::quaternion* q,
  const cvm::real* grad,
  cvm::real* d_ft,
  cvm::real* h_ft,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

class rotation_derivative_gpu;

int calc_Jacobian_derivative_rmsd(
  const bool center,
  const bool rotate,
  const size_t* h_best_perm_index,
  const cvm::real* ref_pos,
  const cvm::quaternion* q,
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const cvm::real* h_rmsd,
  cvm::real* d_jd,
  cvm::real* h_jd,
  unsigned int num_atoms,
  unsigned int num_ref_pos,
  unsigned int* tbcount,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);
}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#endif // COLVARCOMP_DISTANCE_KERNEL_H
