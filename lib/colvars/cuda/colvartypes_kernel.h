#ifndef COLVARTYPES_KERNEL_H
#define COLVARTYPES_KERNEL_H

#include "colvarmodule.h"
#include "colvar_gpu_support.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {

int build_overlapping_matrix(
  const cvm::real* pos1_x,
  const cvm::real* pos1_y,
  const cvm::real* pos1_z,
  const cvm::real* pos2_x,
  const cvm::real* pos2_y,
  const cvm::real* pos2_z,
  cvm::real* S,
  cvm::real* S_eigvec,
  cvm::rmatrix* h_C,
  unsigned int* tbcount,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

int jacobi_4x4(
  double* S_eigvec,
  double* S_eigval,
  int* max_reached,
  cvm::quaternion* q,
  bool monitor_crossings,
  cvm::real crossing_threshold,
  cvm::quaternion* q_old,
  int* discontinuous_rotation,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)

#endif // COLVARTYPES_KERNEL_H
