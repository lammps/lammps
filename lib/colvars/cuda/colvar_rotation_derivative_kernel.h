#ifndef COLVAR_ROTATION_DERIVATIVE_KERNEL_H
#define COLVAR_ROTATION_DERIVATIVE_KERNEL_H

#include "colvarmodule.h"
#include "colvar_gpu_support.h"
#include "colvar_rotation_derivative.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {

int prepare_derivative(
  rotation_derivative_dldq dldq,
  const cvm::real* S_eigval,
  const cvm::real* S_eigvec,
  cvm::real* tmp_Q0Q0,
  cvm::real* tmp_Q0Q0_L,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies);

}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
#endif // COLVAR_ROTATION_DERIVATIVE_KERNEL_H
