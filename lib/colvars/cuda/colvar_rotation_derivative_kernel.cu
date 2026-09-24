#include "colvar_rotation_derivative_kernel.h"

#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

namespace colvars_gpu {

__global__ void prepare_derivative_kernel(
  const rotation_derivative_dldq require_dl_dq,
  const cvm::real* S_eigval,
  const cvm::real* S_eigvec,
  cvm::real* tmp_Q0Q0,
  cvm::real* tmp_Q0Q0_L) {
  const int idx = threadIdx.x + blockIdx.x * blockDim.x;
  __shared__ cvm::real Q[4][4];
  __shared__ cvm::real L[4];
  constexpr int max_eigenvalue_index = 0;
  if (require_dl_dq & rotation_derivative_dldq::use_dl) {
    if (threadIdx.x == 0) {
      Q[0][0] = S_eigvec[max_eigenvalue_index*4+0];
      Q[0][1] = S_eigvec[max_eigenvalue_index*4+1];
      Q[0][2] = S_eigvec[max_eigenvalue_index*4+2];
      Q[0][3] = S_eigvec[max_eigenvalue_index*4+3];
    }
  }
  if (require_dl_dq & rotation_derivative_dldq::use_dq) {
    if (threadIdx.x == 0) {
      Q[0][0] = S_eigvec[max_eigenvalue_index*4+0];
      Q[0][1] = S_eigvec[max_eigenvalue_index*4+1];
      Q[0][2] = S_eigvec[max_eigenvalue_index*4+2];
      Q[0][3] = S_eigvec[max_eigenvalue_index*4+3];
      Q[1][0] = S_eigvec[(max_eigenvalue_index+1)*4+0];
      Q[1][1] = S_eigvec[(max_eigenvalue_index+1)*4+1];
      Q[1][2] = S_eigvec[(max_eigenvalue_index+1)*4+2];
      Q[1][3] = S_eigvec[(max_eigenvalue_index+1)*4+3];
      Q[2][0] = S_eigvec[(max_eigenvalue_index+2)*4+0];
      Q[2][1] = S_eigvec[(max_eigenvalue_index+2)*4+1];
      Q[2][2] = S_eigvec[(max_eigenvalue_index+2)*4+2];
      Q[2][3] = S_eigvec[(max_eigenvalue_index+2)*4+3];
      Q[3][0] = S_eigvec[(max_eigenvalue_index+3)*4+0];
      Q[3][1] = S_eigvec[(max_eigenvalue_index+3)*4+1];
      Q[3][2] = S_eigvec[(max_eigenvalue_index+3)*4+2];
      Q[3][3] = S_eigvec[(max_eigenvalue_index+3)*4+3];
      L[0] = S_eigval[max_eigenvalue_index+0];
      L[1] = S_eigval[max_eigenvalue_index+1];
      L[2] = S_eigval[max_eigenvalue_index+2];
      L[3] = S_eigval[max_eigenvalue_index+3];
    }
  }
  __syncthreads();
  if (idx < 16) {
    if (require_dl_dq & rotation_derivative_dldq::use_dl) {
      const int i = idx / 4;
      const int j = idx % 4;
      tmp_Q0Q0[i*4+j] = Q[0][i] * Q[0][j];
    }
  }
  __syncthreads();
  if (idx < 64) {
    if (require_dl_dq & rotation_derivative_dldq::use_dq) {
      const int i = idx / 16;
      const int j = (idx % 16) / 4;
      const int k = (idx % 16) % 4;
      tmp_Q0Q0_L[i*16+j*4+k] = (Q[1][j] * Q[0][k]) / (L[0] - L[1]) * Q[1][i] +
                               (Q[2][j] * Q[0][k]) / (L[0] - L[2]) * Q[2][i] +
                               (Q[3][j] * Q[0][k]) / (L[0] - L[3]) * Q[3][i];
    }
  }
}

int prepare_derivative(
  rotation_derivative_dldq dldq,
  const cvm::real* S_eigval,
  const cvm::real* S_eigvec,
  cvm::real* tmp_Q0Q0,
  cvm::real* tmp_Q0Q0_L,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = (dldq & rotation_derivative_dldq::use_dq) ? 64 : 16;
  const int num_blocks = 1;
  void* args[] = {
    &dldq,
    &S_eigval,
    &S_eigvec,
    &tmp_Q0Q0,
    &tmp_Q0Q0_L};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)prepare_derivative_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

}

#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
