#include "colvar_gpu_support.h"
#include "colvartypes_kernel.h"
#include "colvartypes.h"

#if defined(COLVARS_CUDA)
#include <cub/block/block_reduce.cuh>
#endif

#if defined (COLVARS_HIP)
#include <hipcub/block/block_reduce.hpp>
#define cub hipcub
#endif

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)

template <int BLOCK_SIZE>
__global__ void build_overlapping_matrix_kernel(
  const cvm::real* __restrict pos1_x,
  const cvm::real* __restrict pos1_y,
  const cvm::real* __restrict pos1_z,
  const cvm::real* __restrict pos2_x,
  const cvm::real* __restrict pos2_y,
  const cvm::real* __restrict pos2_z,
  cvm::real* __restrict S,
  cvm::real* __restrict S_eigvec,
  cvm::rmatrix* __restrict h_C,
  unsigned int* __restrict tbcount,
  int num_atoms) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::rmatrix C;
  C.reset();
  while (i < num_atoms) {
    C.xx += pos1_x[i] * pos2_x[i];
    C.xy += pos1_x[i] * pos2_y[i];
    C.xz += pos1_x[i] * pos2_z[i];
    C.yx += pos1_y[i] * pos2_x[i];
    C.yy += pos1_y[i] * pos2_y[i];
    C.yz += pos1_y[i] * pos2_z[i];
    C.zx += pos1_z[i] * pos2_x[i];
    C.zy += pos1_z[i] * pos2_y[i];
    C.zz += pos1_z[i] * pos2_z[i];
    i += gridSize;
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  C.xx = BlockReduce(temp_storage).Sum(C.xx); __syncthreads();
  C.xy = BlockReduce(temp_storage).Sum(C.xy); __syncthreads();
  C.xz = BlockReduce(temp_storage).Sum(C.xz); __syncthreads();
  C.yx = BlockReduce(temp_storage).Sum(C.yx); __syncthreads();
  C.yy = BlockReduce(temp_storage).Sum(C.yy); __syncthreads();
  C.yz = BlockReduce(temp_storage).Sum(C.yz); __syncthreads();
  C.zx = BlockReduce(temp_storage).Sum(C.zx); __syncthreads();
  C.zy = BlockReduce(temp_storage).Sum(C.zy); __syncthreads();
  C.zz = BlockReduce(temp_storage).Sum(C.zz); __syncthreads();
  if (threadIdx.x == 0) {
    // S is 4 x 4 so I can use it as a temporary buffer
    atomicAdd(&(S[0]), C.xx);
    atomicAdd(&(S[1]), C.xy);
    atomicAdd(&(S[2]), C.xz);
    atomicAdd(&(S[3]), C.yx);
    atomicAdd(&(S[4]), C.yy);
    atomicAdd(&(S[5]), C.yz);
    atomicAdd(&(S[6]), C.zx);
    atomicAdd(&(S[7]), C.zy);
    atomicAdd(&(S[8]), C.zz);
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      C.xx = S[0];
      C.xy = S[1];
      C.xz = S[2];
      C.yx = S[3];
      C.yy = S[4];
      C.yz = S[5];
      C.zx = S[6];
      C.zy = S[7];
      C.zz = S[8];
      // Now we can use the first thread of the last block to set S
      S[0*4+0] =   C.xx + C.yy + C.zz;
      S[1*4+0] =   C.yz - C.zy;
      S[0*4+1] =   C.yz - C.zy;
      S[2*4+0] = - C.xz + C.zx ;
      S[0*4+2] = - C.xz + C.zx ;
      S[3*4+0] =   C.xy - C.yx;
      S[0*4+3] =   C.xy - C.yx;
      S[1*4+1] =   C.xx - C.yy - C.zz;
      S[2*4+1] =   C.xy + C.yx;
      S[1*4+2] =   C.xy + C.yx;
      S[3*4+1] =   C.xz + C.zx;
      S[1*4+3] =   C.xz + C.zx;
      S[2*4+2] = - C.xx + C.yy - C.zz;
      S[3*4+2] =   C.yz + C.zy;
      S[2*4+3] =   C.yz + C.zy;
      S[3*4+3] = - C.xx - C.yy + C.zz;
      S_eigvec[0*4+0] =   C.xx + C.yy + C.zz;
      S_eigvec[1*4+0] =   C.yz - C.zy;
      S_eigvec[0*4+1] =   C.yz - C.zy;
      S_eigvec[2*4+0] = - C.xz + C.zx ;
      S_eigvec[0*4+2] = - C.xz + C.zx ;
      S_eigvec[3*4+0] =   C.xy - C.yx;
      S_eigvec[0*4+3] =   C.xy - C.yx;
      S_eigvec[1*4+1] =   C.xx - C.yy - C.zz;
      S_eigvec[2*4+1] =   C.xy + C.yx;
      S_eigvec[1*4+2] =   C.xy + C.yx;
      S_eigvec[3*4+1] =   C.xz + C.zx;
      S_eigvec[1*4+3] =   C.xz + C.zx;
      S_eigvec[2*4+2] = - C.xx + C.yy - C.zz;
      S_eigvec[3*4+2] =   C.yz + C.zy;
      S_eigvec[2*4+3] =   C.yz + C.zy;
      S_eigvec[3*4+3] = - C.xx - C.yy + C.zz;
      // Save the data to host memory
      memcpy(h_C, &C, sizeof(cvm::rmatrix));
      // Reset the atomic counter
      tbcount[0] = 0;
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  // if (num_atoms == 0) return;
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  void* args[] = {
    &pos1_x, &pos1_y, &pos1_z,
    &pos2_x, &pos2_y, &pos2_z,
    &S, &S_eigvec, &h_C, &tbcount, &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           =
    (void*)build_overlapping_matrix_kernel<block_size>;
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

#define JACOBI_MAX_ITERATION 50
#define JACOBI_TOLERANCE 1e-8
template <int p, int q>
__inline__ __device__ void apply_jacobi(
  double* __restrict A,
  const double c,
  const double s,
  const double c2,
  const double s2,
  const double cs) {
  #pragma unroll
  for (int i = 0; i < 4; ++i) {
    const double oip = A[i*4+p];
    const double oiq = A[i*4+q];
    if (i != p && i != q) {
      A[i*4+p] = c * oip - s * oiq;
      A[p*4+i] = A[i*4+p];
      A[i*4+q] = c * oiq + s * oip;
      A[q*4+i] = A[i*4+q];
    }
  }
  const double opp = A[p*4+p];
  const double oqq = A[q*4+q];
  const double opq = A[p*4+q];
  A[p*4+p] = c2 * opp + s2 * oqq - 2.0 * cs * opq;
  A[q*4+q] = s2 * opp + c2 * oqq + 2.0 * cs * opq;
  A[p*4+q] = 0;
  A[q*4+p] = 0;
}

template <int p, int q>
__inline__ __device__ void multiply_jacobi(
  double* __restrict V, double c, double s) {
  #pragma unroll
  for (int i = 0; i < 4; ++i) {
    const double oip = V[i*4+p];
    const double oiq = V[i*4+q];
    V[i*4+p] = c * oip - s * oiq;
    V[i*4+q] = s * oip + c * oiq;
  }
}

__inline__ __device__ void compute_c_s(
  double a_pq, double a_pp, double a_qq, double& c, double& s, double& c2, double& s2, double& cs) {
  const double theta = 0.5 * (a_qq - a_pp) / a_pq;
  const double t = 1 / (sqrt(theta * theta + 1.0) + fabs(theta));
  c = rsqrt(t * t + 1.0);
  s = theta < 0 ? -t * c : t * c;
  c2 = c * c;
  s2 = s * s;
  cs = c * s;
}

// Use exactly 2 threads
__global__ void jacobi_4x4_kernel(
  double* A_in, double* eigvals, int* max_reached,
  cvm::quaternion* q,
  bool monitor_crossings,
  cvm::real crossing_threshold,
  cvm::quaternion* q_old,
  int* discontinuous_rotation) {
  __shared__ double A[4*4];
  __shared__ double V[4*4];
  const int idx = threadIdx.x;
  if (max_reached && idx == 0) {
    max_reached[0] = 0;
    __threadfence();
  }
  if (idx == 0) {
    memset(V, 0, sizeof(double)*4*4);
    V[0*4+0] = 1;
    V[1*4+1] = 1;
    V[2*4+2] = 1;
    V[3*4+3] = 1;
    A[0] = A_in[0];
    A[1] = A_in[1];
    A[2] = A_in[2];
    A[3] = A_in[3];
    A[4] = A_in[4];
    A[5] = A_in[5];
    A[6] = A_in[6];
    A[7] = A_in[7];
    A[8] = A_in[8];
    A[9] = A_in[9];
    A[10] = A_in[10];
    A[11] = A_in[11];
    A[12] = A_in[12];
    A[13] = A_in[13];
    A[14] = A_in[14];
    A[15] = A_in[15];
  }
  __syncthreads();
  constexpr int p_ids[] = {0, 2, 0, 1, 0, 1};
  constexpr int q_ids[] = {1, 3, 2, 3, 3, 2};
  double off_diag_sum =
    fabs(A[0*4+1]) + fabs(A[0*4+2]) + fabs(A[0*4+3]) +
    fabs(A[1*4+2]) + fabs(A[1*4+3]) +
    fabs(A[2*4+3]);
  int iteration = 0;
  while (off_diag_sum > JACOBI_TOLERANCE) {
    double c = 0, s = 0, c2 = 0, s2 = 0, cs = 0;
    bool rotate = false;
    int p = p_ids[idx];
    int q = q_ids[idx];
    double a_pq = A[p*4+q];
    if (fabs(a_pq) > 0) {
      rotate = true;
      const double a_pp = A[p*4+p];
      const double a_qq = A[q*4+q];
      compute_c_s(a_pq, a_pp, a_qq, c, s, c2, s2, cs);
    }
    COLVARS_SYNC_WARP;
    if (idx == 0 && rotate) {
      apply_jacobi<0, 1>(A, c, s, c2, s2, cs);
      multiply_jacobi<0, 1>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    if (idx == 1 && rotate) {
      apply_jacobi<2, 3>(A, c, s, c2, s2, cs);
      multiply_jacobi<2, 3>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    rotate = false;
    p = p_ids[idx+2];
    q = q_ids[idx+2];
    a_pq = A[p*4+q];
    if (fabs(a_pq) > 0) {
      rotate = true;
      const double a_pp = A[p*4+p];
      const double a_qq = A[q*4+q];
      compute_c_s(a_pq, a_pp, a_qq, c, s, c2, s2, cs);
    }
    COLVARS_SYNC_WARP;
    if (idx == 0 && rotate) {
      apply_jacobi<0, 2>(A, c, s, c2, s2, cs);
      multiply_jacobi<0, 2>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    if (idx == 1 && rotate) {
      apply_jacobi<1, 3>(A, c, s, c2, s2, cs);
      multiply_jacobi<1, 3>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    rotate = false;
    p = p_ids[idx+4];
    q = q_ids[idx+4];
    a_pq = A[p*4+q];
    if (fabs(a_pq) > 0) {
      rotate = true;
      const double a_pp = A[p*4+p];
      const double a_qq = A[q*4+q];
      compute_c_s(a_pq, a_pp, a_qq, c, s, c2, s2, cs);
    }
    COLVARS_SYNC_WARP;
    if (idx == 0 && rotate) {
      apply_jacobi<0, 3>(A, c, s, c2, s2, cs);
      multiply_jacobi<0, 3>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    if (idx == 1 && rotate) {
      apply_jacobi<1, 2>(A, c, s, c2, s2, cs);
      multiply_jacobi<1, 2>(V, c, s);
    }
    COLVARS_SYNC_WARP;
    off_diag_sum =
      fabs(A[0*4+1]) + fabs(A[0*4+2]) + fabs(A[0*4+3]) +
      fabs(A[1*4+2]) + fabs(A[1*4+3]) +
      fabs(A[2*4+3]);
    // Check the number of iterations
    ++iteration;
    if (iteration > JACOBI_MAX_ITERATION) {
      if (idx == 0 && max_reached) atomicAdd(max_reached, 1);
      break;
    }
  }
  // Sort
  double p;
  if (idx == 0) {
    int k;
    // We only require to sort the leading eigenvector to the first column
    k = 0;
    p = A[0*4+0];
    // Find the index of the leading eigenvalue
    #pragma unroll
    for (int i0 = 1; i0 < 4; ++i0) {
      if (A[i0*4+i0] > p) {
        p = A[i0*4+i0];
        k = i0;
      }
    }
    if (k != 0) {
      // Move the leading eigenvalue
      A[k*4+k] = A[0*4+0];
      A[0*4+0] = p;
      // Move the corresponding leading eigenvector
      #pragma unroll
      for (int j0 = 0; j0 < 4; ++j0) {
        p = V[j0*4+0];
        V[j0*4+0] = V[j0*4+k];
        V[j0*4+k] = p;
      }
    }
    A_in[0] = V[0];
    A_in[1] = V[4];
    A_in[2] = V[8];
    A_in[3] = V[12];
    A_in[4] = V[1];
    A_in[5] = V[5];
    A_in[6] = V[9];
    A_in[7] = V[13];
    A_in[8] = V[2];
    A_in[9] = V[6];
    A_in[10] = V[10];
    A_in[11] = V[14];
    A_in[12] = V[3];
    A_in[13] = V[7];
    A_in[14] = V[11];
    A_in[15] = V[15];
    eigvals[0] = A[0*4+0];
    eigvals[1] = A[1*4+1];
    eigvals[2] = A[2*4+2];
    eigvals[3] = A[3*4+3];
    // printf("eigvals = %lf %lf %lf %lf\n", eigvals[0], eigvals[1], eigvals[2], eigvals[3]);
    q->q0 = V[0];
    q->q1 = V[4];
    q->q2 = V[8];
    q->q3 = V[12];
    if (monitor_crossings) {
      if (q_old->norm2() > 0) {
        q->match(*q_old);
        if (q_old->inner(*q) < (1.0 - crossing_threshold)) {
          atomicAdd(discontinuous_rotation, 1);
        }
      }
    }
  }
}

#undef JACOBI_MAX_ITERATION
#undef JACOBI_TOLERANCE

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  void* args[] = {
    &S_eigvec, &S_eigval, &max_reached,
    &q, &monitor_crossings,
    &crossing_threshold, &q_old,
    &discontinuous_rotation};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)jacobi_4x4_kernel;
  kernelNodeParams.gridDim        = dim3(1, 1, 1);
  kernelNodeParams.blockDim       = dim3(2, 1, 1);
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

#elif defined(COLVARS_SYCL)
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
}
