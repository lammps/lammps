#include "colvar_gpu_support.h"
#include "colvarcomp_distance_kernel.h"
#include "colvar_rotation_derivative.h"
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
__global__ void calc_value_rmsd_kernel(
  const cvm::real* __restrict ref_pos_x,
  const cvm::real* __restrict ref_pos_y,
  const cvm::real* __restrict ref_pos_z,
  const cvm::real* __restrict pos_x,
  const cvm::real* __restrict pos_y,
  const cvm::real* __restrict pos_z,
  cvm::real* __restrict d_msd,
  cvm::real* __restrict h_msd,
  const unsigned int num_atoms,
  unsigned int* __restrict tbcount) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned int gridSize = blockDim.x * gridDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::real msd = 0;
  while (i < num_atoms) {
    const cvm::real diff_x = pos_x[i] - ref_pos_x[i];
    const cvm::real diff_y = pos_y[i] - ref_pos_y[i];
    const cvm::real diff_z = pos_z[i] - ref_pos_z[i];
    msd += diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    i += gridSize;
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  msd = BlockReduce(temp_storage).Sum(msd); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(d_msd, msd);
    __threadfence();
    unsigned int value = atomicInc(&(tbcount[0]), gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      *h_msd = *d_msd;
      *tbcount = 0;
      *d_msd = 0;
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  if (num_atoms == 0) return COLVARS_OK;
  int error_code = COLVARS_OK;
  const cvm::real* ag_pos_x = ag_pos;
  const cvm::real* ag_pos_y = ag_pos_x + num_atoms;
  const cvm::real* ag_pos_z = ag_pos_y + num_atoms;
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  for (unsigned int ip = 0; ip < num_permutations; ++ip) {
    const cvm::real* ref_pos_x = ref_pos + ip * num_atoms;
    const cvm::real* ref_pos_y = ref_pos_x + num_ref_pos;
    const cvm::real* ref_pos_z = ref_pos_y + num_ref_pos;
    cvm::real* d_permutation_msd = d_permutation_msds + ip;
    cvm::real* h_permutation_msd = h_permutation_msds + ip;
    unsigned int* tbcount = tbcounts + ip;
    void* args[] = {
      &ref_pos_x,
      &ref_pos_y,
      &ref_pos_z,
      &ag_pos_x,
      &ag_pos_y,
      &ag_pos_z,
      &d_permutation_msd,
      &h_permutation_msd,
      &num_atoms,
      &tbcount};
    cudaKernelNodeParams kernelNodeParams = {0};
    kernelNodeParams.func           = (void*)calc_value_rmsd_kernel<block_size>;
    kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
    kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
    kernelNodeParams.sharedMemBytes = 0;
    kernelNodeParams.kernelParams   = args;
    kernelNodeParams.extra          = NULL;
    cudaGraphNode_t node;
    error_code |= checkGPUError(cudaGraphAddKernelNode(
      &node, graph, dependencies.data(),
      dependencies.size(), &kernelNodeParams));
    nodes[ip] = node;
  }
  return error_code;
}

__global__ void calc_gradients_rmsd_kernel(
  const cvm::real* __restrict h_rmsd,
  const size_t* __restrict h_best_perm_index,
  const cvm::real* __restrict ref_pos_x,
  const cvm::real* __restrict ref_pos_y,
  const cvm::real* __restrict ref_pos_z,
  const cvm::real* __restrict pos_x,
  const cvm::real* __restrict pos_y,
  const cvm::real* __restrict pos_z,
  cvm::real* __restrict grad_x,
  cvm::real* __restrict grad_y,
  cvm::real* __restrict grad_z,
  unsigned int num_atoms) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned int gridSize = blockDim.x * gridDim.x;
  const cvm::real drmsddx2 = (*h_rmsd) > 0 ? 0.5 / ((*h_rmsd) * num_atoms) : 0;
  const size_t start = (*h_best_perm_index) * num_atoms;
  while (i < num_atoms) {
    // NOTE: I want the atom group be reusable by multiple CVCs so I use atomicAdd
    atomicAdd(&grad_x[i], drmsddx2 * 2.0 * (pos_x[i] - ref_pos_x[i + start]));
    atomicAdd(&grad_y[i], drmsddx2 * 2.0 * (pos_y[i] - ref_pos_y[i + start]));
    atomicAdd(&grad_z[i], drmsddx2 * 2.0 * (pos_z[i] - ref_pos_z[i + start]));
    i += gridSize;
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  if (num_atoms == 0) return COLVARS_OK;
  const cvm::real* ag_pos_x = ag_pos;
  const cvm::real* ag_pos_y = ag_pos_x + num_atoms;
  const cvm::real* ag_pos_z = ag_pos_y + num_atoms;
  cvm::real* grad_x = grad;
  cvm::real* grad_y = grad_x + num_atoms;
  cvm::real* grad_z = grad_y + num_atoms;
  const cvm::real* ref_pos_x = ref_pos;
  const cvm::real* ref_pos_y = ref_pos_x + num_ref_pos;
  const cvm::real* ref_pos_z = ref_pos_y + num_ref_pos;
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  void* args[] = {
    &h_rmsd, &h_best_perm_index,
    &ref_pos_x,
    &ref_pos_y,
    &ref_pos_z,
    &ag_pos_x,
    &ag_pos_y,
    &ag_pos_z,
    &grad_x,
    &grad_y,
    &grad_z,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)calc_gradients_rmsd_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

template <bool ag_rotate, int BLOCK_SIZE>
__global__ void calc_force_invgrads_rmsd_kernel(
  const int* __restrict atoms_proxy_index,
  const cvm::real* __restrict atoms_total_force_x_proxy,
  const cvm::real* __restrict atoms_total_force_y_proxy,
  const cvm::real* __restrict atoms_total_force_z_proxy,
  const cvm::quaternion* __restrict q,
  const cvm::real* __restrict grad_x,
  const cvm::real* __restrict grad_y,
  const cvm::real* __restrict grad_z,
  cvm::real* __restrict d_ft,
  cvm::real* __restrict h_ft,
  unsigned int num_atoms,
  unsigned int* tbcount) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned int gridSize = blockDim.x * gridDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::rmatrix rot_mat;
  if (ag_rotate) {
    rot_mat = q->rotation_matrix();
  }
  cvm::real sum_ft = 0;
  // fused kernel of reading total forces and computing ft
  while (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    cvm::real fx = atoms_total_force_x_proxy[proxy_index];
    cvm::real fy = atoms_total_force_y_proxy[proxy_index];
    cvm::real fz = atoms_total_force_z_proxy[proxy_index];
    if (ag_rotate) {
      const cvm::real old_fx = fx;
      const cvm::real old_fy = fy;
      const cvm::real old_fz = fz;
      fx = rot_mat.xx * old_fx + rot_mat.xy * old_fy + rot_mat.xz * old_fz;
      fy = rot_mat.yx * old_fx + rot_mat.yy * old_fy + rot_mat.yz * old_fz;
      fz = rot_mat.zx * old_fx + rot_mat.zy * old_fy + rot_mat.zz * old_fz;
    }
    sum_ft += grad_x[i] * fx + grad_y[i] * fy + grad_z[i] * fz;
    i += gridSize;
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  sum_ft = BlockReduce(temp_storage).Sum(sum_ft); __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(d_ft, sum_ft);
    __threadfence();
    unsigned int value = atomicInc(&(tbcount[0]), gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      *h_ft = (*d_ft) * num_atoms;
      *tbcount = 0;
      *d_ft = 0;
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  if (num_atoms == 0) return COLVARS_OK;
  const cvm::real* grad_x = grad;
  const cvm::real* grad_y = grad_x + num_atoms;
  const cvm::real* grad_z = grad_y + num_atoms;
  const cvm::real* atoms_total_force_x_proxy = atoms_total_force_proxy;
  const cvm::real* atoms_total_force_y_proxy = atoms_total_force_x_proxy + proxy_stride;
  const cvm::real* atoms_total_force_z_proxy = atoms_total_force_y_proxy + proxy_stride;
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  void* args[] = {
    &atoms_proxy_index,
    &atoms_total_force_x_proxy,
    &atoms_total_force_y_proxy,
    &atoms_total_force_z_proxy,
    &q,
    &grad_x,
    &grad_y,
    &grad_z,
    &d_ft,
    &h_ft,
    &num_atoms,
    &tbcount};
  cudaKernelNodeParams kernelNodeParams = {0};
  if (rotate) {
    kernelNodeParams.func           = (void*)calc_force_invgrads_rmsd_kernel<true, block_size>;
  } else {
    kernelNodeParams.func           = (void*)calc_force_invgrads_rmsd_kernel<false, block_size>;
  }
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

template <bool ag_center, bool ag_rotate, int BLOCK_SIZE>
__global__ void calc_Jacobian_derivative_rmsd_kernel(
  const size_t* __restrict h_best_perm_index,
  const cvm::real* __restrict ref_pos_x,
  const cvm::real* __restrict ref_pos_y,
  const cvm::real* __restrict ref_pos_z,
  const cvm::quaternion* __restrict q_ptr,
  const colvars_gpu::rotation_derivative_gpu* __restrict rot_deriv,
  const cvm::real* __restrict h_rmsd,
  cvm::real* __restrict d_jd,
  cvm::real* __restrict h_jd,
  unsigned int num_atoms,
  unsigned int* tbcount) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned int gridSize = blockDim.x * gridDim.x;
  const cvm::real rmsd = *h_rmsd;
  const cvm::real translation_term = ag_center ? 3.0 : 0.0;
  cvm::real rotation_term = 0.0;
  if (ag_rotate) {
    __shared__ bool isLastBlockDone;
    if (threadIdx.x == 0) {
      isLastBlockDone = false;
    }
    __syncthreads();
    cvm::real q[4] = {(*q_ptr)[0], (*q_ptr)[1], (*q_ptr)[2], (*q_ptr)[3]};
    // gradients of products of 2 quaternion components
    cvm::rvector dq[4];
    cvm::rvector grad_rot_mat[3][3];
    const size_t start = (*h_best_perm_index) * num_atoms;
    while (i < num_atoms) {
      // Gradient of optimal quaternion wrt current Cartesian position
      rot_deriv->calc_derivative_wrt_group1<false, true>(i, nullptr, dq);
      const cvm::rvector g11 = 2.0 * q[1] * dq[1];
      const cvm::rvector g22 = 2.0 * q[2] * dq[2];
      const cvm::rvector g33 = 2.0 * q[3] * dq[3];
      const cvm::rvector g01 = q[0] * dq[1] + q[1] * dq[0];
      const cvm::rvector g02 = q[0] * dq[2] + q[2] * dq[0];
      const cvm::rvector g03 = q[0] * dq[3] + q[3] * dq[0];
      const cvm::rvector g12 = q[1] * dq[2] + q[2] * dq[1];
      const cvm::rvector g13 = q[1] * dq[3] + q[3] * dq[1];
      const cvm::rvector g23 = q[2] * dq[3] + q[3] * dq[2];
      // Gradient of the rotation matrix wrt current Cartesian position
      grad_rot_mat[0][0] = -2.0 * (g22 + g33);
      grad_rot_mat[1][0] =  2.0 * (g12 + g03);
      grad_rot_mat[2][0] =  2.0 * (g13 - g02);
      grad_rot_mat[0][1] =  2.0 * (g12 - g03);
      grad_rot_mat[1][1] = -2.0 * (g11 + g33);
      grad_rot_mat[2][1] =  2.0 * (g01 + g23);
      grad_rot_mat[0][2] =  2.0 * (g02 + g13);
      grad_rot_mat[1][2] =  2.0 * (g23 - g01);
      grad_rot_mat[2][2] = -2.0 * (g11 + g22);
      const cvm::atom_pos y(
        ref_pos_x[i + start],
        ref_pos_y[i + start],
        ref_pos_z[i + start]);
      #pragma unroll
      for (size_t alpha = 0; alpha < 3; alpha++) {
        #pragma unroll
        for (size_t beta = 0; beta < 3; beta++) {
          rotation_term += grad_rot_mat[beta][alpha][alpha] * y[beta];
        // Note: equation was derived for inverse rotation (see colvars paper)
        // so here the matrix is transposed
        // (eq would give   divergence += grad_rot_mat[alpha][beta][alpha] * y[beta];)
        }
      }
      i += gridSize;
    }
    // Block reduction: sum the rotation term
    __syncthreads();
    typedef cub::BlockReduce<double, BLOCK_SIZE> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;
    rotation_term = BlockReduce(temp_storage).Sum(rotation_term); __syncthreads();
    // Reuse the d_jd buffer to sum the rotation_term
    if (threadIdx.x == 0) {
      atomicAdd(d_jd, rotation_term);
      __threadfence();
      unsigned int value = atomicInc(&(tbcount[0]), gridDim.x);
      isLastBlockDone = (value == (gridDim.x - 1));
    }
    __syncthreads();
    if (isLastBlockDone) {
      if (threadIdx.x == 0) {
        rotation_term = *d_jd;
        *d_jd = rmsd > 0 ? (3.0 * num_atoms - 1.0 - translation_term - rotation_term) / rmsd : 0.0;
        *h_jd = *d_jd;
        *d_jd = 0;
      }
    }
  } else {
    if (i == 0) {
      *d_jd = rmsd > 0 ? (3.0 * num_atoms - 1.0 - translation_term - rotation_term) / rmsd : 0.0;
      *h_jd = *d_jd;
      *d_jd = 0;
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  int error_code = COLVARS_OK;
  if (num_atoms == 0) return COLVARS_OK;
  const cvm::real* ref_pos_x = ref_pos;
  const cvm::real* ref_pos_y = ref_pos_x + num_ref_pos;
  const cvm::real* ref_pos_z = ref_pos_y + num_ref_pos;
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  void* args[] = {
    &h_best_perm_index,
    &ref_pos_x,
    &ref_pos_y,
    &ref_pos_z,
    &q,
    &rot_deriv,
    &h_rmsd,
    &d_jd,
    &h_jd,
    &num_atoms,
    &tbcount};
  cudaKernelNodeParams kernelNodeParams = {0};
  if (rotate) {
    if (center) {
      kernelNodeParams.func = (void*)calc_Jacobian_derivative_rmsd_kernel<true, true, block_size>;
    } else {
      kernelNodeParams.func = (void*)calc_Jacobian_derivative_rmsd_kernel<false, true, block_size>;
    }
    kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
    kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  } else {
    if (center) {
      kernelNodeParams.func = (void*)calc_Jacobian_derivative_rmsd_kernel<true, false, block_size>;
    } else {
      kernelNodeParams.func = (void*)calc_Jacobian_derivative_rmsd_kernel<false, false, block_size>;
    }
    // We need only a single thread to compute the translation_term
    kernelNodeParams.gridDim        = dim3(1, 1, 1);
    kernelNodeParams.blockDim       = dim3(1, 1, 1);
  }
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
  return error_code;
}
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
}
