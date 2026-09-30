#include "colvar_gpu_support.h"
#include "colvaratoms_kernel.h"
#include "colvartypes.h"

#if defined(COLVARS_CUDA)
#include <cub/block/block_reduce.cuh>
#include <cub/warp/warp_reduce.cuh>
#endif

#if defined (COLVARS_HIP)
#include <hipcub/block/block_reduce.hpp>
#include <hipcub/warp/warp_reduce.hpp>
#define cub hipcub
#endif

namespace colvars_gpu {
#if defined(COLVARS_CUDA) || defined(COLVARS_HIP)
__global__ void atoms_pos_from_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  const cvm::real* __restrict atoms_pos_x_proxy,
  const cvm::real* __restrict atoms_pos_y_proxy,
  const cvm::real* __restrict atoms_pos_z_proxy,
  cvm::real* __restrict atoms_pos_x_ag,
  cvm::real* __restrict atoms_pos_y_ag,
  cvm::real* __restrict atoms_pos_z_ag,
  unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    atoms_pos_x_ag[i] = atoms_pos_x_proxy[proxy_index];
    atoms_pos_y_ag[i] = atoms_pos_y_proxy[proxy_index];
    atoms_pos_z_ag[i] = atoms_pos_z_proxy[proxy_index];
  }
}

int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_proxy = atoms_pos_proxy;
  const cvm::real* atoms_pos_y_proxy = atoms_pos_x_proxy + proxy_stride;
  const cvm::real* atoms_pos_z_proxy = atoms_pos_y_proxy + proxy_stride;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_pos_x_proxy,
     &atoms_pos_y_proxy,
     &atoms_pos_z_proxy,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)atoms_pos_from_proxy_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (cvm::debug()) {
    cvm::log_static("Add atoms_pos_from_proxy node.\n");
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

int atoms_pos_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_pos_proxy,
  cvm::real* atoms_pos_ag,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaStream_t stream) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_pos_x_proxy = atoms_pos_proxy;
  const cvm::real* atoms_pos_y_proxy = atoms_pos_x_proxy + proxy_stride;
  const cvm::real* atoms_pos_z_proxy = atoms_pos_y_proxy + proxy_stride;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_pos_x_proxy,
     &atoms_pos_y_proxy,
     &atoms_pos_z_proxy,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &num_atoms};
  if (cvm::debug()) {
    cvm::log_static("Run atoms_pos_from_proxy.\n");
  }
  return checkGPUError(cudaLaunchKernel((void*)atoms_pos_from_proxy_kernel,
    num_blocks, block_size, args, 0, stream));
}

__global__ void change_one_coordinate_kernel(
  cvm::real* __restrict atoms_pos_ag,
  size_t array_id, cvm::real step_size) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i == 0) {
    atoms_pos_ag[array_id] += step_size;
  }
}

int change_one_coordinate(
  cvm::real* atoms_pos_ag, size_t atom_id_in_group, int xyz,
  cvm::real step_size, unsigned int num_atoms, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (xyz >= 0 && xyz < 3) {
    size_t array_id = num_atoms * xyz + atom_id_in_group;
    void* args[] = {&atoms_pos_ag, &array_id, &step_size};
    if (cvm::debug()) {
      cvm::log_static("Run change_one_coordinate.\n");
    }
    error_code |= checkGPUError(cudaLaunchKernel(
      (void*)change_one_coordinate_kernel,
      1, 1, args, 0, stream));
  }
  return error_code;
}

template <int BLOCK_SIZE, bool b_cog, bool b_com, bool save_cog>
__global__ void atoms_calc_cog_com_kernel(
  const cvm::real* __restrict atoms_mass,
  const cvm::real* __restrict atoms_pos_x_ag,
  const cvm::real* __restrict atoms_pos_y_ag,
  const cvm::real* __restrict atoms_pos_z_ag,
  cvm::rvector* __restrict cog_tmp,
  cvm::rvector* __restrict cog_out,
  cvm::rvector* __restrict cog_saved,
  cvm::rvector* __restrict com_tmp,
  cvm::rvector* __restrict com_out,
  cvm::rvector* __restrict h_cog_out,
  cvm::rvector* __restrict h_com_out,
  cvm::real total_mass,
  unsigned int* __restrict tbcount,
  unsigned int num_atoms) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::rvector cog{0, 0, 0};
  cvm::rvector com{0, 0, 0};
  while (i < num_atoms) {
    if (b_cog) {
      cog.x += atoms_pos_x_ag[i];
      cog.y += atoms_pos_y_ag[i];
      cog.z += atoms_pos_z_ag[i];
    }
    if (b_com) {
      com.x += atoms_mass[i] * atoms_pos_x_ag[i];
      com.y += atoms_mass[i] * atoms_pos_y_ag[i];
      com.z += atoms_mass[i] * atoms_pos_z_ag[i];
    }
    i += gridSize;
  }
  __syncthreads();
  typedef cub::BlockReduce<double, BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  if (b_cog) {
    cog.x = BlockReduce(temp_storage).Sum(cog.x); __syncthreads();
    cog.y = BlockReduce(temp_storage).Sum(cog.y); __syncthreads();
    cog.z = BlockReduce(temp_storage).Sum(cog.z); __syncthreads();
  }
  if (b_com) {
    com.x = BlockReduce(temp_storage).Sum(com.x); __syncthreads();
    com.y = BlockReduce(temp_storage).Sum(com.y); __syncthreads();
    com.z = BlockReduce(temp_storage).Sum(com.z); __syncthreads();
  }
  if (threadIdx.x == 0) {
    if (b_cog) {
      atomicAdd(&(cog_tmp->x), cog.x);
      atomicAdd(&(cog_tmp->y), cog.y);
      atomicAdd(&(cog_tmp->z), cog.z);
    }
    if (b_com) {
      atomicAdd(&(com_tmp->x), com.x);
      atomicAdd(&(com_tmp->y), com.y);
      atomicAdd(&(com_tmp->z), com.z);
    }
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    if (threadIdx.x == 0) {
      if (b_cog) {
        cog.x = cog_tmp->x / num_atoms;
        cog.y = cog_tmp->y / num_atoms;
        cog.z = cog_tmp->z / num_atoms;
        cog_out->x = cog.x;
        cog_out->y = cog.y;
        cog_out->z = cog.z;
        h_cog_out->x = cog.x;
        h_cog_out->y = cog.y;
        h_cog_out->z = cog.z;
        if (save_cog) {
          cog_saved->x = cog.x;
          cog_saved->y = cog.y;
          cog_saved->z = cog.z;
        }
        // Clear the temporary vector
        cog_tmp->x = 0;
        cog_tmp->y = 0;
        cog_tmp->z = 0;
        // printf("main = %d, (calc_cog) atom group = %p, COG = (%lf, %lf, %lf)\n", int(b_cog && b_com), atoms_pos_x_ag, cog.x, cog.y, cog.z);
      }
      if (b_com) {
        com.x = com_tmp->x / total_mass;
        com.y = com_tmp->y / total_mass;
        com.z = com_tmp->z / total_mass;
        com_out->x = com.x;
        com_out->y = com.y;
        com_out->z = com.z;
        com_tmp->x = 0;
        com_tmp->y = 0;
        com_tmp->z = 0;
        h_com_out->x = com.x;
        h_com_out->y = com.y;
        h_com_out->z = com.z;
      }
      // printf("tbcount = %p\n", (void*)tbcount);
      tbcount[0] = 0;
      __threadfence();
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  const cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  const cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  const cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] =
    {&atoms_mass,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &d_cog_tmp,
     &d_cog_out,
     &d_cog_origin,
     &d_com_tmp,
     &d_com_out,
     &h_cog_out,
     &h_com_out,
     &total_mass,
     &tbcount,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  if (d_cog_origin == nullptr) {
    kernelNodeParams.func =
      (void*)atoms_calc_cog_com_kernel<block_size, true, true, false>;
  } else {
    kernelNodeParams.func =
      (void*)atoms_calc_cog_com_kernel<block_size, true, true, true>;
  }
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
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  const cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  const cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  const cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  const cvm::real* atoms_mass = nullptr;
  const cvm::rvector* d_com_tmp = nullptr;
  const cvm::rvector* d_com_out = nullptr;
  const cvm::rvector* h_com_out = nullptr;
  cvm::real total_mass = 0.0;
  void* args[] =
    {&atoms_mass,
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &d_cog_tmp,
     &d_cog_out,
     &d_cog_origin,
     &d_com_tmp,
     &d_com_out,
     &h_cog_out,
     &h_com_out,
     &total_mass,
     &tbcount,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  if (d_cog_origin == nullptr) {
    kernelNodeParams.func = (void*)atoms_calc_cog_com_kernel<block_size, true, false, false>;
  } else {
    kernelNodeParams.func = (void*)atoms_calc_cog_com_kernel<block_size, true, false, true>;
  }
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

template <bool ag_rotate>
__global__ void atoms_total_force_from_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  const cvm::real* __restrict atoms_total_force_x_proxy,
  const cvm::real* __restrict atoms_total_force_y_proxy,
  const cvm::real* __restrict atoms_total_force_z_proxy,
  cvm::real* __restrict atoms_total_force_x_ag,
  cvm::real* __restrict atoms_total_force_y_ag,
  cvm::real* __restrict atoms_total_force_z_ag,
  const cvm::quaternion* __restrict q,
  unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  cvm::rmatrix rot_mat;
  if (ag_rotate) {
    rot_mat = q->rotation_matrix();
  }
  if (i < num_atoms) {
    const int proxy_index = atoms_proxy_index[i];
    const cvm::real fx = atoms_total_force_x_proxy[proxy_index];
    const cvm::real fy = atoms_total_force_y_proxy[proxy_index];
    const cvm::real fz = atoms_total_force_z_proxy[proxy_index];
    if (ag_rotate) {
      atoms_total_force_x_ag[i] = rot_mat.xx * fx + rot_mat.xy * fy + rot_mat.xz * fz;
      atoms_total_force_y_ag[i] = rot_mat.yx * fx + rot_mat.yy * fy + rot_mat.yz * fz;
      atoms_total_force_z_ag[i] = rot_mat.zx * fx + rot_mat.zy * fy + rot_mat.zz * fz;
    } else {
      atoms_total_force_x_ag[i] = fx;
      atoms_total_force_y_ag[i] = fy;
      atoms_total_force_z_ag[i] = fz;
    }
  }
}

int atoms_total_force_from_proxy(
  const int* atoms_proxy_index,
  const cvm::real* atoms_total_force_proxy,
  cvm::real* atoms_total_force_ag,
  bool rotate,
  const cvm::quaternion* q,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaStream_t stream) {
  if (num_atoms == 0) return COLVARS_OK;
  const unsigned int block_size = default_block_size;
  const unsigned int grid = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_total_force_x_proxy = atoms_total_force_proxy;
  const cvm::real* atoms_total_force_y_proxy = atoms_total_force_x_proxy + proxy_stride;
  const cvm::real* atoms_total_force_z_proxy = atoms_total_force_y_proxy + proxy_stride;
  cvm::real* atoms_total_force_x_ag = atoms_total_force_ag;
  cvm::real* atoms_total_force_y_ag = atoms_total_force_x_ag + num_atoms;
  cvm::real* atoms_total_force_z_ag = atoms_total_force_y_ag + num_atoms;
  void* args[] =
    {&atoms_proxy_index,
     &atoms_total_force_x_proxy,
     &atoms_total_force_y_proxy,
     &atoms_total_force_z_proxy,
     &atoms_total_force_x_ag,
     &atoms_total_force_y_ag,
     &atoms_total_force_z_ag,
     &q,
     &num_atoms};
  if (cvm::debug()) {
    cvm::log_static("Run " + cvm::to_str(__func__) + " kernel.\n");
  }
  if (rotate) {
    return checkGPUError(cudaLaunchKernel(
      (void*)atoms_total_force_from_proxy_kernel<true>,
      grid, block_size, args, 0, stream));
  } else {
    return checkGPUError(cudaLaunchKernel(
      (void*)atoms_total_force_from_proxy_kernel<false>,
      grid, block_size, args, 0, stream));
  }
}

template <bool ag_rotate>
__global__ void apply_colvar_force_to_proxy_kernel(
  const int* __restrict atoms_proxy_index,
  cvm::real* __restrict atoms_applied_force_x_proxy,
  cvm::real* __restrict atoms_applied_force_y_proxy,
  cvm::real* __restrict atoms_applied_force_z_proxy,
  const cvm::real* __restrict grad_x,
  const cvm::real* __restrict grad_y,
  const cvm::real* __restrict grad_z,
  cvm::real* __restrict force_ptr,
  const cvm::quaternion* __restrict q,
  unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const cvm::real force = (*force_ptr);
  cvm::rmatrix rot_inv;
  if (ag_rotate) {
    rot_inv = q->conjugate().rotation_matrix();
  }
  if (i < num_atoms) {
    const unsigned int proxy_index = atoms_proxy_index[i];
    cvm::real fx, fy, fz;
    if (ag_rotate) {
      fx = force * (rot_inv.xx * grad_x[i] +
                    rot_inv.xy * grad_y[i] +
                    rot_inv.xz * grad_z[i]);
      fy = force * (rot_inv.yx * grad_x[i] +
                    rot_inv.yy * grad_y[i] +
                    rot_inv.yz * grad_z[i]);
      fz = force * (rot_inv.zx * grad_x[i] +
                    rot_inv.zy * grad_y[i] +
                    rot_inv.zz * grad_z[i]);
    } else {
      fx = force * grad_x[i];
      fy = force * grad_y[i];
      fz = force * grad_z[i];
    }
    atomicAdd(&(atoms_applied_force_x_proxy[proxy_index]), fx);
    atomicAdd(&(atoms_applied_force_y_proxy[proxy_index]), fy);
    atomicAdd(&(atoms_applied_force_z_proxy[proxy_index]), fz);
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  // if (num_atoms == 0) return;
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_applied_force_x_proxy = atoms_applied_force_proxy;
  cvm::real* atoms_applied_force_y_proxy = atoms_applied_force_x_proxy + proxy_stride;
  cvm::real* atoms_applied_force_z_proxy = atoms_applied_force_y_proxy + proxy_stride;
  const cvm::real* grad_x = atoms_grad_ag;
  const cvm::real* grad_y = grad_x + num_atoms;
  const cvm::real* grad_z = grad_y + num_atoms;
  void* args[] = {
    &atoms_proxy_index,
    &atoms_applied_force_x_proxy,
    &atoms_applied_force_y_proxy,
    &atoms_applied_force_z_proxy,
    &grad_x,
    &grad_y,
    &grad_z,
    &colvar_force,
    &q,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  if (rotate) {
    kernelNodeParams.func = (void*)apply_colvar_force_to_proxy_kernel<true>;
  } else {
    kernelNodeParams.func = (void*)apply_colvar_force_to_proxy_kernel<false>;
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

int apply_fitting_colvar_force_to_proxy(
  const int* atoms_proxy_index,
  cvm::real* atoms_applied_force_proxy,
  const cvm::real* atoms_grad_ag,
  cvm::real* colvar_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_applied_force_x_proxy = atoms_applied_force_proxy;
  cvm::real* atoms_applied_force_y_proxy = atoms_applied_force_x_proxy + proxy_stride;
  cvm::real* atoms_applied_force_z_proxy = atoms_applied_force_y_proxy + proxy_stride;
  const cvm::real* grad_x = atoms_grad_ag;
  const cvm::real* grad_y = grad_x + num_atoms;
  const cvm::real* grad_z = grad_y + num_atoms;
  const cvm::quaternion* q = nullptr;
  void* args[] = {
    &atoms_proxy_index,
    &atoms_applied_force_x_proxy,
    &atoms_applied_force_y_proxy,
    &atoms_applied_force_z_proxy,
    &grad_x,
    &grad_y,
    &grad_z,
    &colvar_force,
    &q,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  kernelNodeParams.func = (void*)apply_colvar_force_to_proxy_kernel<false>;
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

/**
 * @brief Calculate the gradients of rotated positions wrt quaternion
 */
template <bool B_ag_center,
          bool B_ag_rotate,
          int BLOCK_SIZE>
__global__ void calc_fit_forces_impl_loop1_kernel(
  const cvm::real* __restrict atoms_grad_or_force_x,
  const cvm::real* __restrict atoms_grad_or_force_y,
  const cvm::real* __restrict atoms_grad_or_force_z,
  const cvm::real* __restrict atoms_pos_unrotated_x,
  const cvm::real* __restrict atoms_pos_unrotated_y,
  const cvm::real* __restrict atoms_pos_unrotated_z,
  const colvars_gpu::rotation_derivative_gpu* __restrict rot_deriv,
  const cvm::quaternion* __restrict q,
  const unsigned int main_group_size,
  const unsigned int group_for_fit_size,
  double3* __restrict atom_grad,
  double* __restrict sum_dxdq,
  cvm::rmatrix* __restrict dxdC,
  unsigned int* __restrict tbcount) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  __shared__ bool isLastBlockDone;
  if (threadIdx.x == 0) {
    isLastBlockDone = false;
  }
  __syncthreads();
  cvm::real C[3][3] = {{0}};
  cvm::rvector main_grad{0, 0, 0};
  while (i < main_group_size) {
    const cvm::rvector g{atoms_grad_or_force_x[i],
                         atoms_grad_or_force_y[i],
                         atoms_grad_or_force_z[i]};
    if (B_ag_center || B_ag_rotate) {
      main_grad += g;
    }
    if (B_ag_rotate) {
      C[0][0] += g.x * atoms_pos_unrotated_x[i];
      C[0][1] += g.x * atoms_pos_unrotated_y[i];
      C[0][2] += g.x * atoms_pos_unrotated_z[i];
      C[1][0] += g.y * atoms_pos_unrotated_x[i];
      C[1][1] += g.y * atoms_pos_unrotated_y[i];
      C[1][2] += g.y * atoms_pos_unrotated_z[i];
      C[2][0] += g.z * atoms_pos_unrotated_x[i];
      C[2][1] += g.z * atoms_pos_unrotated_y[i];
      C[2][2] += g.z * atoms_pos_unrotated_z[i];
    }
    i += gridSize;
  }
  __syncthreads();
  typedef cub::BlockReduce<cvm::real, BLOCK_SIZE> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  main_grad.x = BlockReduce(temp_storage).Sum(main_grad.x); __syncthreads();
  main_grad.y = BlockReduce(temp_storage).Sum(main_grad.y); __syncthreads();
  main_grad.z = BlockReduce(temp_storage).Sum(main_grad.z); __syncthreads();
  if (B_ag_rotate) {
    C[0][0] = BlockReduce(temp_storage).Sum(C[0][0]); __syncthreads();
    C[0][1] = BlockReduce(temp_storage).Sum(C[0][1]); __syncthreads();
    C[0][2] = BlockReduce(temp_storage).Sum(C[0][2]); __syncthreads();
    C[1][0] = BlockReduce(temp_storage).Sum(C[1][0]); __syncthreads();
    C[1][1] = BlockReduce(temp_storage).Sum(C[1][1]); __syncthreads();
    C[1][2] = BlockReduce(temp_storage).Sum(C[1][2]); __syncthreads();
    C[2][0] = BlockReduce(temp_storage).Sum(C[2][0]); __syncthreads();
    C[2][1] = BlockReduce(temp_storage).Sum(C[2][1]); __syncthreads();
    C[2][2] = BlockReduce(temp_storage).Sum(C[2][2]); __syncthreads();
  }
  if (threadIdx.x == 0) {
    if (B_ag_rotate) {
      colvars_gpu::array1d<cvm::real, 4> partial_dxdq;
      partial_dxdq =
        q->derivative_element_wise_product_sum<decltype(partial_dxdq)>(C);
      atomicAdd(&(sum_dxdq[0]), partial_dxdq[0]);
      atomicAdd(&(sum_dxdq[1]), partial_dxdq[1]);
      atomicAdd(&(sum_dxdq[2]), partial_dxdq[2]);
      atomicAdd(&(sum_dxdq[3]), partial_dxdq[3]);
    }
    atomicAdd(&(atom_grad->x), main_grad.x);
    atomicAdd(&(atom_grad->y), main_grad.y);
    atomicAdd(&(atom_grad->z), main_grad.z);
    __threadfence();
    unsigned int value = atomicInc(tbcount, gridDim.x);
    isLastBlockDone = (value == (gridDim.x - 1));
  }
  __syncthreads();
  if (isLastBlockDone) {
    // Compute dxdC in a single warp
    const unsigned int warpID = threadIdx.x / warpSize;
    if (warpID == 0) {
      const unsigned int tid = threadIdx.x;
      constexpr const int valid_items = 4;
      cvm::rmatrix dxdq_dqdC;
      dxdq_dqdC.reset();
      if (tid < valid_items) {
        dxdq_dqdC += rot_deriv->project_force_to_C_from_dxdqi(tid, sum_dxdq[tid]);
      }
      COLVARS_SYNC_WARP;
      using WarpReduce = cub::WarpReduce<cvm::real, valid_items>;
      __shared__ typename WarpReduce::TempStorage warp_temp_storage;
      dxdq_dqdC.xx = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.xx, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.xy = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.xy, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.xz = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.xz, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.yx = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.yx, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.yy = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.yy, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.yz = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.yz, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.zx = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.zx, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.zy = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.zy, valid_items); COLVARS_SYNC_WARP;
      dxdq_dqdC.zz = WarpReduce(warp_temp_storage).Sum(dxdq_dqdC.zz, valid_items); COLVARS_SYNC_WARP;
      if (tid == 0) {
        dxdC->xx = dxdq_dqdC.xx;
        dxdC->xy = dxdq_dqdC.xy;
        dxdC->xz = dxdq_dqdC.xz;
        dxdC->yx = dxdq_dqdC.yx;
        dxdC->yy = dxdq_dqdC.yy;
        dxdC->yz = dxdq_dqdC.yz;
        dxdC->zx = dxdq_dqdC.zx;
        dxdC->zy = dxdq_dqdC.zy;
        dxdC->zz = dxdq_dqdC.zz;
        // Clear the sum_dxdq array for the next reduction
        sum_dxdq[0] = 0;
        sum_dxdq[1] = 0;
        sum_dxdq[2] = 0;
        sum_dxdq[3] = 0;
      }
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      if (B_ag_center) {
        main_grad.x = atom_grad->x;
        main_grad.y = atom_grad->y;
        main_grad.z = atom_grad->z;
        if (B_ag_rotate) {
          const cvm::rmatrix rot_inv = q->conjugate().rotation_matrix();
          const cvm::real x = main_grad.x * rot_inv.xx + main_grad.y * rot_inv.xy + main_grad.z * rot_inv.xz;
          const cvm::real y = main_grad.x * rot_inv.yx + main_grad.y * rot_inv.yy + main_grad.z * rot_inv.yz;
          const cvm::real z = main_grad.x * rot_inv.zx + main_grad.y * rot_inv.zy + main_grad.z * rot_inv.zz;
          main_grad.x = x;
          main_grad.y = y;
          main_grad.z = z;
        }
        main_grad.x *= -1.0 / group_for_fit_size;
        main_grad.y *= -1.0 / group_for_fit_size;
        main_grad.z *= -1.0 / group_for_fit_size;
        atom_grad->x = main_grad.x;
        atom_grad->y = main_grad.y;
        atom_grad->z = main_grad.z;
      }
      tbcount[0] = 0;
    }
  }
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms_main + block_size - 1) / block_size;
  num_blocks = std::min(default_reduce_max_num_blocks, num_blocks);
  const cvm::real* pos_x = pos_unrotated;
  const cvm::real* pos_y = pos_x + num_atoms_main;
  const cvm::real* pos_z = pos_y + num_atoms_main;
  const cvm::real* grad_x = main_grad;
  const cvm::real* grad_y = grad_x + num_atoms_main;
  const cvm::real* grad_z = grad_y + num_atoms_main;
  void* args[] = {
    &grad_x, &grad_y, &grad_z,
    &pos_x, &pos_y, &pos_z,
    &rot_deriv, &q,
    &num_atoms_main, &num_atoms_fitting,
    &atom_grad, &sum_dxdq, &dxdC, &tbcount};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<true, true, block_size>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<true, false, block_size>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop1_kernel<false, true, block_size>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

// loop 2: iterate over the fitting group
template <bool B_ag_center,
          bool B_ag_rotate,
          int BLOCK_SIZE,
          bool do_fit_gradients>
__global__ void calc_fit_forces_impl_loop2_kernel(
  cvm::real* __restrict fit_grad_x,
  cvm::real* __restrict fit_grad_y,
  cvm::real* __restrict fit_grad_z,
  const colvars_gpu::rotation_derivative_gpu* __restrict rot_deriv,
  const double3* __restrict atom_grad,
  const cvm::rmatrix* __restrict dxdC,
  const int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const unsigned int group_for_fit_size) {
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned int gridSize = blockDim.x * gridDim.x;
  while (i < group_for_fit_size) {
    cvm::rvector fitting_force_grad{0, 0, 0};
    if (B_ag_center) {
      fitting_force_grad.x += atom_grad->x;
      fitting_force_grad.y += atom_grad->y;
      fitting_force_grad.z += atom_grad->z;
    }
    if (B_ag_rotate) {
      fitting_force_grad += rot_deriv->project_force_to_group1(i, *dxdC);
    }
    // accessor_fitting(i, fitting_force_grad);
    if (do_fit_gradients) {
      fit_grad_x[i] = fitting_force_grad.x;
      fit_grad_y[i] = fitting_force_grad.y;
      fit_grad_z[i] = fitting_force_grad.z;
    } else {
      const int pid = atoms_proxy_index[i];
      atomicAdd(&(proxy_new_force_x[pid]), fitting_force_grad.x);
      atomicAdd(&(proxy_new_force_y[pid]), fitting_force_grad.y);
      atomicAdd(&(proxy_new_force_z[pid]), fitting_force_grad.z);
    }
    i += gridSize;
  }
}

int calc_fit_gradients_impl_loop2(
  cvm::real* fit_grad,
  const colvars_gpu::rotation_derivative_gpu* rot_deriv,
  const double3* atom_grad,
  const cvm::rmatrix* dxdC,
  unsigned int group_for_fit_size,
  bool ag_center, bool ag_rotate,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (group_for_fit_size + block_size - 1) / block_size;
  cvm::real* fit_grad_x = fit_grad;
  cvm::real* fit_grad_y = fit_grad_x + group_for_fit_size;
  cvm::real* fit_grad_z = fit_grad_y + group_for_fit_size;
  const int* atoms_proxy_index = nullptr;
  cvm::real* proxy_new_force_x = nullptr;
  cvm::real* proxy_new_force_y = nullptr;
  cvm::real* proxy_new_force_z = nullptr;
  void* args[] = {
    &fit_grad_x, &fit_grad_y, &fit_grad_z,
    &rot_deriv, &atom_grad, &dxdC,
    &atoms_proxy_index,
    &proxy_new_force_x,
    &proxy_new_force_y,
    &proxy_new_force_z,
    &group_for_fit_size};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, true, block_size, true>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, false, block_size, true>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<false, true, block_size, true>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void apply_translation_kernel(
  cvm::real* __restrict atoms_pos_x_ag,
  cvm::real* __restrict atoms_pos_y_ag,
  cvm::real* __restrict atoms_pos_z_ag,
  cvm::real translation_vector_factor,
  const cvm::rvector* __restrict translation_vector,
  unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    atoms_pos_x_ag[i] += translation_vector_factor * translation_vector->x;
    atoms_pos_y_ag[i] += translation_vector_factor * translation_vector->y;
    atoms_pos_z_ag[i] += translation_vector_factor * translation_vector->z;
  }
}

int apply_translation(
  cvm::real* atoms_pos_ag,
  cvm::real translation_vector_factor,
  const cvm::rvector* translation_vector,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] = {
     &atoms_pos_x_ag,
     &atoms_pos_y_ag,
     &atoms_pos_z_ag,
     &translation_vector_factor,
     &translation_vector,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_translation_kernel;
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
    cvm::log_static("x_ptr = " + cvm::to_str((void*)atoms_pos_x_ag) + "\n");
    cvm::log_static("y_ptr = " + cvm::to_str((void*)atoms_pos_y_ag) + "\n");
    cvm::log_static("z_ptr = " + cvm::to_str((void*)atoms_pos_z_ag) + "\n");
    cvm::log_static("pos = " + cvm::to_str((void*)translation_vector) + "\n");
    cvm::log_static("factor = " + cvm::to_str(translation_vector_factor) + "\n");
    cvm::log_static("num_atoms = " + cvm::to_str((int)num_atoms) + "\n");
  }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

__global__ void rotate_with_quaternion_kernel(
  cvm::real* __restrict pos_x,
  cvm::real* __restrict pos_y,
  cvm::real* __restrict pos_z,
  const cvm::quaternion* __restrict q, int num_atoms) {
  const auto rot_mat = q->rotation_matrix();
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  while (i < num_atoms) {
  // if (i < num_atoms) {
    const cvm::real new_x = rot_mat.xx * pos_x[i] +
                            rot_mat.xy * pos_y[i] +
                            rot_mat.xz * pos_z[i];
    const cvm::real new_y = rot_mat.yx * pos_x[i] +
                            rot_mat.yy * pos_y[i] +
                            rot_mat.yz * pos_z[i];
    const cvm::real new_z = rot_mat.zx * pos_x[i] +
                            rot_mat.zy * pos_y[i] +
                            rot_mat.zz * pos_z[i];
    pos_x[i] = new_x;
    pos_y[i] = new_y;
    pos_z[i] = new_z;
    i += gridSize;
  }
}

int rotate_with_quaternion(
  cvm::real* atoms_pos_ag,
  const cvm::quaternion* q,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  cvm::real* atoms_pos_x_ag = atoms_pos_ag;
  cvm::real* atoms_pos_y_ag = atoms_pos_x_ag + num_atoms;
  cvm::real* atoms_pos_z_ag = atoms_pos_y_ag + num_atoms;
  void* args[] = {
    &atoms_pos_x_ag,
    &atoms_pos_y_ag,
    &atoms_pos_z_ag,
    &q, &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)rotate_with_quaternion_kernel;
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

__global__ void accumulate_cpu_force_kernel(
  const cvm::real* __restrict h_atoms_force_x,
  const cvm::real* __restrict h_atoms_force_y,
  const cvm::real* __restrict h_atoms_force_z,
  cvm::real* __restrict d_atoms_force_x,
  cvm::real* __restrict d_atoms_force_y,
  cvm::real* __restrict d_atoms_force_z,
  const unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    // Because the same atom group may be referenced by multiple CVCs,
    // so I have to use atomicAdd to accumulate the forces
    atomicAdd(&(d_atoms_force_x[i]), h_atoms_force_x[i]);
    atomicAdd(&(d_atoms_force_y[i]), h_atoms_force_y[i]);
    atomicAdd(&(d_atoms_force_z[i]), h_atoms_force_z[i]);
  }
}

int accumulate_cpu_force(
  const cvm::real* h_atoms_force,
  cvm::real* d_atoms_force,
  unsigned int num_atoms,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const cvm::real* h_atoms_force_x = h_atoms_force;
  const cvm::real* h_atoms_force_y = h_atoms_force_x + num_atoms;
  const cvm::real* h_atoms_force_z = h_atoms_force_y + num_atoms;
  cvm::real* d_atoms_force_x = d_atoms_force;
  cvm::real* d_atoms_force_y = d_atoms_force_x + num_atoms;
  cvm::real* d_atoms_force_z = d_atoms_force_y + num_atoms;
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  void* args[] = {
    &h_atoms_force_x,
    &h_atoms_force_y,
    &h_atoms_force_z,
    &d_atoms_force_x,
    &d_atoms_force_y,
    &d_atoms_force_z,
    &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)accumulate_cpu_force_kernel;
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

__global__ void apply_force_with_inverse_rotation_kernel(
  const cvm::real* __restrict atoms_force_x,
  const cvm::real* __restrict atoms_force_y,
  const cvm::real* __restrict atoms_force_z,
  const cvm::quaternion* q,
  int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const unsigned int num_atoms) {
  const cvm::rmatrix rot_inv = q->conjugate().rotation_matrix();
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const cvm::rvector f_ia{
      rot_inv.xx * atoms_force_x[i] +
      rot_inv.xy * atoms_force_y[i] +
      rot_inv.xz * atoms_force_z[i],
      rot_inv.yx * atoms_force_x[i] +
      rot_inv.yy * atoms_force_y[i] +
      rot_inv.yz * atoms_force_z[i],
      rot_inv.zx * atoms_force_x[i] +
      rot_inv.zy * atoms_force_y[i] +
      rot_inv.zz * atoms_force_z[i],
    };
    const int pid = atoms_proxy_index[i];
    atomicAdd(&(proxy_new_force_x[pid]), f_ia.x);
    atomicAdd(&(proxy_new_force_y[pid]), f_ia.y);
    atomicAdd(&(proxy_new_force_z[pid]), f_ia.z);
  }
}

int apply_force_with_inverse_rotation(
  const cvm::real* atoms_force,
  const cvm::quaternion* q,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_force_x = atoms_force;
  const cvm::real* atoms_force_y = atoms_force_x + num_atoms;
  const cvm::real* atoms_force_z = atoms_force_y + num_atoms;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] =
    {&atoms_force_x,
     &atoms_force_y,
     &atoms_force_z,
     &q,
     &atoms_proxy_index,
     &proxy_new_force_x,
     &proxy_new_force_y,
     &proxy_new_force_z,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_force_with_inverse_rotation_kernel;
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

__global__ void apply_force_kernel(
  const cvm::real* __restrict atoms_force_x,
  const cvm::real* __restrict atoms_force_y,
  const cvm::real* __restrict atoms_force_z,
  int* __restrict atoms_proxy_index,
  cvm::real* __restrict proxy_new_force_x,
  cvm::real* __restrict proxy_new_force_y,
  cvm::real* __restrict proxy_new_force_z,
  const unsigned int num_atoms) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < num_atoms) {
    const int pid = atoms_proxy_index[i];
    atomicAdd(&(proxy_new_force_x[pid]), atoms_force_x[i]);
    atomicAdd(&(proxy_new_force_y[pid]), atoms_force_y[i]);
    atomicAdd(&(proxy_new_force_z[pid]), atoms_force_z[i]);
  }
}

int apply_force(
  const cvm::real* atoms_force,
  const int* atoms_proxy_index,
  cvm::real* proxy_new_force,
  unsigned int num_atoms,
  unsigned int proxy_stride,
  cudaGraphNode_t& node,
  cudaGraph_t& graph,
  const std::vector<cudaGraphNode_t>& dependencies) {
  const int block_size = default_block_size;
  const int num_blocks = (num_atoms + block_size - 1) / block_size;
  const cvm::real* atoms_force_x = atoms_force;
  const cvm::real* atoms_force_y = atoms_force_x + num_atoms;
  const cvm::real* atoms_force_z = atoms_force_y + num_atoms;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] =
    {&atoms_force_x,
     &atoms_force_y,
     &atoms_force_z,
     &atoms_proxy_index,
     &proxy_new_force_x,
     &proxy_new_force_y,
     &proxy_new_force_z,
     &num_atoms};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.func           = (void*)apply_force_kernel;
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
  const std::vector<cudaGraphNode_t>& dependencies) {
  if (cvm::debug()) {
    cvm::log_static("Add " + cvm::to_str(__func__) + " node.\n");
  }
  return calc_fit_gradients_impl_loop1(
    pos_unrotated, main_force, rot_deriv, q, num_atoms_main,
    num_atoms_fitting, atom_grad, sum_dxdq, dxdC, tbcount,
    ag_center, ag_rotate, node, graph, dependencies);
}

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
  const std::vector<cudaGraphNode_t>& dependencies) {
  const unsigned int block_size = default_block_size;
  const unsigned int num_blocks = (group_for_fit_size + block_size - 1) / block_size;
  cvm::real* fit_grad_x = nullptr;
  cvm::real* fit_grad_y = nullptr;
  cvm::real* fit_grad_z = nullptr;
  cvm::real* proxy_new_force_x = proxy_new_force;
  cvm::real* proxy_new_force_y = proxy_new_force_x + proxy_stride;
  cvm::real* proxy_new_force_z = proxy_new_force_y + proxy_stride;
  void* args[] = {
    &fit_grad_x, &fit_grad_y, &fit_grad_z,
    &rot_deriv, &atom_grad, &dxdC,
    &atoms_proxy_index,
    &proxy_new_force_x,
    &proxy_new_force_y,
    &proxy_new_force_z,
    &group_for_fit_size};
  cudaKernelNodeParams kernelNodeParams = {0};
  kernelNodeParams.gridDim        = dim3(num_blocks, 1, 1);
  kernelNodeParams.blockDim       = dim3(block_size, 1, 1);
  kernelNodeParams.sharedMemBytes = 0;
  kernelNodeParams.kernelParams   = args;
  kernelNodeParams.extra          = NULL;
  if (ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, true, block_size, false>;
  }
  if (ag_center && !ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<true, false, block_size, false>;
  }
  if (!ag_center && ag_rotate) {
    kernelNodeParams.func =
      (void*)calc_fit_forces_impl_loop2_kernel<false, true, block_size, false>;
  }
  if (!ag_center && !ag_rotate) {
    return COLVARS_OK;
  }
  // cvmodule pointer not available here
  // if (cvmodule->debug()) {
  //   cvmodule->log("Add " + cvm::to_str(__func__) + " node.\n");
  // }
  return checkGPUError(cudaGraphAddKernelNode(
    &node, graph, dependencies.data(),
    dependencies.size(), &kernelNodeParams));
}

#elif defined(COLVARS_SYCL)
#endif // defined(COLVARS_CUDA) || defined(COLVARS_HIP)
};
