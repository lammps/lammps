//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_CUDAEXEC_HPP
#define KOKKOS_CUDAEXEC_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <mutex>
#include <cstdint>
#include <cmath>
#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>
#include <Cuda/Kokkos_Cuda_GraphNodeKernel.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief  Access to constant memory on the device */
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE

__device__ __constant__ extern unsigned long
    kokkos_impl_cuda_constant_memory_buffer[];

#else

__device__ __constant__ unsigned long kokkos_impl_cuda_constant_memory_buffer
    [Kokkos::Impl::CudaTraits::ConstantMemoryUsage / sizeof(unsigned long)];

#endif

template <typename T>
inline __device__ T* kokkos_impl_cuda_shared_memory() {
  extern __shared__ Kokkos::CudaSpace::size_type sh[];
  return (T*)sh;
}

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------

template <class DriverType>
__global__ static void cuda_parallel_launch_constant_memory() {
  const DriverType& driver =
      *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB, minBperSM) static void cuda_parallel_launch_constant_memory() {
  const DriverType& driver =
      *((const DriverType*)kokkos_impl_cuda_constant_memory_buffer);

  driver();
}

template <class DriverType>
__global__ static void cuda_parallel_launch_local_memory(
    const DriverType driver) {
  driver();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void cuda_parallel_launch_local_memory(const DriverType
                                                                 driver) {
  driver();
}

template <class DriverType>
__global__ static void cuda_parallel_launch_global_memory(
    const DriverType* driver) {
  driver->operator()();
}

template <class DriverType, unsigned int maxTperB, unsigned int minBperSM>
__global__ __launch_bounds__(
    maxTperB,
    minBperSM) static void cuda_parallel_launch_global_memory(const DriverType*
                                                                  driver) {
  driver->operator()();
}

//==============================================================================
// <editor-fold desc="Some helper functions for launch code readability"> {{{1

inline bool is_empty_launch(dim3 const& grid, dim3 const& block) {
  return (grid.x == 0) || ((block.x * block.y * block.z) == 0);
}

inline void check_shmem_request(CudaInternal const* cuda_instance, int shmem) {
  int const maxShmemPerBlock = cuda_instance->m_deviceProp.sharedMemPerBlock;
  if (maxShmemPerBlock < shmem) {
    Kokkos::Impl::throw_runtime_exception(
        "CudaParallelLaunch (or graph node creation) FAILED: shared memory "
        "request is too large");
  }
}

// These functions need to be templated on DriverType and LaunchBounds
// so that the static bool is unique for each type combo
// KernelFuncPtr does not necessarily contain that type information.
template <class DriverType, class LaunchBounds, class KernelFuncPtr>
const cudaFuncAttributes& get_cuda_kernel_func_attributes(
    int cuda_device, const KernelFuncPtr& func) {
  // Only call cudaFuncGetAttributes once for each unique kernel
  // by leveraging static variable initialization rules
  static std::map<int, cudaFuncAttributes> func_attr;
  if (func_attr.find(cuda_device) == func_attr.end()) {
    cudaFuncAttributes attr;
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaSetDevice(cuda_device));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFuncGetAttributes(&attr, func));
    func_attr.emplace(cuda_device, attr);
  }
  return func_attr[cuda_device];
}

template <class DriverType, class LaunchBounds, class KernelFuncPtr>
inline void configure_shmem_preference(const int cuda_device,
                                       const KernelFuncPtr& func,
                                       const cudaDeviceProp& device_props,
                                       const size_t block_size, int& shmem,
                                       const size_t occupancy) {
#ifndef KOKKOS_ARCH_KEPLER

  const auto& func_attr =
      get_cuda_kernel_func_attributes<DriverType, LaunchBounds>(cuda_device,
                                                                func);

  // Compute limits for number of blocks due to registers/SM
  const size_t regs_per_sm     = device_props.regsPerMultiprocessor;
  const size_t regs_per_thread = func_attr.numRegs;
  // The granularity of register allocation is chunks of 256 registers per warp
  // -> 8 registers per thread
  const size_t allocated_regs_per_thread = 8 * ((regs_per_thread + 8 - 1) / 8);
  size_t max_blocks_regs =
      regs_per_sm / (allocated_regs_per_thread * block_size);

  // Compute the maximum number of warps as a function of the number of
  // registers
  const size_t max_warps_per_sm_registers =
      cuda_max_warps_per_sm_registers(device_props, func_attr);

  // Correct the number of blocks to respect the maximum number of warps per
  // SM, which is constrained to be a multiple of the warp allocation
  // granularity defined in `cuda_warp_per_sm_allocation_granularity`.
  while ((max_blocks_regs * block_size / device_props.warpSize) >
         max_warps_per_sm_registers)
    max_blocks_regs--;

  // Compute how many threads per sm we actually want
  const size_t max_threads_per_sm = device_props.maxThreadsPerMultiProcessor;
  // only allocate multiples of warp size
  const size_t num_threads_desired =
      ((max_threads_per_sm * occupancy / 100 + 31) / 32) * 32;
  // Get close to the desired occupancy,
  // don't undershoot by much but also don't allocate a whole new block just
  // because one is a few threads over otherwise.
  size_t num_blocks_desired =
      (num_threads_desired + block_size * 0.8) / block_size;
  num_blocks_desired = ::std::min(max_blocks_regs, num_blocks_desired);
  if (num_blocks_desired == 0) num_blocks_desired = 1;

  // Calculate how much shared memory we need per block
  size_t shmem_per_block = shmem + func_attr.sharedSizeBytes;

  // The minimum shared memory allocation we can have in total per SM is 8kB.
  // If we want to lower occupancy we have to make sure we request at least that
  // much in aggregate over all blocks, so that shared memory actually becomes a
  // limiting factor for occupancy
  constexpr size_t min_shmem_size_per_sm = 8192;
  if ((occupancy < 100) &&
      (shmem_per_block * num_blocks_desired < min_shmem_size_per_sm)) {
    shmem_per_block = min_shmem_size_per_sm / num_blocks_desired;
    // Need to set the caller's shmem variable so that the
    // kernel launch uses the correct dynamic shared memory request
    shmem = shmem_per_block - func_attr.sharedSizeBytes;
  }

  // Compute the carveout fraction we need based on occupancy
  // Use multiples of 8kB
  const size_t max_shmem_per_sm = device_props.sharedMemPerMultiprocessor;
  size_t carveout               = shmem_per_block == 0
                        ? 0
                        : 100 *
                              (((num_blocks_desired * shmem_per_block +
                                 min_shmem_size_per_sm - 1) /
                                min_shmem_size_per_sm) *
                               min_shmem_size_per_sm) /
                              max_shmem_per_sm;
  if (carveout > 100) carveout = 100;

  // Set the carveout, but only call it once per kernel or when it changes
  // FIXME_CUDA_MULTIPLE_DEVICES
  auto set_cache_config = [&] {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (CudaInternal::singleton().cuda_func_set_attribute_wrapper(
            func, cudaFuncAttributePreferredSharedMemoryCarveout, carveout)));
    return carveout;
  };
  // Store the value in a static variable so we only reset if needed
  static size_t cache_config_preference_cached = set_cache_config();
  if (cache_config_preference_cached != carveout) {
    cache_config_preference_cached = set_cache_config();
  }
#else
  // Use the parameters so we don't get a warning
  (void)func;
  (void)device_props;
  (void)block_size;
  (void)occupancy;
#endif
}

// </editor-fold> end Some helper functions for launch code readability }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="DeduceCudaLaunchMechanism"> {{{2

// Use local memory up to ConstantMemoryUseThreshold
// Use global memory above ConstantMemoryUsage
// In between use ConstantMemory

template <class DriverType>
struct DeduceCudaLaunchMechanism {
  constexpr static const Kokkos::Experimental::WorkItemProperty::
      HintLightWeight_t light_weight =
          Kokkos::Experimental::WorkItemProperty::HintLightWeight;
  constexpr static const Kokkos::Experimental::WorkItemProperty::
      HintHeavyWeight_t heavy_weight =
          Kokkos::Experimental::WorkItemProperty::HintHeavyWeight;
  constexpr static const typename DriverType::Policy::work_item_property
      property = typename DriverType::Policy::work_item_property();

  static constexpr const Experimental::CudaLaunchMechanism
      valid_launch_mechanism =
          // BuildValidMask
      (sizeof(DriverType) < CudaTraits::KernelArgumentLimit
           ? Experimental::CudaLaunchMechanism::LocalMemory
           : Experimental::CudaLaunchMechanism::Default) |
      (sizeof(DriverType) < CudaTraits::ConstantMemoryUsage
           ? Experimental::CudaLaunchMechanism::ConstantMemory
           : Experimental::CudaLaunchMechanism::Default) |
      Experimental::CudaLaunchMechanism::GlobalMemory;

  static constexpr const Experimental::CudaLaunchMechanism
      requested_launch_mechanism =
          (((property & light_weight) == light_weight)
               ? Experimental::CudaLaunchMechanism::LocalMemory
               : Experimental::CudaLaunchMechanism::ConstantMemory) |
          Experimental::CudaLaunchMechanism::GlobalMemory;

  static constexpr const Experimental::CudaLaunchMechanism
      default_launch_mechanism =
          // BuildValidMask
      (sizeof(DriverType) < CudaTraits::ConstantMemoryUseThreshold)
          ? Experimental::CudaLaunchMechanism::LocalMemory
          : ((sizeof(DriverType) < CudaTraits::ConstantMemoryUsage)
                 ? Experimental::CudaLaunchMechanism::ConstantMemory
                 : Experimental::CudaLaunchMechanism::GlobalMemory);

  //              None                LightWeight    HeavyWeight
  // F<UseT       LCG LCG L  L        LCG  LG L  L    LCG  CG L  C
  // UseT<F<KAL   LCG LCG C  C        LCG  LG C  L    LCG  CG C  C
  // Kal<F<CMU     CG LCG C  C         CG  LG C  G     CG  CG C  C
  // CMU<F          G LCG G  G          G  LG G  G      G  CG G  G
  static constexpr const Experimental::CudaLaunchMechanism launch_mechanism =
      ((property & light_weight) == light_weight)
          ? (sizeof(DriverType) < CudaTraits::KernelArgumentLimit
                 ? Experimental::CudaLaunchMechanism::LocalMemory
                 : Experimental::CudaLaunchMechanism::GlobalMemory)
          : (((property & heavy_weight) == heavy_weight)
                 ? (sizeof(DriverType) < CudaTraits::ConstantMemoryUsage
                        ? Experimental::CudaLaunchMechanism::ConstantMemory
                        : Experimental::CudaLaunchMechanism::GlobalMemory)
                 : (default_launch_mechanism));
};

// </editor-fold> end DeduceCudaLaunchMechanism }}}2
//==============================================================================

//==============================================================================
// <editor-fold desc="CudaParallelLaunchKernelInvoker"> {{{1

// Base classes that summarize the differences between the different launch
// mechanisms

template <class DriverType, class LaunchBounds,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunchKernelFunc;

template <class DriverType, class LaunchBounds,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunchKernelInvoker;

//------------------------------------------------------------------------------
// <editor-fold desc="Local memory"> {{{2

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::LocalMemory> {
  static std::decay_t<decltype(cuda_parallel_launch_local_memory<
                               DriverType, MaxThreadsPerBlock, MinBlocksPerSM>)>
  get_kernel_func() {
    return cuda_parallel_launch_local_memory<DriverType, MaxThreadsPerBlock,
                                             MinBlocksPerSM>;
  }
};

template <class DriverType>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<0, 0>,
    Experimental::CudaLaunchMechanism::LocalMemory> {
  static std::decay_t<decltype(cuda_parallel_launch_local_memory<DriverType>)>
  get_kernel_func() {
    return cuda_parallel_launch_local_memory<DriverType>;
  }
};

//------------------------------------------------------------------------------

template <class DriverType, class LaunchBounds>
struct CudaParallelLaunchKernelInvoker<
    DriverType, LaunchBounds, Experimental::CudaLaunchMechanism::LocalMemory>
    : CudaParallelLaunchKernelFunc<
          DriverType, LaunchBounds,
          Experimental::CudaLaunchMechanism::LocalMemory> {
  using base_t = CudaParallelLaunchKernelFunc<
      DriverType, LaunchBounds, Experimental::CudaLaunchMechanism::LocalMemory>;
  static_assert(sizeof(DriverType) < CudaTraits::KernelArgumentLimit,
                "Kokkos Error: Requested CudaLaunchLocalMemory with a Functor "
                "larger than 4096 bytes.");

  static void invoke_kernel(DriverType const& driver, dim3 const& grid,
                            dim3 const& block, int shmem,
                            CudaInternal const* cuda_instance) {
    (base_t::get_kernel_func())<<<grid, block, shmem,
                                  cuda_instance->get_stream()>>>(driver);
  }

  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance) {
    //----------------------------------------
    auto const& graph = Impl::get_cuda_graph_from_kernel(driver);
    KOKKOS_EXPECTS(bool(graph));
    auto& graph_node = Impl::get_cuda_graph_node_from_kernel(driver);
    // Expect node not yet initialized
    KOKKOS_EXPECTS(!bool(graph_node));

    if (!Impl::is_empty_launch(grid, block)) {
      Impl::check_shmem_request(cuda_instance, shmem);
      if constexpr (DriverType::Policy::
                        experimental_contains_desired_occupancy) {
        int desired_occupancy =
            driver.get_policy().impl_get_desired_occupancy().value();
        size_t block_size = block.x * block.y * block.z;
        Impl::configure_shmem_preference<DriverType, LaunchBounds>(
            cuda_instance->m_cudaDev, base_t::get_kernel_func(),
            cuda_instance->m_deviceProp, block_size, shmem, desired_occupancy);
      }

      void const* args[] = {&driver};

      cudaKernelNodeParams params = {};

      params.blockDim       = block;
      params.gridDim        = grid;
      params.sharedMemBytes = shmem;
      params.func           = (void*)base_t::get_kernel_func();
      params.kernelParams   = (void**)args;
      params.extra          = nullptr;

      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (cuda_instance->cuda_graph_add_kernel_node_wrapper(
              &graph_node, graph, /* dependencies = */ nullptr,
              /* numDependencies = */ 0, &params)));
    } else {
      // We still need an empty node for the dependency structure
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (cuda_instance->cuda_graph_add_empty_node_wrapper(
              &graph_node, graph,
              /* dependencies = */ nullptr,
              /* numDependencies = */ 0)));
    }
    KOKKOS_ENSURES(bool(graph_node))
  }
};

// </editor-fold> end local memory }}}2
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// <editor-fold desc="Global Memory"> {{{2

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::GlobalMemory> {
  static void* get_kernel_func() {
    return cuda_parallel_launch_global_memory<DriverType, MaxThreadsPerBlock,
                                              MinBlocksPerSM>;
  }
};

template <class DriverType>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<0, 0>,
    Experimental::CudaLaunchMechanism::GlobalMemory> {
  static std::decay_t<decltype(cuda_parallel_launch_global_memory<DriverType>)>
  get_kernel_func() {
    return cuda_parallel_launch_global_memory<DriverType>;
  }
};

//------------------------------------------------------------------------------

template <class DriverType, class LaunchBounds>
struct CudaParallelLaunchKernelInvoker<
    DriverType, LaunchBounds, Experimental::CudaLaunchMechanism::GlobalMemory>
    : CudaParallelLaunchKernelFunc<
          DriverType, LaunchBounds,
          Experimental::CudaLaunchMechanism::GlobalMemory> {
  using base_t = CudaParallelLaunchKernelFunc<
      DriverType, LaunchBounds,
      Experimental::CudaLaunchMechanism::GlobalMemory>;

  static void invoke_kernel(DriverType const& driver, dim3 const& grid,
                            dim3 const& block, int shmem,
                            CudaInternal const* cuda_instance) {
    DriverType* driver_ptr = reinterpret_cast<DriverType*>(
        cuda_instance->scratch_functor(sizeof(DriverType)));

    KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_instance->cuda_memcpy_async_wrapper(
        driver_ptr, &driver, sizeof(DriverType), cudaMemcpyDefault)));
    (base_t::get_kernel_func())<<<grid, block, shmem,
                                  cuda_instance->get_stream()>>>(driver_ptr);
  }

  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance) {
    //----------------------------------------
    auto const& graph = Impl::get_cuda_graph_from_kernel(driver);
    KOKKOS_EXPECTS(bool(graph));
    auto& graph_node = Impl::get_cuda_graph_node_from_kernel(driver);
    // Expect node not yet initialized
    KOKKOS_EXPECTS(!bool(graph_node));

    if (!Impl::is_empty_launch(grid, block)) {
      Impl::check_shmem_request(cuda_instance, shmem);
      if constexpr (DriverType::Policy::
                        experimental_contains_desired_occupancy) {
        int desired_occupancy =
            driver.get_policy().impl_get_desired_occupancy().value();
        size_t block_size = block.x * block.y * block.z;
        Impl::configure_shmem_preference<DriverType, LaunchBounds>(
            cuda_instance->m_cudaDev, base_t::get_kernel_func(),
            cuda_instance->m_deviceProp, block_size, shmem, desired_occupancy);
      }

      auto* driver_ptr = Impl::allocate_driver_storage_for_kernel(driver);

      // Unlike in the non-graph case, we can get away with doing an async copy
      // here because the `DriverType` instance is held in the GraphNodeImpl
      // which is guaranteed to be alive until the graph instance itself is
      // destroyed, where there should be a fence ensuring that the allocation
      // associated with this kernel on the device side isn't deleted.
      KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_instance->cuda_memcpy_async_wrapper(
          driver_ptr, &driver, sizeof(DriverType), cudaMemcpyDefault)));

      void const* args[] = {&driver_ptr};

      cudaKernelNodeParams params = {};

      params.blockDim       = block;
      params.gridDim        = grid;
      params.sharedMemBytes = shmem;
      params.func           = (void*)base_t::get_kernel_func();
      params.kernelParams   = (void**)args;
      params.extra          = nullptr;

      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (cuda_instance->cuda_graph_add_kernel_node_wrapper(
              &graph_node, graph, /* dependencies = */ nullptr,
              /* numDependencies = */ 0, &params)));
    } else {
      // We still need an empty node for the dependency structure
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (cuda_instance->cuda_graph_add_empty_node_wrapper(
              &graph_node, graph,
              /* dependencies = */ nullptr,
              /* numDependencies = */ 0)));
    }
    KOKKOS_ENSURES(bool(graph_node))
  }
};

// </editor-fold> end Global Memory }}}2
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// <editor-fold desc="Constant Memory"> {{{2

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    Experimental::CudaLaunchMechanism::ConstantMemory> {
  static std::decay_t<decltype(cuda_parallel_launch_constant_memory<
                               DriverType, MaxThreadsPerBlock, MinBlocksPerSM>)>
  get_kernel_func() {
    return cuda_parallel_launch_constant_memory<DriverType, MaxThreadsPerBlock,
                                                MinBlocksPerSM>;
  }
};

template <class DriverType>
struct CudaParallelLaunchKernelFunc<
    DriverType, Kokkos::LaunchBounds<0, 0>,
    Experimental::CudaLaunchMechanism::ConstantMemory> {
  static std::decay_t<
      decltype(cuda_parallel_launch_constant_memory<DriverType>)>
  get_kernel_func() {
    return cuda_parallel_launch_constant_memory<DriverType>;
  }
};

//------------------------------------------------------------------------------

template <class DriverType, class LaunchBounds>
struct CudaParallelLaunchKernelInvoker<
    DriverType, LaunchBounds, Experimental::CudaLaunchMechanism::ConstantMemory>
    : CudaParallelLaunchKernelFunc<
          DriverType, LaunchBounds,
          Experimental::CudaLaunchMechanism::ConstantMemory> {
  using base_t = CudaParallelLaunchKernelFunc<
      DriverType, LaunchBounds,
      Experimental::CudaLaunchMechanism::ConstantMemory>;
  static_assert(sizeof(DriverType) < CudaTraits::ConstantMemoryUsage,
                "Kokkos Error: Requested CudaLaunchConstantMemory with a "
                "Functor larger than 32kB.");

  static void invoke_kernel(DriverType const& driver, dim3 const& grid,
                            dim3 const& block, int shmem,
                            CudaInternal const* cuda_instance) {
    int cuda_device = cuda_instance->m_cudaDev;
    // Wait until the previous kernel that uses the constant buffer is done
    std::lock_guard<std::mutex> lock(
        CudaInternal::constantMemMutexPerDevice[cuda_device]);
    KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_instance->cuda_event_synchronize_wrapper(
        CudaInternal::constantMemReusablePerDevice[cuda_device])));

    // Copy functor (synchronously) to staging buffer in pinned host memory
    unsigned long* staging =
        cuda_instance->constantMemHostStagingPerDevice[cuda_device];
    memcpy(staging, &driver, sizeof(DriverType));

    // Copy functor asynchronously from there to constant memory on the device
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (cuda_instance->cuda_memcpy_to_symbol_async_wrapper(
            kokkos_impl_cuda_constant_memory_buffer, staging,
            sizeof(DriverType), 0, cudaMemcpyHostToDevice)));

    // Invoke the driver function on the device
    (base_t::get_kernel_func())<<<grid, block, shmem,
                                  cuda_instance->get_stream()>>>();

    // Record an event that says when the constant buffer can be reused
    KOKKOS_IMPL_CUDA_SAFE_CALL((cuda_instance->cuda_event_record_wrapper(
        CudaInternal::constantMemReusablePerDevice[cuda_device])));
  }

  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance) {
    // Just use global memory; coordinating through events to share constant
    // memory with the non-graph interface is not really reasonable since
    // events don't work with Graphs directly, and this would anyway require
    // a much more complicated structure that finds previous nodes in the
    // dependency structure of the graph and creates an implicit dependence
    // based on the need for constant memory (which we would then have to
    // somehow go and prove was not creating a dependency cycle, and I don't
    // even know if there's an efficient way to do that, let alone in the
    // structure we currenty have).
    using global_launch_impl_t = CudaParallelLaunchKernelInvoker<
        DriverType, LaunchBounds,
        Experimental::CudaLaunchMechanism::GlobalMemory>;
    global_launch_impl_t::create_parallel_launch_graph_node(
        driver, grid, block, shmem, cuda_instance);
  }
};

// </editor-fold> end Constant Memory }}}2
//------------------------------------------------------------------------------

// </editor-fold> end CudaParallelLaunchKernelInvoker }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="CudaParallelLaunchImpl"> {{{1

template <class DriverType, class LaunchBounds,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunchImpl;

template <class DriverType, unsigned int MaxThreadsPerBlock,
          unsigned int MinBlocksPerSM,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunchImpl<
    DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
    LaunchMechanism>
    : CudaParallelLaunchKernelInvoker<
          DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
          LaunchMechanism> {
  using base_t = CudaParallelLaunchKernelInvoker<
      DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>,
      LaunchMechanism>;

  inline static void launch_kernel(const DriverType& driver, const dim3& grid,
                                   const dim3& block, int shmem,
                                   const CudaInternal* cuda_instance) {
    if (!Impl::is_empty_launch(grid, block)) {
      // Prevent multiple threads to simultaneously set the cache configuration
      // preference and launch the same kernel
      static std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);

      Impl::check_shmem_request(cuda_instance, shmem);

      if constexpr (DriverType::Policy::
                        experimental_contains_desired_occupancy) {
        int desired_occupancy =
            driver.get_policy().impl_get_desired_occupancy().value();
        size_t block_size = block.x * block.y * block.z;
        Impl::configure_shmem_preference<
            DriverType,
            Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>>(
            cuda_instance->m_cudaDev, base_t::get_kernel_func(),
            cuda_instance->m_deviceProp, block_size, shmem, desired_occupancy);
      }

      desul::ensure_cuda_lock_arrays_on_device();

      // Invoke the driver function on the device
      base_t::invoke_kernel(driver, grid, block, shmem, cuda_instance);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetLastError());
      cuda_instance->fence(
          "Kokkos::Impl::launch_kernel: Debug Only Check for Execution Error");
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes(int cuda_device) {
    return get_cuda_kernel_func_attributes<
        DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>>(
        cuda_device, base_t::get_kernel_func());
  }
};

// </editor-fold> end CudaParallelLaunchImpl }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="CudaParallelLaunch"> {{{1

template <class DriverType, class LaunchBounds = Kokkos::LaunchBounds<>,
          Experimental::CudaLaunchMechanism LaunchMechanism =
              DeduceCudaLaunchMechanism<DriverType>::launch_mechanism,
          bool DoGraph = DriverType::Policy::is_graph_kernel::value>
struct CudaParallelLaunch;

// General launch mechanism
template <class DriverType, class LaunchBounds,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunch<DriverType, LaunchBounds, LaunchMechanism,
                          /* DoGraph = */ false>
    : CudaParallelLaunchImpl<DriverType, LaunchBounds, LaunchMechanism> {
  using base_t =
      CudaParallelLaunchImpl<DriverType, LaunchBounds, LaunchMechanism>;
  template <class... Args>
  CudaParallelLaunch(Args&&... args) {
    base_t::launch_kernel((Args &&) args...);
  }
};

// Launch mechanism for creating graph nodes
template <class DriverType, class LaunchBounds,
          Experimental::CudaLaunchMechanism LaunchMechanism>
struct CudaParallelLaunch<DriverType, LaunchBounds, LaunchMechanism,
                          /* DoGraph = */ true>
    : CudaParallelLaunchImpl<DriverType, LaunchBounds, LaunchMechanism> {
  using base_t =
      CudaParallelLaunchImpl<DriverType, LaunchBounds, LaunchMechanism>;
  template <class... Args>
  CudaParallelLaunch(Args&&... args) {
    base_t::create_parallel_launch_graph_node((Args &&) args...);
  }
};

// </editor-fold> end CudaParallelLaunch }}}1
//==============================================================================

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDAEXEC_HPP */
