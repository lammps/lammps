/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDAEXEC_HPP
#define KOKKOS_CUDAEXEC_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <mutex>
#include <string>
#include <cstdint>
#include <cmath>
#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>
#include <Cuda/Kokkos_Cuda_Locks.hpp>
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
// Maximize L1 cache and minimize shared memory:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferL1 );
// For 2.0 capability: 48 KB L1 and 16 KB shared
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
  if (cuda_instance->m_maxShmemPerBlock < shmem) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("CudaParallelLaunch (or graph node creation) FAILED: shared"
                    " memory request is too large"));
  }
}

// This function needs to be template on DriverType and LaunchBounds
// so that the static bool is unique for each type combo
// KernelFuncPtr does not necessarily contain that type information.
template <class DriverType, class LaunchBounds, class KernelFuncPtr>
inline void configure_shmem_preference(KernelFuncPtr const& func,
                                       bool prefer_shmem) {
#ifndef KOKKOS_ARCH_KEPLER
  // On Kepler the L1 has no benefit since it doesn't cache reads
  auto set_cache_config = [&] {
    CUDA_SAFE_CALL(cudaFuncSetCacheConfig(
        func,
        (prefer_shmem ? cudaFuncCachePreferShared : cudaFuncCachePreferL1)));
    return prefer_shmem;
  };
  static bool cache_config_preference_cached = set_cache_config();
  if (cache_config_preference_cached != prefer_shmem) {
    cache_config_preference_cached = set_cache_config();
  }
#else
  // Use the parameters so we don't get a warning
  (void)func;
  (void)prefer_shmem;
#endif
}

template <class Policy>
std::enable_if_t<Policy::experimental_contains_desired_occupancy>
modify_launch_configuration_if_desired_occupancy_is_specified(
    Policy const& policy, cudaDeviceProp const& properties,
    cudaFuncAttributes const& attributes, dim3 const& block, int& shmem,
    bool& prefer_shmem) {
  int const block_size        = block.x * block.y * block.z;
  int const desired_occupancy = policy.impl_get_desired_occupancy().value();

  size_t const shmem_per_sm_prefer_l1 = get_shmem_per_sm_prefer_l1(properties);
  size_t const static_shmem           = attributes.sharedSizeBytes;

  // round to nearest integer and avoid division by zero
  int active_blocks = std::max(
      1, static_cast<int>(std::round(
             static_cast<double>(properties.maxThreadsPerMultiProcessor) /
             block_size * desired_occupancy / 100)));
  int const dynamic_shmem =
      shmem_per_sm_prefer_l1 / active_blocks - static_shmem;

  if (dynamic_shmem > shmem) {
    shmem        = dynamic_shmem;
    prefer_shmem = false;
  }
}

template <class Policy>
std::enable_if_t<!Policy::experimental_contains_desired_occupancy>
modify_launch_configuration_if_desired_occupancy_is_specified(
    Policy const&, cudaDeviceProp const&, cudaFuncAttributes const&,
    dim3 const& /*block*/, int& /*shmem*/, bool& /*prefer_shmem*/) {}

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
    (base_t::
         get_kernel_func())<<<grid, block, shmem, cuda_instance->m_stream>>>(
        driver);
  }

#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance, bool prefer_shmem) {
    //----------------------------------------
    auto const& graph = Impl::get_cuda_graph_from_kernel(driver);
    KOKKOS_EXPECTS(bool(graph));
    auto& graph_node = Impl::get_cuda_graph_node_from_kernel(driver);
    // Expect node not yet initialized
    KOKKOS_EXPECTS(!bool(graph_node));

    if (!Impl::is_empty_launch(grid, block)) {
      Impl::check_shmem_request(cuda_instance, shmem);
      Impl::configure_shmem_preference<DriverType, LaunchBounds>(
          base_t::get_kernel_func(), prefer_shmem);

      void const* args[] = {&driver};

      cudaKernelNodeParams params = {};

      params.blockDim       = block;
      params.gridDim        = grid;
      params.sharedMemBytes = shmem;
      params.func           = (void*)base_t::get_kernel_func();
      params.kernelParams   = (void**)args;
      params.extra          = nullptr;

      CUDA_SAFE_CALL(cudaGraphAddKernelNode(
          &graph_node, graph, /* dependencies = */ nullptr,
          /* numDependencies = */ 0, &params));
    } else {
      // We still need an empty node for the dependency structure
      CUDA_SAFE_CALL(cudaGraphAddEmptyNode(&graph_node, graph,
                                           /* dependencies = */ nullptr,
                                           /* numDependencies = */ 0));
    }
    KOKKOS_ENSURES(bool(graph_node))
  }
#endif
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

    cudaMemcpyAsync(driver_ptr, &driver, sizeof(DriverType), cudaMemcpyDefault,
                    cuda_instance->m_stream);
    (base_t::
         get_kernel_func())<<<grid, block, shmem, cuda_instance->m_stream>>>(
        driver_ptr);
  }

#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance, bool prefer_shmem) {
    //----------------------------------------
    auto const& graph = Impl::get_cuda_graph_from_kernel(driver);
    KOKKOS_EXPECTS(bool(graph));
    auto& graph_node = Impl::get_cuda_graph_node_from_kernel(driver);
    // Expect node not yet initialized
    KOKKOS_EXPECTS(!bool(graph_node));

    if (!Impl::is_empty_launch(grid, block)) {
      Impl::check_shmem_request(cuda_instance, shmem);
      Impl::configure_shmem_preference<DriverType, LaunchBounds>(
          base_t::get_kernel_func(), prefer_shmem);

      auto* driver_ptr = Impl::allocate_driver_storage_for_kernel(driver);

      // Unlike in the non-graph case, we can get away with doing an async copy
      // here because the `DriverType` instance is held in the GraphNodeImpl
      // which is guaranteed to be alive until the graph instance itself is
      // destroyed, where there should be a fence ensuring that the allocation
      // associated with this kernel on the device side isn't deleted.
      cudaMemcpyAsync(driver_ptr, &driver, sizeof(DriverType),
                      cudaMemcpyDefault, cuda_instance->m_stream);

      void const* args[] = {&driver_ptr};

      cudaKernelNodeParams params = {};

      params.blockDim       = block;
      params.gridDim        = grid;
      params.sharedMemBytes = shmem;
      params.func           = (void*)base_t::get_kernel_func();
      params.kernelParams   = (void**)args;
      params.extra          = nullptr;

      CUDA_SAFE_CALL(cudaGraphAddKernelNode(
          &graph_node, graph, /* dependencies = */ nullptr,
          /* numDependencies = */ 0, &params));
    } else {
      // We still need an empty node for the dependency structure
      CUDA_SAFE_CALL(cudaGraphAddEmptyNode(&graph_node, graph,
                                           /* dependencies = */ nullptr,
                                           /* numDependencies = */ 0));
    }
    KOKKOS_ENSURES(bool(graph_node))
  }
#endif
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
    // Wait until the previous kernel that uses the constant buffer is done
    CUDA_SAFE_CALL(cudaEventSynchronize(cuda_instance->constantMemReusable));

    // Copy functor (synchronously) to staging buffer in pinned host memory
    unsigned long* staging = cuda_instance->constantMemHostStaging;
    memcpy(staging, &driver, sizeof(DriverType));

    // Copy functor asynchronously from there to constant memory on the device
    cudaMemcpyToSymbolAsync(kokkos_impl_cuda_constant_memory_buffer, staging,
                            sizeof(DriverType), 0, cudaMemcpyHostToDevice,
                            cudaStream_t(cuda_instance->m_stream));

    // Invoke the driver function on the device
    (base_t::
         get_kernel_func())<<<grid, block, shmem, cuda_instance->m_stream>>>();

    // Record an event that says when the constant buffer can be reused
    CUDA_SAFE_CALL(cudaEventRecord(cuda_instance->constantMemReusable,
                                   cudaStream_t(cuda_instance->m_stream)));
  }

#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
  inline static void create_parallel_launch_graph_node(
      DriverType const& driver, dim3 const& grid, dim3 const& block, int shmem,
      CudaInternal const* cuda_instance, bool prefer_shmem) {
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
        driver, grid, block, shmem, cuda_instance, prefer_shmem);
  }
#endif
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
                                   const CudaInternal* cuda_instance,
                                   bool prefer_shmem) {
    if (!Impl::is_empty_launch(grid, block)) {
      // Prevent multiple threads to simultaneously set the cache configuration
      // preference and launch the same kernel
      static std::mutex mutex;
      std::lock_guard<std::mutex> lock(mutex);

      Impl::check_shmem_request(cuda_instance, shmem);

      // If a desired occupancy is specified, we compute how much shared memory
      // to ask for to achieve that occupancy, assuming that the cache
      // configuration is `cudaFuncCachePreferL1`.  If the amount of dynamic
      // shared memory computed is actually smaller than `shmem` we overwrite
      // `shmem` and set `prefer_shmem` to `false`.
      modify_launch_configuration_if_desired_occupancy_is_specified(
          driver.get_policy(), cuda_instance->m_deviceProp,
          get_cuda_func_attributes(), block, shmem, prefer_shmem);

      Impl::configure_shmem_preference<
          DriverType, Kokkos::LaunchBounds<MaxThreadsPerBlock, MinBlocksPerSM>>(
          base_t::get_kernel_func(), prefer_shmem);

      KOKKOS_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();

      // Invoke the driver function on the device
      base_t::invoke_kernel(driver, grid, block, shmem, cuda_instance);

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
      CUDA_SAFE_CALL(cudaGetLastError());
      cuda_instance->fence();
#endif
    }
  }

  static cudaFuncAttributes get_cuda_func_attributes() {
    // Race condition inside of cudaFuncGetAttributes if the same address is
    // given requires using a local variable as input instead of a static Rely
    // on static variable initialization to make sure only one thread executes
    // the code and the result is visible.
    auto wrap_get_attributes = []() -> cudaFuncAttributes {
      cudaFuncAttributes attr_tmp;
      CUDA_SAFE_CALL(
          cudaFuncGetAttributes(&attr_tmp, base_t::get_kernel_func()));
      return attr_tmp;
    };
    static cudaFuncAttributes attr = wrap_get_attributes();
    return attr;
  }
};

// </editor-fold> end CudaParallelLaunchImpl }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="CudaParallelLaunch"> {{{1

template <class DriverType, class LaunchBounds = Kokkos::LaunchBounds<>,
          Experimental::CudaLaunchMechanism LaunchMechanism =
              DeduceCudaLaunchMechanism<DriverType>::launch_mechanism,
          bool DoGraph = DriverType::Policy::is_graph_kernel::value
#ifndef KOKKOS_CUDA_ENABLE_GRAPHS
                         && false
#endif
          >
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

#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
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
#endif

// </editor-fold> end CudaParallelLaunch }}}1
//==============================================================================

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDAEXEC_HPP */
