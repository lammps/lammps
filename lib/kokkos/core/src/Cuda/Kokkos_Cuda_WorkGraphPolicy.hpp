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

#ifndef KOKKOS_CUDA_WORKGRAPHPOLICY_HPP
#define KOKKOS_CUDA_WORKGRAPHPOLICY_HPP

#include <Cuda/Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::WorkGraphPolicy<Traits...>,
                  Kokkos::Cuda> {
 public:
  using Policy = Kokkos::WorkGraphPolicy<Traits...>;
  using Self   = ParallelFor<FunctorType, Policy, Kokkos::Cuda>;

 private:
  Policy m_policy;
  FunctorType m_functor;

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_one(
      const std::int32_t w) const noexcept {
    m_functor(w);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_one(
      const std::int32_t w) const noexcept {
    const TagType t{};
    m_functor(t, w);
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  __device__ inline void operator()() const noexcept {
    // The following makes most threads idle,
    // which helps significantly with throughput due to reducing conflict rates
    // on the work acquisition, updated based on perf experiments of the
    // static Fibonacci experiment on Volta
    if (0 == (threadIdx.y % 4)) {
      // Spin until COMPLETED_TOKEN.
      // END_TOKEN indicates no work is currently available.

      for (std::int32_t w = Policy::END_TOKEN;
           Policy::COMPLETED_TOKEN != (w = m_policy.pop_work());) {
        if (Policy::END_TOKEN != w) {
          exec_one<typename Policy::work_tag>(w);
          m_policy.completed_work(w);
        }
// On pre-volta architectures we need a __syncwarp here to prevent
// infinite loops depending on the scheduling order above
#if defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL) || \
    defined(KOKKOS_ARCH_PASCAL)
        __syncwarp(__activemask());
#endif
      }
    }
  }

  inline void execute() {
    const int warps_per_block = 4;
    const int multi_processor_count =
        m_policy.space().cuda_device_prop().multiProcessorCount;
    const dim3 grid(multi_processor_count, 1, 1);
    const dim3 block(1, Kokkos::Impl::CudaTraits::WarpSize, warps_per_block);
    const int shared = 0;

    Kokkos::Impl::CudaParallelLaunch<Self>(
        *this, grid, block, shared, Cuda().impl_internal_space_instance());
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_policy(arg_policy), m_functor(arg_functor) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_CUDA_WORKGRAPHPOLICY_HPP */
