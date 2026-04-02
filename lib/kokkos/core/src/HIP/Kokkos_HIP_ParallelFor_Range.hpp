// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_PARALLEL_FOR_RANGE_HPP
#define KOKKOS_HIP_PARALLEL_FOR_RANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::HIP> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member          = typename Policy::member_type;
  using WorkTag         = typename Policy::work_tag;
  using LaunchBounds    = typename Policy::launch_bounds;
  using StaticBatchSize = typename Policy::static_batch_size;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline __device__ std::enable_if_t<std::is_void_v<TagType>> exec_range(
      const Member i) const {
    m_functor(i);
  }

  template <class TagType>
  inline __device__ std::enable_if_t<!std::is_void_v<TagType>> exec_range(
      const Member i) const {
    m_functor(TagType(), i);
  }

 public:
  using functor_type = FunctorType;

  ParallelFor() = delete;

  inline __device__ void operator()() const {
    constexpr auto batch_size = Member(StaticBatchSize::batch_size);
    const auto work_stride    = Member(blockDim.y) * gridDim.x;
    const Member work_end     = m_policy.end();

    for (Member iwork = m_policy.begin() + threadIdx.y +
                        static_cast<Member>(blockDim.y) * blockIdx.x;
         iwork < work_end;
         iwork =
             iwork < static_cast<Member>(work_end - work_stride * batch_size)
                 ? iwork + work_stride * batch_size
                 : work_end) {
      for (Member i = 0; i < static_cast<Member>(work_stride * batch_size) &&
                         i < work_end - iwork;
           i = (i < static_cast<Member>(work_end - work_stride - iwork))
                   ? i + work_stride
                   : work_end - iwork) {
        this->template exec_range<WorkTag>(iwork + i);
      }
    }
  }

  inline void execute() const {
    constexpr typename Policy::index_type batch_size =
        StaticBatchSize::batch_size;
    const typename Policy::index_type nwork =
        (m_policy.end() - m_policy.begin()) / batch_size +
        ((m_policy.end() - m_policy.begin()) % batch_size == 0 ? 0 : 1);

    using DriverType = ParallelFor<FunctorType, Policy, Kokkos::HIP>;
    const int block_size =
        Kokkos::Impl::get_preferred_blocksize_for_range<DriverType,
                                                        LaunchBounds>(
            m_policy.space().impl_internal_space_instance(), nwork);

    if (block_size == 0) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelFor< HIP > could not find a "
                      "valid execution configuration."));
    }
    const dim3 block(1, block_size, 1);
    const int maxGridSizeX = m_policy.space().hip_device_prop().maxGridSize[0];
    const dim3 grid(
        std::min(typename Policy::index_type((nwork + block.y - 1) / block.y),
                 typename Policy::index_type(maxGridSizeX)),
        1, 1);

    Kokkos::Impl::hip_parallel_launch<DriverType, LaunchBounds>(
        *this, grid, block, 0, m_policy.space().impl_internal_space_instance(),
        false);
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
