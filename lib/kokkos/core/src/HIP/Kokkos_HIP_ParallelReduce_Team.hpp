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

#ifndef KOKKOS_HIP_PARALLEL_REDUCE_TEAM_HPP
#define KOKKOS_HIP_PARALLEL_REDUCE_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_Team.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP_TeamPolicyInternal.hpp>

namespace Kokkos {
namespace Impl {

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, HIP> {
 public:
  using Policy      = TeamPolicyInternal<HIP, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using member_type   = typename Policy::member_type;
  using work_tag      = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

 public:
  using functor_type = FunctorType;
  using size_type    = HIP::size_type;

  // static int constexpr UseShflReduction = false;
  // FIXME_HIP This should be disabled unconditionally for best performance, but
  // it currently causes tests to fail.
  static constexpr int UseShflReduction =
      (ReducerType::static_value_size() != 0);

 private:
  struct ShflReductionTag {};
  struct SHMEMReductionTag {};

  // Algorithmic constraints: blockDim.y is a power of two AND
  // blockDim.y == blockDim.z == 1 shared memory utilization:
  //
  //  [ global reduce space ]
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type m_team_begin;
  size_type m_shmem_begin;
  size_type m_shmem_size;
  void* m_scratch_ptr[2];
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  size_t m_num_scratch_locks;
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      member_type const& member, reference_type update) const {
    m_functor_reducer.get_functor()(member, update);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      member_type const& member, reference_type update) const {
    m_functor_reducer.get_functor()(TagType(), member, update);
  }

  __device__ inline void iterate_through_league(int const threadid,
                                                reference_type value) const {
    int const int_league_size = static_cast<int>(m_league_size);
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<work_tag>(
          member_type(
              kokkos_impl_hip_shared_memory<char>() + m_team_begin,
              m_shmem_begin, m_shmem_size,
              reinterpret_cast<void*>(
                  reinterpret_cast<char*>(m_scratch_ptr[1]) +
                  static_cast<ptrdiff_t>(threadid / (blockDim.x * blockDim.y)) *
                      m_scratch_size[1]),
              m_scratch_size[1], league_rank, m_league_size),
          value);
    }
  }

  int compute_block_count() const {
    constexpr auto light_weight =
        Kokkos::Experimental::WorkItemProperty::HintLightWeight;
    constexpr typename Policy::work_item_property property;
    // Numbers were tuned on MI210 using dot product and yAx benchmarks
    constexpr int block_max =
        (property & light_weight) == light_weight ? 2097152 : 65536;
    constexpr int preferred_block_min = 1024;
    int block_count                   = m_league_size;
    if (block_count < preferred_block_min) {
      // keep blocks as is, already low parallelism
    } else if (block_count >= block_max) {
      block_count = block_max;

    } else {
      int nwork = m_league_size * m_team_size;
      int items_per_thread =
          (nwork + block_count * m_team_size - 1) / (block_count * m_team_size);
      if (items_per_thread < 4) {
        int ratio = std::min(
            (block_count + preferred_block_min - 1) / preferred_block_min,
            (4 + items_per_thread - 1) / items_per_thread);
        block_count /= ratio;
      }
    }

    return block_count;
  }

 public:
  __device__ inline void operator()() const {
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = hip_get_scratch_index(m_league_size, m_scratch_locks,
                                       m_num_scratch_locks);
    }

    using ReductionTag = std::conditional_t<UseShflReduction, ShflReductionTag,
                                            SHMEMReductionTag>;
    run(ReductionTag{}, threadid);

    if (m_scratch_size[1] > 0) {
      hip_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  __device__ inline void run(SHMEMReductionTag, int const threadid) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    integral_nonzero_constant<size_type, ReducerType::static_value_size() /
                                             sizeof(size_type)> const
        word_count(reducer.value_size() / sizeof(size_type));

    reference_type value =
        reducer.init(kokkos_impl_hip_shared_memory<size_type>() +
                     threadIdx.y * word_count.value);
    // Iterate this block through the league
    iterate_through_league(threadid, value);

    // Reduce with final value at blockDim.y - 1 location.
    bool do_final_reduce = (m_league_size == 0);
    if (!do_final_reduce)
      do_final_reduce =
          hip_single_inter_block_reduce_scan<false, FunctorType, work_tag>(
              reducer, blockIdx.x, gridDim.x,
              kokkos_impl_hip_shared_memory<size_type>(), m_scratch_space,
              m_scratch_flags);
    if (do_final_reduce) {
      // This is the final block with the final result at the final threads'
      // location

      size_type* const shared = kokkos_impl_hip_shared_memory<size_type>() +
                                (blockDim.y - 1) * word_count.value;
      size_type* const global = m_result_ptr_device_accessible
                                    ? reinterpret_cast<size_type*>(m_result_ptr)
                                    : m_scratch_space;

      if (threadIdx.y == 0) {
        reducer.final(reinterpret_cast<value_type*>(shared));
      }

      if (HIPTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  __device__ inline void run(ShflReductionTag, int const threadid) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    value_type value;
    reducer.init(&value);

    // Iterate this block through the league
    iterate_through_league(threadid, value);

    pointer_type const result =
        m_result_ptr_device_accessible
            ? m_result_ptr
            : reinterpret_cast<pointer_type>(m_scratch_space);

    value_type init;
    reducer.init(&init);
    if (m_league_size == 0) {
      reducer.final(&value);
      *result = value;
    } else if (Impl::hip_inter_block_shuffle_reduction(
                   value, init, reducer, m_scratch_space, result,
                   m_scratch_flags, blockDim.y)) {
      unsigned int const id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        reducer.final(&value);
        *result = value;
      }
    }
  }

  inline void execute() {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    const bool is_empty_range  = m_league_size == 0 || m_team_size == 0;
    const bool need_device_set = ReducerType::has_init_member_function() ||
                                 ReducerType::has_final_member_function() ||
                                 !m_result_ptr_host_accessible ||
                                 Policy::is_graph_kernel::value ||
                                 !std::is_same<ReducerType, InvalidType>::value;
    if (!is_empty_range || need_device_set) {
      int const block_count = compute_block_count();

      m_scratch_space = hip_internal_scratch_space(
          m_policy.space(), reducer.value_size() * block_count);
      m_scratch_flags =
          hip_internal_scratch_flags(m_policy.space(), sizeof(size_type));

      dim3 block(m_vector_size, m_team_size, 1);
      dim3 grid(block_count, 1, 1);
      if (is_empty_range) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }
      const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

      Impl::hip_parallel_launch<ParallelReduce, launch_bounds>(
          *this, grid, block, shmem_size_total,
          m_policy.space().impl_internal_space_instance(),
          true);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().impl_internal_space_instance()->fence();

        if (m_result_ptr) {
          const int size = reducer.value_size();
          DeepCopy<HostSpace, HIPSpace, HIP>(m_policy.space(), m_result_ptr,
                                             m_scratch_space, size);
        }
      }
    } else {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(CombinedFunctorReducerType const& arg_functor_reducer,
                 Policy const& arg_policy, ViewType const& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<HIPSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_team_begin(0),
        m_shmem_begin(0),
        m_shmem_size(0),
        m_scratch_ptr{nullptr, nullptr},
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    auto internal_space_instance =
        m_policy.space().impl_internal_space_instance();
    m_team_size = m_team_size >= 0 ? m_team_size
                                   : arg_policy.team_size_recommended(
                                         arg_functor_reducer.get_functor(),
                                         arg_functor_reducer.get_reducer(),
                                         ParallelReduceTag());

    m_team_begin =
        UseShflReduction
            ? 0
            : hip_single_inter_block_reduce_scan_shmem<false, work_tag,
                                                       value_type>(
                  arg_functor_reducer.get_functor(), m_team_size);
    m_shmem_begin = sizeof(double) * (m_team_size + 2);
    m_shmem_size  = m_policy.scratch_size(0, m_team_size) +
                   FunctorTeamShmemSize<FunctorType>::value(
                       arg_functor_reducer.get_functor(), m_team_size);
    m_scratch_size[0]   = m_shmem_size;
    m_scratch_size[1]   = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks     = internal_space_instance->m_scratch_locks;
    m_num_scratch_locks = internal_space_instance->m_num_scratch_locks;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      m_scratch_pool_id = internal_space_instance->acquire_team_scratch_space();
      m_scratch_ptr[1]  = internal_space_instance->resize_team_scratch_space(
          m_scratch_pool_id,
          static_cast<std::int64_t>(m_scratch_size[1]) *
              (std::min(
                  static_cast<std::int64_t>(HIP().concurrency() /
                                            (m_team_size * m_vector_size)),
                  static_cast<std::int64_t>(m_league_size))));
    }

    // The global parallel_reduce does not support vector_length other than 1 at
    // the moment
    if ((arg_policy.impl_vector_length() > 1) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a vector length of "
          "greater than 1 is not currently supported for HIP for dynamic "
          "sized reduction types.");

    if ((m_team_size < HIPTraits::WarpSize) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller "
          "than 64 is not currently supported with HIP for dynamic sized "
          "reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

    if (!Kokkos::Impl::is_integral_power_of_two(m_team_size) &&
        !UseShflReduction) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > bad team size"));
    }

    if (internal_space_instance->m_maxShmemPerBlock < shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > requested too much "
                      "L0 scratch memory"));
    }

    size_t max_size = arg_policy.team_size_max(
        arg_functor_reducer.get_functor(), arg_functor_reducer.get_reducer(),
        ParallelReduceTag());
    if (static_cast<int>(m_team_size) > static_cast<int>(max_size)) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > requested too "
                      "large team size."));
    }
  }

  ~ParallelReduce() {
    if (m_scratch_pool_id >= 0) {
      m_policy.space()
          .impl_internal_space_instance()
          ->release_team_scratch_space(m_scratch_pool_id);
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
