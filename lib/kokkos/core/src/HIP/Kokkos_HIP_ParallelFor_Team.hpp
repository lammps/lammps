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

#ifndef KOKKOS_HIP_PARALLEL_FOR_TEAM_HPP
#define KOKKOS_HIP_PARALLEL_FOR_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_Team.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <HIP/Kokkos_HIP_TeamPolicyInternal.hpp>

namespace Kokkos {
namespace Impl {

template <typename FunctorType, typename... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>, HIP> {
 public:
  using Policy       = TeamPolicy<Properties...>;
  using functor_type = FunctorType;
  using size_type    = HIP::size_type;

 private:
  using member_type   = typename Policy::member_type;
  using work_tag      = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  // Algorithmic constraints: blockDim.y is a power of two AND
  // blockDim.y  == blockDim.z == 1 shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]

  FunctorType const m_functor;
  Policy const m_policy;
  size_type const m_league_size;
  int m_team_size;
  size_type const m_vector_size;
  int m_shmem_begin;
  int m_shmem_size;
  void* m_scratch_ptr[2];
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  size_t m_num_scratch_locks;

  template <typename TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const member_type& member) const {
    m_functor(member);
  }

  template <typename TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const member_type& member) const {
    m_functor(TagType(), member);
  }

 public:
  ParallelFor()                   = delete;
  ParallelFor(ParallelFor const&) = default;
  ParallelFor& operator=(ParallelFor const&) = delete;

  __device__ inline void operator()() const {
    // Iterate this block through the league
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = hip_get_scratch_index(m_league_size, m_scratch_locks,
                                       m_num_scratch_locks);
    }

    int const int_league_size = static_cast<int>(m_league_size);
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<work_tag>(typename Policy::member_type(
          kokkos_impl_hip_shared_memory<void>(), m_shmem_begin, m_shmem_size,
          static_cast<void*>(static_cast<char*>(m_scratch_ptr[1]) +
                             ptrdiff_t(threadid / (blockDim.x * blockDim.y)) *
                                 m_scratch_size[1]),
          m_scratch_size[1], league_rank, m_league_size));
    }
    if (m_scratch_size[1] > 0) {
      hip_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  inline void execute() const {
    int64_t const shmem_size_total = m_shmem_begin + m_shmem_size;
    dim3 const grid(static_cast<int>(m_league_size), 1, 1);
    dim3 const block(static_cast<int>(m_vector_size),
                     static_cast<int>(m_team_size), 1);

    using closure_type =
        ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>, HIP>;
    Impl::hip_parallel_launch<closure_type, launch_bounds>(
        *this, grid, block, shmem_size_total,
        m_policy.space().impl_internal_space_instance(),
        true);  // copy to device and execute
  }

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    auto internal_space_instance =
        m_policy.space().impl_internal_space_instance();
    m_team_size = m_team_size >= 0 ? m_team_size
                                   : arg_policy.team_size_recommended(
                                         arg_functor, ParallelForTag());

    m_shmem_begin = (sizeof(double) * (m_team_size + 2));
    m_shmem_size =
        (m_policy.scratch_size(0, m_team_size) +
         FunctorTeamShmemSize<FunctorType>::value(m_functor, m_team_size));
    m_scratch_size[0]   = m_policy.scratch_size(0, m_team_size);
    m_scratch_size[1]   = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks     = internal_space_instance->m_scratch_locks;
    m_num_scratch_locks = internal_space_instance->m_num_scratch_locks;

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    m_scratch_ptr[0] = nullptr;
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

    int const shmem_size_total = m_shmem_begin + m_shmem_size;
    if (internal_space_instance->m_maxShmemPerBlock < shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > insufficient shared memory"));
    }

    size_t max_size = arg_policy.team_size_max(arg_functor, ParallelForTag());
    if (static_cast<int>(m_team_size) > static_cast<int>(max_size)) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > requested too large team size."));
    }
  }

  ~ParallelFor() {
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
