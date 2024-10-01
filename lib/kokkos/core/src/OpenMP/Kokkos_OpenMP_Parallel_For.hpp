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

#ifndef KOKKOS_OPENMP_PARALLEL_FOR_HPP
#define KOKKOS_OPENMP_PARALLEL_FOR_HPP

#include <omp.h>
#include <OpenMP/Kokkos_OpenMP_Instance.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOS_PRAGMA_IVDEP_IF_ENABLED
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#undef KOKKOS_PRAGMA_IVDEP_IF_ENABLED
#define KOKKOS_PRAGMA_IVDEP_IF_ENABLED _Pragma("ivdep")
#endif

#ifndef KOKKOS_COMPILER_NVHPC
#define KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE , m_policy.chunk_size()
#else
#define KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE
#endif

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::OpenMP> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  OpenMPInternal* m_instance;
  const FunctorType m_functor;
  const Policy m_policy;

  inline static void exec_range(const FunctorType& functor, const Member ibeg,
                                const Member iend) {
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (auto iwork = ibeg; iwork < iend; ++iwork) {
      exec_work(functor, iwork);
    }
  }

  template <class Enable = WorkTag>
  inline static std::enable_if_t<std::is_void<WorkTag>::value &&
                                 std::is_same<Enable, WorkTag>::value>
  exec_work(const FunctorType& functor, const Member iwork) {
    functor(iwork);
  }

  template <class Enable = WorkTag>
  inline static std::enable_if_t<!std::is_void<WorkTag>::value &&
                                 std::is_same<Enable, WorkTag>::value>
  exec_work(const FunctorType& functor, const Member iwork) {
    functor(WorkTag{}, iwork);
  }

  template <class Policy>
  std::enable_if_t<std::is_same<typename Policy::schedule_type::type,
                                Kokkos::Dynamic>::value>
  execute_parallel() const {
    // prevent bug in NVHPC 21.9/CUDA 11.4 (entering zero iterations loop)
    if (m_policy.begin() >= m_policy.end()) return;
#pragma omp parallel for schedule(dynamic KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE) \
    num_threads(m_instance->thread_pool_size())
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (auto iwork = m_policy.begin(); iwork < m_policy.end(); ++iwork) {
      exec_work(m_functor, iwork);
    }
  }

  template <class Policy>
  std::enable_if_t<!std::is_same<typename Policy::schedule_type::type,
                                 Kokkos::Dynamic>::value>
  execute_parallel() const {
// Specifying an chunksize with GCC compiler leads to performance regression
// with static schedule.
#ifdef KOKKOS_COMPILER_GNU
#pragma omp parallel for schedule(static) \
    num_threads(m_instance->thread_pool_size())
#else
#pragma omp parallel for schedule(static KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE) \
    num_threads(m_instance->thread_pool_size())
#endif
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (auto iwork = m_policy.begin(); iwork < m_policy.end(); ++iwork) {
      exec_work(m_functor, iwork);
    }
  }

 public:
  inline void execute() const {
    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);
    if (execute_in_serial(m_policy.space())) {
      exec_range(m_functor, m_policy.begin(), m_policy.end());
      return;
    }

#ifndef KOKKOS_INTERNAL_DISABLE_NATIVE_OPENMP
    execute_parallel<Policy>();
#else
    constexpr bool is_dynamic =
        std::is_same<typename Policy::schedule_type::type,
                     Kokkos::Dynamic>::value;
#pragma omp parallel num_threads(m_instance->thread_pool_size())
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());

      data.set_work_partition(m_policy.end() - m_policy.begin(),
                              m_policy.chunk_size());

      if (is_dynamic) {
        // Make sure work partition is set before stealing
        if (data.pool_rendezvous()) data.pool_rendezvous_release();
      }

      std::pair<int64_t, int64_t> range(0, 0);

      do {
        range = is_dynamic ? data.get_work_stealing_chunk()
                           : data.get_work_partition();

        exec_range(m_functor, range.first + m_policy.begin(),
                   range.second + m_policy.begin());

      } while (is_dynamic && 0 <= range.first);
    }
#endif
  }

  inline ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_instance(nullptr), m_functor(arg_functor), m_policy(arg_policy) {
    m_instance = arg_policy.space().impl_internal_space_instance();
  }
};

// MDRangePolicy impl
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::OpenMP> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using WorkTag       = typename MDRangePolicy::work_tag;

  using Member = typename Policy::member_type;

  using index_type   = typename Policy::index_type;
  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

  OpenMPInternal* m_instance;
  const iterate_type m_iter;

  inline void exec_range(const Member ibeg, const Member iend) const {
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      m_iter(iwork);
    }
  }

  template <class Policy>
  typename std::enable_if_t<std::is_same<typename Policy::schedule_type::type,
                                         Kokkos::Dynamic>::value>
  execute_parallel() const {
#pragma omp parallel for schedule(dynamic, 1) \
    num_threads(m_instance->thread_pool_size())
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (index_type iwork = 0; iwork < m_iter.m_rp.m_num_tiles; ++iwork) {
      m_iter(iwork);
    }
  }

  template <class Policy>
  typename std::enable_if<!std::is_same<typename Policy::schedule_type::type,
                                        Kokkos::Dynamic>::value>::type
  execute_parallel() const {
#pragma omp parallel for schedule(static, 1) \
    num_threads(m_instance->thread_pool_size())
    KOKKOS_PRAGMA_IVDEP_IF_ENABLED
    for (index_type iwork = 0; iwork < m_iter.m_rp.m_num_tiles; ++iwork) {
      m_iter(iwork);
    }
  }

 public:
  inline void execute() const {
    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);

#ifndef KOKKOS_COMPILER_INTEL
    if (execute_in_serial(m_iter.m_rp.space())) {
      exec_range(0, m_iter.m_rp.m_num_tiles);
      return;
    }
#endif

#ifndef KOKKOS_INTERNAL_DISABLE_NATIVE_OPENMP
    execute_parallel<Policy>();
#else
    constexpr bool is_dynamic =
        std::is_same<typename Policy::schedule_type::type,
                     Kokkos::Dynamic>::value;

#pragma omp parallel num_threads(m_instance->thread_pool_size())
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());

      data.set_work_partition(m_iter.m_rp.m_num_tiles, 1);

      if (is_dynamic) {
        // Make sure work partition is set before stealing
        if (data.pool_rendezvous()) data.pool_rendezvous_release();
      }

      std::pair<int64_t, int64_t> range(0, 0);

      do {
        range = is_dynamic ? data.get_work_stealing_chunk()
                           : data.get_work_partition();

        exec_range(range.first, range.second);

      } while (is_dynamic && 0 <= range.first);
    }
    // END #pragma omp parallel
#endif
  }

  inline ParallelFor(const FunctorType& arg_functor, MDRangePolicy arg_policy)
      : m_instance(nullptr), m_iter(arg_policy, arg_functor) {
    m_instance = arg_policy.space().impl_internal_space_instance();
  }

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::OpenMP> {
 private:
  enum { TEAM_REDUCE_SIZE = 512 };

  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::OpenMP, Properties...>;
  using WorkTag  = typename Policy::work_tag;
  using SchedTag = typename Policy::schedule_type::type;
  using Member   = typename Policy::member_type;

  OpenMPInternal* m_instance;
  const FunctorType m_functor;
  const Policy m_policy;
  const size_t m_shmem_size;

  template <class TagType>
  inline static std::enable_if_t<(std::is_void<TagType>::value)> exec_team(
      const FunctorType& functor, HostThreadTeamData& data,
      const int league_rank_begin, const int league_rank_end,
      const int league_size) {
    for (int r = league_rank_begin; r < league_rank_end;) {
      functor(Member(data, r, league_size));

      if (++r < league_rank_end) {
        // Don't allow team members to lap one another
        // so that they don't overwrite shared memory.
        if (data.team_rendezvous()) {
          data.team_rendezvous_release();
        }
      }
    }
  }

  template <class TagType>
  inline static std::enable_if_t<(!std::is_void<TagType>::value)> exec_team(
      const FunctorType& functor, HostThreadTeamData& data,
      const int league_rank_begin, const int league_rank_end,
      const int league_size) {
    const TagType t{};

    for (int r = league_rank_begin; r < league_rank_end;) {
      functor(t, Member(data, r, league_size));

      if (++r < league_rank_end) {
        // Don't allow team members to lap one another
        // so that they don't overwrite shared memory.
        if (data.team_rendezvous()) {
          data.team_rendezvous_release();
        }
      }
    }
  }

 public:
  inline void execute() const {
    enum { is_dynamic = std::is_same<SchedTag, Kokkos::Dynamic>::value };

    const size_t pool_reduce_size  = 0;  // Never shrinks
    const size_t team_reduce_size  = TEAM_REDUCE_SIZE * m_policy.team_size();
    const size_t team_shared_size  = m_shmem_size;
    const size_t thread_local_size = 0;  // Never shrinks

    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);

    m_instance->resize_thread_data(pool_reduce_size, team_reduce_size,
                                   team_shared_size, thread_local_size);

    if (execute_in_serial(m_policy.space())) {
      ParallelFor::template exec_team<WorkTag>(
          m_functor, *(m_instance->get_thread_data()), 0,
          m_policy.league_size(), m_policy.league_size());

      return;
    }

#pragma omp parallel num_threads(m_instance->thread_pool_size())
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());

      const int active = data.organize_team(m_policy.team_size());

      if (active) {
        data.set_work_partition(
            m_policy.league_size(),
            (0 < m_policy.chunk_size() ? m_policy.chunk_size()
                                       : m_policy.team_iter()));
      }

      if (is_dynamic) {
        // Must synchronize to make sure each team has set its
        // partition before beginning the work stealing loop.
        if (data.pool_rendezvous()) data.pool_rendezvous_release();
      }

      if (active) {
        std::pair<int64_t, int64_t> range(0, 0);

        do {
          range = is_dynamic ? data.get_work_stealing_chunk()
                             : data.get_work_partition();

          ParallelFor::template exec_team<WorkTag>(m_functor, data, range.first,
                                                   range.second,
                                                   m_policy.league_size());

        } while (is_dynamic && 0 <= range.first);
      }

      data.disband_team();
    }
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_instance(nullptr),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {
    m_instance = arg_policy.space().impl_internal_space_instance();
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef KOKKOS_PRAGMA_IVDEP_IF_ENABLED
#undef KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE

#endif /* KOKKOS_OPENMP_PARALLEL_FOR_HPP */
