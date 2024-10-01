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

#ifndef KOKKOS_OPENMP_PARALLEL_REDUCE_HPP
#define KOKKOS_OPENMP_PARALLEL_REDUCE_HPP

#include <omp.h>
#include <OpenMP/Kokkos_OpenMP_Instance.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::OpenMP> {
 private:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  OpenMPInternal* m_instance;
  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update) {
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(iwork, update);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update) {
    const TagType t{};
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(t, iwork, update);
    }
  }

 public:
  inline void execute() const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    if (m_policy.end() <= m_policy.begin()) {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
      return;
    }
    enum {
      is_dynamic = std::is_same<typename Policy::schedule_type::type,
                                Kokkos::Dynamic>::value
    };

    const size_t pool_reduce_bytes = reducer.value_size();

    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);

    m_instance->resize_thread_data(pool_reduce_bytes, 0  // team_reduce_bytes
                                   ,
                                   0  // team_shared_bytes
                                   ,
                                   0  // thread_local_bytes
    );

    if (execute_in_serial(m_policy.space())) {
      const pointer_type ptr =
          m_result_ptr
              ? m_result_ptr
              : pointer_type(
                    m_instance->get_thread_data(0)->pool_reduce_local());

      reference_type update = reducer.init(ptr);

      ParallelReduce::template exec_range<WorkTag>(
          m_functor_reducer.get_functor(), m_policy.begin(), m_policy.end(),
          update);

      reducer.final(ptr);

      return;
    }
    const int pool_size = m_instance->thread_pool_size();
#pragma omp parallel num_threads(pool_size)
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());

      data.set_work_partition(m_policy.end() - m_policy.begin(),
                              m_policy.chunk_size());

      if (is_dynamic) {
        // Make sure work partition is set before stealing
        if (data.pool_rendezvous()) data.pool_rendezvous_release();
      }

      reference_type update = reducer.init(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()));

      std::pair<int64_t, int64_t> range(0, 0);

      do {
        range = is_dynamic ? data.get_work_stealing_chunk()
                           : data.get_work_partition();

        ParallelReduce::template exec_range<WorkTag>(
            m_functor_reducer.get_functor(), range.first + m_policy.begin(),
            range.second + m_policy.begin(), update);

      } while (is_dynamic && 0 <= range.first);
    }

    // Reduction:

    const pointer_type ptr =
        pointer_type(m_instance->get_thread_data(0)->pool_reduce_local());

    for (int i = 1; i < pool_size; ++i) {
      reducer.join(ptr,
                   reinterpret_cast<pointer_type>(
                       m_instance->get_thread_data(i)->pool_reduce_local()));
    }

    reducer.final(ptr);

    if (m_result_ptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = ptr[j];
      }
    }
  }

  //----------------------------------------

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                        Policy arg_policy, const ViewType& arg_view)
      : m_instance(nullptr),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_view.data()) {
    m_instance = arg_policy.space().impl_internal_space_instance();
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP reduce result must be a View accessible from "
        "HostSpace");
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// MDRangePolicy impl
template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>, Kokkos::OpenMP> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename MDRangePolicy::work_tag;
  using Member  = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, CombinedFunctorReducerType, WorkTag, reference_type>;

  OpenMPInternal* m_instance;
  const iterate_type m_iter;
  const pointer_type m_result_ptr;

  inline void exec_range(const Member ibeg, const Member iend,
                         reference_type update) const {
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      m_iter(iwork, update);
    }
  }

 public:
  inline void execute() const {
    const ReducerType& reducer     = m_iter.m_func.get_reducer();
    const size_t pool_reduce_bytes = reducer.value_size();

    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);

    m_instance->resize_thread_data(pool_reduce_bytes, 0  // team_reduce_bytes
                                   ,
                                   0  // team_shared_bytes
                                   ,
                                   0  // thread_local_bytes
    );

#ifndef KOKKOS_COMPILER_INTEL
    if (execute_in_serial(m_iter.m_rp.space())) {
      const pointer_type ptr =
          m_result_ptr
              ? m_result_ptr
              : pointer_type(
                    m_instance->get_thread_data(0)->pool_reduce_local());

      reference_type update = reducer.init(ptr);

      ParallelReduce::exec_range(0, m_iter.m_rp.m_num_tiles, update);

      reducer.final(ptr);

      return;
    }
#endif

    enum {
      is_dynamic = std::is_same<typename Policy::schedule_type::type,
                                Kokkos::Dynamic>::value
    };

    const int pool_size = m_instance->thread_pool_size();
#pragma omp parallel num_threads(pool_size)
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());

      data.set_work_partition(m_iter.m_rp.m_num_tiles, 1);

      if (is_dynamic) {
        // Make sure work partition is set before stealing
        if (data.pool_rendezvous()) data.pool_rendezvous_release();
      }

      reference_type update = reducer.init(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()));

      std::pair<int64_t, int64_t> range(0, 0);

      do {
        range = is_dynamic ? data.get_work_stealing_chunk()
                           : data.get_work_partition();

        ParallelReduce::exec_range(range.first, range.second, update);

      } while (is_dynamic && 0 <= range.first);
    }
    // END #pragma omp parallel

    // Reduction:

    const pointer_type ptr =
        pointer_type(m_instance->get_thread_data(0)->pool_reduce_local());

    for (int i = 1; i < pool_size; ++i) {
      reducer.join(ptr,
                   reinterpret_cast<pointer_type>(
                       m_instance->get_thread_data(i)->pool_reduce_local()));
    }

    reducer.final(ptr);

    if (m_result_ptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = ptr[j];
      }
    }
  }

  //----------------------------------------

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 MDRangePolicy arg_policy, const ViewType& arg_view)
      : m_instance(nullptr),
        m_iter(arg_policy, arg_functor_reducer),
        m_result_ptr(arg_view.data()) {
    m_instance = arg_policy.space().impl_internal_space_instance();
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP reduce result must be a View accessible from "
        "HostSpace");
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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, Kokkos::OpenMP> {
 private:
  enum { TEAM_REDUCE_SIZE = 512 };

  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::OpenMP, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag  = typename Policy::work_tag;
  using SchedTag = typename Policy::schedule_type::type;
  using Member   = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  OpenMPInternal* m_instance;
  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const int m_shmem_size;

  template <class TagType>
  inline static std::enable_if_t<(std::is_void<TagType>::value)> exec_team(
      const FunctorType& functor, HostThreadTeamData& data,
      reference_type& update, const int league_rank_begin,
      const int league_rank_end, const int league_size) {
    for (int r = league_rank_begin; r < league_rank_end;) {
      functor(Member(data, r, league_size), update);

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
      reference_type& update, const int league_rank_begin,
      const int league_rank_end, const int league_size) {
    const TagType t{};

    for (int r = league_rank_begin; r < league_rank_end;) {
      functor(t, Member(data, r, league_size), update);

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

    const ReducerType& reducer = m_functor_reducer.get_reducer();

    if (m_policy.league_size() == 0 || m_policy.team_size() == 0) {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
      return;
    }

    const size_t pool_reduce_size = reducer.value_size();

    const size_t team_reduce_size  = TEAM_REDUCE_SIZE * m_policy.team_size();
    const size_t team_shared_size  = m_shmem_size + m_policy.scratch_size(1);
    const size_t thread_local_size = 0;  // Never shrinks

    // Serialize kernels on the same execution space instance
    std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);

    m_instance->resize_thread_data(pool_reduce_size, team_reduce_size,
                                   team_shared_size, thread_local_size);

    if (execute_in_serial(m_policy.space())) {
      HostThreadTeamData& data = *(m_instance->get_thread_data());
      pointer_type ptr =
          m_result_ptr ? m_result_ptr : pointer_type(data.pool_reduce_local());
      reference_type update       = reducer.init(ptr);
      const int league_rank_begin = 0;
      const int league_rank_end   = m_policy.league_size();
      ParallelReduce::template exec_team<WorkTag>(
          m_functor_reducer.get_functor(), data, update, league_rank_begin,
          league_rank_end, m_policy.league_size());

      reducer.final(ptr);

      return;
    }

    const int pool_size = m_instance->thread_pool_size();
#pragma omp parallel num_threads(pool_size)
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
        reference_type update = reducer.init(
            reinterpret_cast<pointer_type>(data.pool_reduce_local()));

        std::pair<int64_t, int64_t> range(0, 0);

        do {
          range = is_dynamic ? data.get_work_stealing_chunk()
                             : data.get_work_partition();

          ParallelReduce::template exec_team<WorkTag>(
              m_functor_reducer.get_functor(), data, update, range.first,
              range.second, m_policy.league_size());

        } while (is_dynamic && 0 <= range.first);
      } else {
        reducer.init(reinterpret_cast<pointer_type>(data.pool_reduce_local()));
      }

      data.disband_team();

      //  This thread has updated 'pool_reduce_local()' with its
      //  contributions to the reduction.  The parallel region is
      //  about to terminate and the master thread will load and
      //  reduce each 'pool_reduce_local()' contribution.
      //  Must 'memory_fence()' to guarantee that storing the update to
      //  'pool_reduce_local()' will complete before this thread
      //  exits the parallel region.

      memory_fence();
    }

    // Reduction:

    const pointer_type ptr =
        pointer_type(m_instance->get_thread_data(0)->pool_reduce_local());

    for (int i = 1; i < pool_size; ++i) {
      reducer.join(ptr,
                   reinterpret_cast<pointer_type>(
                       m_instance->get_thread_data(i)->pool_reduce_local()));
    }

    reducer.final(ptr);

    if (m_result_ptr) {
      const int n = reducer.value_count();

      for (int j = 0; j < n; ++j) {
        m_result_ptr[j] = ptr[j];
      }
    }
  }

  //----------------------------------------

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                        const Policy& arg_policy, const ViewType& arg_result)
      : m_instance(nullptr),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_shmem_size(
            arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
            FunctorTeamShmemSize<FunctorType>::value(
                arg_functor_reducer.get_functor(), arg_policy.team_size())) {
    m_instance = arg_policy.space().impl_internal_space_instance();

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP reduce result must be a View accessible from "
        "HostSpace");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* KOKKOS_OPENMP_PARALLEL_REDUCE_HPP */
