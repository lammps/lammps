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

#ifndef KOKKOS_OPENMP_PARALLEL_HPP
#define KOKKOS_OPENMP_PARALLEL_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_OPENMP)

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

inline bool execute_in_serial(OpenMP const& space = OpenMP()) {
  return (OpenMP::in_parallel(space) &&
          !(omp_get_nested() && (omp_get_level() == 1)));
}

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
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

namespace Kokkos {
namespace Impl {

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

    m_instance->acquire_lock();

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

    m_instance->release_lock();
  }

  //----------------------------------------

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                        Policy arg_policy, const ViewType& arg_view)
      : m_instance(nullptr),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_view.data()) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP reduce result must be a View accessible from "
        "HostSpace");
  }
};

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

    m_instance->acquire_lock();

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

      m_instance->release_lock();

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

    m_instance->release_lock();
  }

  //----------------------------------------

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 MDRangePolicy arg_policy, const ViewType& arg_view)
      : m_instance(nullptr),
        m_iter(arg_policy, arg_functor_reducer),
        m_result_ptr(arg_view.data()) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
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

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::OpenMP> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::SCAN, Policy, FunctorType, void>;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  OpenMPInternal* m_instance;
  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update, const bool final) {
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(iwork, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update, const bool final) {
    const TagType t{};
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(t, iwork, update, final);
    }
  }

 public:
  inline void execute() const {
    const int value_count          = Analysis::value_count(m_functor);
    const size_t pool_reduce_bytes = 2 * Analysis::value_size(m_functor);

    m_instance->resize_thread_data(pool_reduce_bytes, 0  // team_reduce_bytes
                                   ,
                                   0  // team_shared_bytes
                                   ,
                                   0  // thread_local_bytes
    );

    if (execute_in_serial(m_policy.space())) {
      typename Analysis::Reducer final_reducer(m_functor);

      reference_type update = final_reducer.init(
          pointer_type(m_instance->get_thread_data(0)->pool_reduce_local()));

      ParallelScan::template exec_range<WorkTag>(m_functor, m_policy.begin(),
                                                 m_policy.end(), update, true);

      return;
    }

#pragma omp parallel num_threads(m_instance->thread_pool_size())
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());
      typename Analysis::Reducer final_reducer(m_functor);

      const WorkRange range(m_policy, omp_get_thread_num(),
                            omp_get_num_threads());

      reference_type update_sum = final_reducer.init(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()));

      ParallelScan::template exec_range<WorkTag>(
          m_functor, range.begin(), range.end(), update_sum, false);

      if (data.pool_rendezvous()) {
        pointer_type ptr_prev = nullptr;

        const int n = omp_get_num_threads();

        for (int i = 0; i < n; ++i) {
          pointer_type ptr =
              (pointer_type)data.pool_member(i)->pool_reduce_local();

          if (i) {
            for (int j = 0; j < value_count; ++j) {
              ptr[j + value_count] = ptr_prev[j + value_count];
            }
            final_reducer.join(ptr + value_count, ptr_prev);
          } else {
            final_reducer.init(ptr + value_count);
          }

          ptr_prev = ptr;
        }

        data.pool_rendezvous_release();
      }

      reference_type update_base = final_reducer.reference(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()) +
          value_count);

      ParallelScan::template exec_range<WorkTag>(
          m_functor, range.begin(), range.end(), update_base, true);
    }
  }

  //----------------------------------------

  inline ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_instance(nullptr), m_functor(arg_functor), m_policy(arg_policy) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
  }

  //----------------------------------------
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::OpenMP> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using Analysis = FunctorAnalysis<FunctorPatternInterface::SCAN, Policy,
                                   FunctorType, ReturnType>;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using value_type     = typename Analysis::value_type;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  OpenMPInternal* m_instance;
  const FunctorType m_functor;
  const Policy m_policy;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update, const bool final) {
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(iwork, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType& functor, const Member ibeg, const Member iend,
      reference_type update, const bool final) {
    const TagType t{};
    for (Member iwork = ibeg; iwork < iend; ++iwork) {
      functor(t, iwork, update, final);
    }
  }

 public:
  inline void execute() const {
    const int value_count          = Analysis::value_count(m_functor);
    const size_t pool_reduce_bytes = 2 * Analysis::value_size(m_functor);

    m_instance->acquire_lock();

    m_instance->resize_thread_data(pool_reduce_bytes, 0  // team_reduce_bytes
                                   ,
                                   0  // team_shared_bytes
                                   ,
                                   0  // thread_local_bytes
    );

    if (execute_in_serial(m_policy.space())) {
      typename Analysis::Reducer final_reducer(m_functor);

      reference_type update = final_reducer.init(
          pointer_type(m_instance->get_thread_data(0)->pool_reduce_local()));

      this->template exec_range<WorkTag>(m_functor, m_policy.begin(),
                                         m_policy.end(), update, true);

      *m_result_ptr = update;

      m_instance->release_lock();

      return;
    }

#pragma omp parallel num_threads(m_instance->thread_pool_size())
    {
      HostThreadTeamData& data = *(m_instance->get_thread_data());
      typename Analysis::Reducer final_reducer(m_functor);

      const WorkRange range(m_policy, omp_get_thread_num(),
                            omp_get_num_threads());
      reference_type update_sum = final_reducer.init(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()));

      ParallelScanWithTotal::template exec_range<WorkTag>(
          m_functor, range.begin(), range.end(), update_sum, false);

      if (data.pool_rendezvous()) {
        pointer_type ptr_prev = nullptr;

        const int n = omp_get_num_threads();

        for (int i = 0; i < n; ++i) {
          pointer_type ptr =
              (pointer_type)data.pool_member(i)->pool_reduce_local();

          if (i) {
            for (int j = 0; j < value_count; ++j) {
              ptr[j + value_count] = ptr_prev[j + value_count];
            }
            final_reducer.join(ptr + value_count, ptr_prev);
          } else {
            final_reducer.init(ptr + value_count);
          }

          ptr_prev = ptr;
        }

        data.pool_rendezvous_release();
      }

      reference_type update_base = final_reducer.reference(
          reinterpret_cast<pointer_type>(data.pool_reduce_local()) +
          value_count);

      ParallelScanWithTotal::template exec_range<WorkTag>(
          m_functor, range.begin(), range.end(), update_base, true);

      if (omp_get_thread_num() == omp_get_num_threads() - 1) {
        *m_result_ptr = update_base;
      }
    }

    m_instance->release_lock();
  }

  //----------------------------------------

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const Policy& arg_policy,
                        const ViewType& arg_result_view)
      : m_instance(nullptr),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP parallel_scan result must be host-accessible!");
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
  }

  //----------------------------------------
};

}  // namespace Impl
}  // namespace Kokkos

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

    m_instance->acquire_lock();

    m_instance->resize_thread_data(pool_reduce_size, team_reduce_size,
                                   team_shared_size, thread_local_size);

    if (execute_in_serial(m_policy.space())) {
      ParallelFor::template exec_team<WorkTag>(
          m_functor, *(m_instance->get_thread_data()), 0,
          m_policy.league_size(), m_policy.league_size());

      m_instance->release_lock();

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

    m_instance->release_lock();
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_instance(nullptr),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif
  }
};

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

    m_instance->acquire_lock();

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

      m_instance->release_lock();

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

    m_instance->release_lock();
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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
    if (t_openmp_instance) {
      m_instance = t_openmp_instance;
    } else {
      m_instance = arg_policy.space().impl_internal_space_instance();
    }
#else
    m_instance = arg_policy.space().impl_internal_space_instance();
#endif

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::OpenMP reduce result must be a View accessible from "
        "HostSpace");
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef KOKKOS_PRAGMA_IVDEP_IF_ENABLED
#undef KOKKOS_OPENMP_OPTIONAL_CHUNK_SIZE

#endif
#endif /* KOKKOS_OPENMP_PARALLEL_HPP */
