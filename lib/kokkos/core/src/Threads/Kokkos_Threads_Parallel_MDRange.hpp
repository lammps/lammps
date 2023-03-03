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

#ifndef KOKKOS_THREADS_PARALLEL_MDRANGE_HPP
#define KOKKOS_THREADS_PARALLEL_MDRANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Threads> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;

  using WorkTag = typename MDRangePolicy::work_tag;

  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

  const iterate_type m_iter;

  inline void exec_range(const Member ibeg, const Member iend) const {
    for (Member i = ibeg; i < iend; ++i) {
      m_iter(i);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    auto const num_tiles = self.m_iter.m_rp.m_num_tiles;
    WorkRange range(Policy(0, num_tiles).set_chunk_size(1), exec.pool_rank(),
                    exec.pool_size());

    self.exec_range(range.begin(), range.end());

    exec.fan_in();
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    auto const num_tiles = self.m_iter.m_rp.m_num_tiles;
    WorkRange range(Policy(0, num_tiles).set_chunk_size(1), exec.pool_rank(),
                    exec.pool_size());

    exec.set_work_range(range.begin(), range.end(), 1);
    exec.reset_steal_target();
    exec.barrier();

    long work_index = exec.get_work_index();

    while (work_index != -1) {
      const Member begin = static_cast<Member>(work_index);
      const Member end   = begin + 1 < num_tiles ? begin + 1 : num_tiles;

      self.exec_range(begin, end);
      work_index = exec.get_work_index();
    }

    exec.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsExec::start(&ParallelFor::exec, this);
    ThreadsExec::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const MDRangePolicy &arg_policy)
      : m_iter(arg_policy, arg_functor) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy &, const Functor &) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Threads> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;

  using WorkTag   = typename MDRangePolicy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                         MDRangePolicy, ReducerTypeFwd>;
  using pointer_type   = typename Analysis::pointer_type;
  using value_type     = typename Analysis::value_type;
  using reference_type = typename Analysis::reference_type;

  using iterate_type =
      typename Kokkos::Impl::HostIterateTile<MDRangePolicy, FunctorType,
                                             WorkTag, reference_type>;

  const iterate_type m_iter;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;

  inline void exec_range(const Member &ibeg, const Member &iend,
                         reference_type update) const {
    for (Member i = ibeg; i < iend; ++i) {
      m_iter(i, update);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);

    const auto num_tiles = self.m_iter.m_rp.m_num_tiles;
    const WorkRange range(Policy(0, num_tiles).set_chunk_size(1),
                          exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer reducer(
        &ReducerConditional::select(self.m_iter.m_func, self.m_reducer));

    self.exec_range(
        range.begin(), range.end(),
        reducer.init(static_cast<pointer_type>(exec.reduce_memory())));

    exec.fan_in_reduce(reducer);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);

    const auto num_tiles = self.m_iter.m_rp.m_num_tiles;
    const WorkRange range(Policy(0, num_tiles).set_chunk_size(1),
                          exec.pool_rank(), exec.pool_size());

    exec.set_work_range(range.begin(), range.end(), 1);
    exec.reset_steal_target();
    exec.barrier();

    long work_index = exec.get_work_index();
    typename Analysis::Reducer reducer(
        &ReducerConditional::select(self.m_iter.m_func, self.m_reducer));

    reference_type update =
        reducer.init(static_cast<pointer_type>(exec.reduce_memory()));
    while (work_index != -1) {
      const Member begin = static_cast<Member>(work_index);
      const Member end   = begin + 1 < num_tiles ? begin + 1 : num_tiles;
      self.exec_range(begin, end, update);
      work_index = exec.get_work_index();
    }

    exec.fan_in_reduce(reducer);
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(Analysis::value_size(ReducerConditional::select(
                                    m_iter.m_func, m_reducer)),
                                0);

    ThreadsExec::start(&ParallelReduce::exec, this);

    ThreadsExec::fence();

    if (m_result_ptr) {
      const pointer_type data =
          (pointer_type)ThreadsExec::root_reduce_scratch();

      const unsigned n = Analysis::value_count(
          ReducerConditional::select(m_iter.m_func, m_reducer));
      for (unsigned i = 0; i < n; ++i) {
        m_result_ptr[i] = data[i];
      }
    }
  }

  template <class HostViewType>
  ParallelReduce(const FunctorType &arg_functor,
                 const MDRangePolicy &arg_policy,
                 const HostViewType &arg_result_view,
                 std::enable_if_t<Kokkos::is_view<HostViewType>::value &&
                                      !Kokkos::is_reducer<ReducerType>::value,
                                  void *> = nullptr)
      : m_iter(arg_policy, arg_functor),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<HostViewType>::value,
                  "Kokkos::Threads reduce result must be a View");

    static_assert(
        std::is_same<typename HostViewType::memory_space, HostSpace>::value,
        "Kokkos::Threads reduce result must be a View in HostSpace");
  }

  inline ParallelReduce(const FunctorType &arg_functor,
                        MDRangePolicy arg_policy, const ReducerType &reducer)
      : m_iter(arg_policy, arg_functor),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {
    /*static_assert( std::is_same< typename ViewType::memory_space
                                    , Kokkos::HostSpace >::value
      , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace"
      );*/
  }

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy &, const Functor &) {
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

#endif
