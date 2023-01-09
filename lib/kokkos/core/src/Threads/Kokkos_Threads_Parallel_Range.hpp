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

#ifndef KOKKOS_THREADS_PARALLEL_RANGE_HPP
#define KOKKOS_THREADS_PARALLEL_RANGE_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member ibeg, const Member iend) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member ibeg, const Member iend) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    ParallelFor::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                              range.end());

    exec.fan_in();
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    exec.set_work_range(range.begin() - self.m_policy.begin(),
                        range.end() - self.m_policy.begin(),
                        self.m_policy.chunk_size());
    exec.reset_steal_target();
    exec.barrier();

    long work_index = exec.get_work_index();

    while (work_index != -1) {
      const Member begin =
          static_cast<Member>(work_index) * self.m_policy.chunk_size() +
          self.m_policy.begin();
      const Member end =
          begin + self.m_policy.chunk_size() < self.m_policy.end()
              ? begin + self.m_policy.chunk_size()
              : self.m_policy.end();
      ParallelFor::template exec_range<WorkTag>(self.m_functor, begin, end);
      work_index = exec.get_work_index();
    }

    exec.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsExec::start(&ParallelFor::exec, this);
    ThreadsExec::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Threads> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using WorkTag   = typename Policy::work_tag;
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
                                         Policy, ReducerTypeFwd>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer reducer(
        &ReducerConditional::select(self.m_functor, self.m_reducer));

    ParallelReduce::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(),
        reducer.init(static_cast<pointer_type>(exec.reduce_memory())));

    exec.fan_in_reduce(reducer);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    exec.set_work_range(range.begin() - self.m_policy.begin(),
                        range.end() - self.m_policy.begin(),
                        self.m_policy.chunk_size());
    exec.reset_steal_target();
    exec.barrier();

    long work_index = exec.get_work_index();
    typename Analysis::Reducer reducer(
        &ReducerConditional::select(self.m_functor, self.m_reducer));

    reference_type update =
        reducer.init(static_cast<pointer_type>(exec.reduce_memory()));
    while (work_index != -1) {
      const Member begin =
          static_cast<Member>(work_index) * self.m_policy.chunk_size() +
          self.m_policy.begin();
      const Member end =
          begin + self.m_policy.chunk_size() < self.m_policy.end()
              ? begin + self.m_policy.chunk_size()
              : self.m_policy.end();
      ParallelReduce::template exec_range<WorkTag>(self.m_functor, begin, end,
                                                   update);
      work_index = exec.get_work_index();
    }

    exec.fan_in_reduce(reducer);
  }

 public:
  inline void execute() const {
    if (m_policy.end() <= m_policy.begin()) {
      if (m_result_ptr) {
        typename Analysis::Reducer final_reducer(
            &ReducerConditional::select(m_functor, m_reducer));
        final_reducer.init(m_result_ptr);
        final_reducer.final(m_result_ptr);
      }
    } else {
      ThreadsExec::resize_scratch(
          Analysis::value_size(
              ReducerConditional::select(m_functor, m_reducer)),
          0);

      ThreadsExec::start(&ParallelReduce::exec, this);

      ThreadsExec::fence();

      if (m_result_ptr) {
        const pointer_type data =
            (pointer_type)ThreadsExec::root_reduce_scratch();

        const unsigned n = Analysis::value_count(
            ReducerConditional::select(m_functor, m_reducer));
        for (unsigned i = 0; i < n; ++i) {
          m_result_ptr[i] = data[i];
        }
      }
    }
  }

  template <class HostViewType>
  ParallelReduce(const FunctorType &arg_functor, const Policy &arg_policy,
                 const HostViewType &arg_result_view,
                 std::enable_if_t<Kokkos::is_view<HostViewType>::value &&
                                      !Kokkos::is_reducer<ReducerType>::value,
                                  void *> = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<HostViewType>::value,
                  "Kokkos::Threads reduce result must be a View");

    static_assert(
        std::is_same<typename HostViewType::memory_space, HostSpace>::value,
        "Kokkos::Threads reduce result must be a View in HostSpace");
  }

  inline ParallelReduce(const FunctorType &arg_functor, Policy arg_policy,
                        const ReducerType &reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {
    /*static_assert( std::is_same< typename ViewType::memory_space
                                    , Kokkos::HostSpace >::value
      , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace"
      );*/
  }
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkRange = typename Policy::WorkRange;
  using WorkTag   = typename Policy::work_tag;
  using Member    = typename Policy::member_type;
  using Analysis  = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         Policy, FunctorType>;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update, final);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelScan &self = *((const ParallelScan *)arg);

    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer final_reducer(&self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(exec.reduce_memory()));

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, false);

    //  exec.template scan_large( final_reducer );
    exec.scan_small(final_reducer);

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, true);

    exec.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsExec::start(&ParallelScan::exec, this);
    ThreadsExec::fence();
  }

  ParallelScan(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkRange = typename Policy::WorkRange;
  using WorkTag   = typename Policy::work_tag;
  using Member    = typename Policy::member_type;

  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         Policy, FunctorType>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  ReturnType &m_returnvalue;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update, final);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelScanWithTotal &self = *((const ParallelScanWithTotal *)arg);

    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer final_reducer(&self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(exec.reduce_memory()));

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, false);

    //  exec.template scan_large(final_reducer);
    exec.scan_small(final_reducer);

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, true);

    exec.fan_in();

    if (exec.pool_rank() == exec.pool_size() - 1) {
      self.m_returnvalue = update;
    }
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsExec::start(&ParallelScanWithTotal::exec, this);
    ThreadsExec::fence();
  }

  ParallelScanWithTotal(const FunctorType &arg_functor,
                        const Policy &arg_policy, ReturnType &arg_returnvalue)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_returnvalue(arg_returnvalue) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
