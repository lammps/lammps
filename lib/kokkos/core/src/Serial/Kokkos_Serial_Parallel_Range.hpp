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

#ifndef KOKKO_SERIAL_PARALLEL_RANGE_HPP
#define KOKKO_SERIAL_PARALLEL_RANGE_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Serial> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  std::enable_if_t<std::is_void<TagType>::value> exec() const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(i);
    }
  }

  template <class TagType>
  std::enable_if_t<!std::is_void<TagType>::value> exec() const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(t, i);
    }
  }

 public:
  inline void execute() const {
    this->template exec<typename Policy::work_tag>();
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

/*--------------------------------------------------------------------------*/

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Serial> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;

  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, WorkTag,
                         void>;

  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy, ReducerTypeFwd>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      reference_type update) const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(i, update);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      reference_type update) const {
    const TagType t{};

    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(t, i, update);
    }
  }

 public:
  inline void execute() const {
    const size_t pool_reduce_size =
        Analysis::value_size(ReducerConditional::select(m_functor, m_reducer));
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    auto* internal_instance = m_policy.space().impl_internal_space_instance();
    // Need to lock resize_thread_team_data
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    pointer_type ptr =
        m_result_ptr
            ? m_result_ptr
            : pointer_type(
                  internal_instance->m_thread_team_data.pool_reduce_local());

    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    reference_type update = final_reducer.init(ptr);

    this->template exec<WorkTag>(update);

    final_reducer.final(ptr);
  }

  template <class HostViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const HostViewType& arg_result_view,
                 std::enable_if_t<Kokkos::is_view<HostViewType>::value &&
                                      !Kokkos::is_reducer<ReducerType>::value,
                                  void*> = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<HostViewType>::value,
                  "Kokkos::Serial reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename HostViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Serial reduce result must be a View in HostSpace");
  }

  inline ParallelReduce(const FunctorType& arg_functor, Policy arg_policy,
                        const ReducerType& reducer)
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

/*--------------------------------------------------------------------------*/

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Serial> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::SCAN, Policy, FunctorType>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      reference_type update) const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(i, update, true);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      reference_type update) const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(t, i, update, true);
    }
  }

 public:
  inline void execute() const {
    const size_t pool_reduce_size  = Analysis::value_size(m_functor);
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    // Need to lock resize_thread_team_data
    auto* internal_instance = m_policy.space().impl_internal_space_instance();
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    typename Analysis::Reducer final_reducer(&m_functor);

    reference_type update = final_reducer.init(pointer_type(
        internal_instance->m_thread_team_data.pool_reduce_local()));

    this->template exec<WorkTag>(update);
  }

  inline ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

/*--------------------------------------------------------------------------*/
template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Serial> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::SCAN, Policy, FunctorType>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  ReturnType& m_returnvalue;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      reference_type update) const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(i, update, true);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      reference_type update) const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(t, i, update, true);
    }
  }

 public:
  inline void execute() {
    const size_t pool_reduce_size  = Analysis::value_size(m_functor);
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    // Need to lock resize_thread_team_data
    auto* internal_instance = m_policy.space().impl_internal_space_instance();
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    typename Analysis::Reducer final_reducer(&m_functor);

    reference_type update = final_reducer.init(pointer_type(
        internal_instance->m_thread_team_data.pool_reduce_local()));

    this->template exec<WorkTag>(update);

    m_returnvalue = update;
  }

  inline ParallelScanWithTotal(const FunctorType& arg_functor,
                               const Policy& arg_policy,
                               ReturnType& arg_returnvalue)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_returnvalue(arg_returnvalue) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
