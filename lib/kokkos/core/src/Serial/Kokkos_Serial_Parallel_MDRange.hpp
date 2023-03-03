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

#ifndef KOKKO_SERIAL_PARALLEL_MDRANGE_HPP
#define KOKKO_SERIAL_PARALLEL_MDRANGE_HPP

#include <Kokkos_Parallel.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Serial> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

  const iterate_type m_iter;

  void exec() const {
    const typename Policy::member_type e = m_iter.m_rp.m_num_tiles;
    for (typename Policy::member_type i = 0; i < e; ++i) {
      m_iter(i);
    }
  }

 public:
  inline void execute() const { this->exec(); }
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline ParallelFor(const FunctorType& arg_functor,
                     const MDRangePolicy& arg_policy)
      : m_iter(arg_policy, arg_functor) {}
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Serial> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;

  using WorkTag = typename MDRangePolicy::work_tag;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, WorkTag,
                         void>;

  using Analysis = FunctorAnalysis<FunctorPatternInterface::REDUCE,
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

  inline void exec(reference_type update) const {
    const typename Policy::member_type e = m_iter.m_rp.m_num_tiles;
    for (typename Policy::member_type i = 0; i < e; ++i) {
      m_iter(i, update);
    }
  }

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline void execute() const {
    const size_t pool_reduce_size = Analysis::value_size(
        ReducerConditional::select(m_iter.m_func, m_reducer));
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    auto* internal_instance =
        m_iter.m_rp.space().impl_internal_space_instance();
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
        &ReducerConditional::select(m_iter.m_func, m_reducer));

    reference_type update = final_reducer.init(ptr);

    this->exec(update);

    final_reducer.final(ptr);
  }

  template <class HostViewType>
  ParallelReduce(const FunctorType& arg_functor,
                 const MDRangePolicy& arg_policy,
                 const HostViewType& arg_result_view,
                 std::enable_if_t<Kokkos::is_view<HostViewType>::value &&
                                      !Kokkos::is_reducer<ReducerType>::value,
                                  void*> = nullptr)
      : m_iter(arg_policy, arg_functor),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<HostViewType>::value,
                  "Kokkos::Serial reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename HostViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Serial reduce result must be a View in HostSpace");
  }

  inline ParallelReduce(const FunctorType& arg_functor,
                        MDRangePolicy arg_policy, const ReducerType& reducer)
      : m_iter(arg_policy, arg_functor),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {
    /*static_assert( std::is_same< typename ViewType::memory_space
                                    , Kokkos::HostSpace >::value
      , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace"
      );*/
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
