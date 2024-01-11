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

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>, Kokkos::Serial> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename MDRangePolicy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, CombinedFunctorReducerType, WorkTag, reference_type>;
  const iterate_type m_iter;
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
    const ReducerType& reducer     = m_iter.m_func.get_reducer();
    const size_t pool_reduce_size  = reducer.value_size();
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

    reference_type update = reducer.init(ptr);

    this->exec(update);

    reducer.final(ptr);
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const MDRangePolicy& arg_policy,
                 const ViewType& arg_result_view)
      : m_iter(arg_policy, arg_functor_reducer),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<ViewType>::value,
                  "Kokkos::Serial reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Serial reduce result must be a View accessible from "
        "HostSpace");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
