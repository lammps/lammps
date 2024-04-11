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

#ifndef KOKKOS_OPENMPTARGET_PARALLELREDUCE_RANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLELREDUCE_RANGE_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel_Common.hpp>

namespace Kokkos {
namespace Impl {

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename Policy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  static constexpr bool FunctorHasJoin = Impl::FunctorAnalysis<
      Impl::FunctorPatternInterface::REDUCE, Policy, FunctorType,
      typename ReducerType::value_type>::Reducer::has_join_member_function();
  static constexpr bool UseReducer =
      !std::is_same_v<FunctorType, typename ReducerType::functor_type>;
  static constexpr bool IsArray = std::is_pointer_v<reference_type>;

  using ParReduceSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy,
                               typename ReducerType::functor_type, pointer_type,
                               typename ReducerType::value_type>;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  bool m_result_ptr_on_device;
  const int m_result_ptr_num_elems;
  // Only let one ParallelReduce instance at a time use the scratch memory.
  // The constructor acquires the mutex which is released in the destructor.
  std::scoped_lock<std::mutex> m_scratch_memory_lock;
  using TagType = typename Policy::work_tag;

 public:
  void execute() const {
    const FunctorType& functor = m_functor_reducer.get_functor();
    if constexpr (FunctorHasJoin) {
      // Enter this loop if the Functor has a init-join.
      ParReduceSpecialize::execute_init_join(functor, m_policy, m_result_ptr,
                                             m_result_ptr_on_device);
    } else if constexpr (UseReducer) {
      // Enter this loop if the Functor is a reducer type.
      ParReduceSpecialize::execute_reducer(functor, m_policy, m_result_ptr,
                                           m_result_ptr_on_device);
    } else if constexpr (IsArray) {
      // Enter this loop if the reduction is on an array and the routine is
      // templated over the size of the array.
      if (m_result_ptr_num_elems <= 2) {
        ParReduceSpecialize::template execute_array<TagType, 2>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 4) {
        ParReduceSpecialize::template execute_array<TagType, 4>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 8) {
        ParReduceSpecialize::template execute_array<TagType, 8>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 16) {
        ParReduceSpecialize::template execute_array<TagType, 16>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 32) {
        ParReduceSpecialize::template execute_array<TagType, 32>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else {
        Kokkos::abort("array reduction length must be <= 32");
      }
    } else {
      // This loop handles the basic scalar reduction.
      ParReduceSpecialize::template execute_array<TagType, 1>(
          functor, m_policy, m_result_ptr, m_result_ptr_on_device);
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_num_elems(arg_result_view.size()),
        m_scratch_memory_lock(OpenMPTargetExec::m_mutex_scratch_ptr) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
