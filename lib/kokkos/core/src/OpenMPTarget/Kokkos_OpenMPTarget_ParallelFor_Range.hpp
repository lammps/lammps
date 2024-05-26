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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_FOR_RANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_FOR_RANGE_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  void execute() const { execute_impl<WorkTag>(); }

  template <class TagType>
  void execute_impl() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = m_policy.begin();
    const auto end   = m_policy.end();

    if (end <= begin) return;

    FunctorType a_functor(m_functor);

#pragma omp target teams distribute parallel for map(to : a_functor)
    for (auto i = begin; i < end; ++i) {
      if constexpr (std::is_void<TagType>::value) {
        a_functor(i);
      } else {
        a_functor(TagType(), i);
      }
    }
  }

  ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
