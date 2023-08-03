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

#ifndef KOKKOS_OPENACC_PARALLEL_FOR_RANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_FOR_RANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_ScheduleType.hpp>
#include <Kokkos_Parallel.hpp>

namespace Kokkos::Experimental::Impl {
template <class IndexType, class Functor>
void OpenACCParallelForRangePolicy(Schedule<Static>, int chunk_size,
                                   IndexType begin, IndexType end,
                                   Functor afunctor, int async_arg) {
  // FIXME_OPENACC FIXME_NVHPC workaround compiler bug (incorrect scope
  // analysis)
  // NVC++-S-1067-Cannot determine bounds for array - functor
  auto const functor(afunctor);
  if (chunk_size >= 1) {
// clang-format off
#pragma acc parallel loop gang(static:chunk_size) vector copyin(functor) async(async_arg)
    // clang-format on
    for (auto i = begin; i < end; ++i) {
      functor(i);
    }
  } else {
// clang-format off
#pragma acc parallel loop gang(static:*) vector copyin(functor) async(async_arg)
    // clang-format on
    for (auto i = begin; i < end; ++i) {
      functor(i);
    }
  }
}

template <class IndexType, class Functor>
void OpenACCParallelForRangePolicy(Schedule<Dynamic>, int chunk_size,
                                   IndexType begin, IndexType end,
                                   Functor afunctor, int async_arg) {
  // FIXME_OPENACC FIXME_NVHPC workaround compiler bug (incorrect scope
  // analysis)
  // NVC++-S-1067-Cannot determine bounds for array - functor
  auto const functor(afunctor);
  if (chunk_size >= 1) {
// clang-format off
#pragma acc parallel loop gang(static:chunk_size) vector copyin(functor) async(async_arg)
    // clang-format on
    for (auto i = begin; i < end; ++i) {
      functor(i);
    }
  } else {
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
    // clang-format on
    for (auto i = begin; i < end; ++i) {
      functor(i);
    }
  }
}
}  // namespace Kokkos::Experimental::Impl

template <class Functor, class... Traits>
class Kokkos::Impl::ParallelFor<Functor, Kokkos::RangePolicy<Traits...>,
                                Kokkos::Experimental::OpenACC> {
  using Policy = Kokkos::RangePolicy<Traits...>;
  Kokkos::Experimental::Impl::FunctorAdapter<
      Functor, Policy, Kokkos::Experimental::Impl::RoutineClause::seq>
      m_functor;
  Policy m_policy;
  using ScheduleType = Kokkos::Experimental::Impl::OpenACCScheduleType<Policy>;

 public:
  ParallelFor(Functor const& functor, Policy const& policy)
      : m_functor(functor), m_policy(policy) {}

  void execute() const {
    auto const begin = m_policy.begin();
    auto const end   = m_policy.end();

    if (end <= begin) {
      return;
    }

    int const async_arg  = m_policy.space().acc_async_queue();
    int const chunk_size = m_policy.chunk_size();

    Kokkos::Experimental::Impl::OpenACCParallelForRangePolicy(
        ScheduleType(), chunk_size, begin, end, m_functor, async_arg);
  }
};

#endif
