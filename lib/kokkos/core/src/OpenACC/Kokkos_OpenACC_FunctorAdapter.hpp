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

#ifndef KOKKOS_OPENACC_FUNCTOR_ADAPTER_HPP
#define KOKKOS_OPENACC_FUNCTOR_ADAPTER_HPP

#include <OpenACC/Kokkos_OpenACC_Macros.hpp>
#include <type_traits>

namespace Kokkos::Experimental::Impl {

enum class RoutineClause { worker, seq };

template <class Functor, class Policy, RoutineClause>
class FunctorAdapter;

#define KOKKOS_IMPL_ACC_FUNCTOR_ADAPTER(CLAUSE)                    \
  template <class Functor, class Policy>                           \
  class FunctorAdapter<Functor, Policy, RoutineClause::CLAUSE> {   \
    Functor m_functor;                                             \
    using WorkTag = typename Policy::work_tag;                     \
                                                                   \
   public:                                                         \
    FunctorAdapter(Functor const &functor) : m_functor(functor) {} \
                                                                   \
    KOKKOS_IMPL_ACC_PRAGMA(routine CLAUSE)                         \
    template <class... Args>                                       \
    KOKKOS_FUNCTION void operator()(Args &&... args) const {       \
      if constexpr (std::is_void_v<WorkTag>) {                     \
        m_functor(static_cast<Args &&>(args)...);                  \
      } else {                                                     \
        m_functor(WorkTag(), static_cast<Args &&>(args)...);       \
      }                                                            \
    }                                                              \
  }

KOKKOS_IMPL_ACC_FUNCTOR_ADAPTER(worker);
KOKKOS_IMPL_ACC_FUNCTOR_ADAPTER(seq);

#undef KOKKOS_IMPL_ACC_FUNCTOR_ADAPTER

}  // namespace Kokkos::Experimental::Impl

#endif
