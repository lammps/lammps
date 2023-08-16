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

#ifndef KOKKOS_STD_ALGORITHMS_FOR_EACH_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_FOR_EACH_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class UnaryFunctorType>
struct StdForEachFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  UnaryFunctorType m_functor;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_functor(m_first[i]); }

  KOKKOS_FUNCTION
  StdForEachFunctor(IteratorType _first, UnaryFunctorType _functor)
      : m_first(std::move(_first)), m_functor(std::move(_functor)) {}
};

template <class ExecutionSpace, class IteratorType, class UnaryFunctorType>
UnaryFunctorType for_each_impl(const std::string& label,
                               const ExecutionSpace& ex, IteratorType first,
                               IteratorType last, UnaryFunctorType functor) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      StdForEachFunctor<IteratorType, UnaryFunctorType>(first, functor));
  ex.fence("Kokkos::for_each: fence after operation");

  return functor;
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class UnaryFunctorType>
IteratorType for_each_n_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, SizeType n,
                             UnaryFunctorType functor) {
  auto last = first + n;
  Impl::static_assert_random_access_and_accessible(ex, first, last);
  Impl::expect_valid_range(first, last);

  if (n == 0) {
    return first;
  }

  for_each_impl(label, ex, first, last, std::move(functor));
  // no neeed to fence since for_each_impl fences already

  return last;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
