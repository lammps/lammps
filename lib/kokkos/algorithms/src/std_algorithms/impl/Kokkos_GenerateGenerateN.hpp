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

#ifndef KOKKOS_STD_ALGORITHMS_GENERATE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_GENERATE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class Generator>
struct StdGenerateFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  Generator m_generator;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_first[i] = m_generator(); }

  KOKKOS_FUNCTION
  StdGenerateFunctor(IteratorType _first, Generator _g)
      : m_first(std::move(_first)), m_generator(std::move(_g)) {}
};

template <class ExecutionSpace, class IteratorType, class Generator>
void generate_impl(const std::string& label, const ExecutionSpace& ex,
                   IteratorType first, IteratorType last, Generator g) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdGenerateFunctor<IteratorType, Generator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, g));
  ex.fence("Kokkos::generate: fence after operation");
}

template <class ExecutionSpace, class IteratorType, class Size, class Generator>
IteratorType generate_n_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, Size count, Generator g) {
  if (count <= 0) {
    return first;
  }

  generate_impl(label, ex, first, first + count, g);
  return first + count;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
