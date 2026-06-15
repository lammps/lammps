// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_REVERSE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REVERSE_IMPL_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator>
struct StdReverseFunctor {
  using index_type = typename InputIterator::difference_type;
  static_assert(std::is_signed_v<index_type>,
                "Kokkos: StdReverseFunctor requires signed index type");

  InputIterator m_first;
  InputIterator m_last;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    ::Kokkos::kokkos_swap(m_first[i], m_last[-i - 1]);
  }

  KOKKOS_FUNCTION
  StdReverseFunctor(InputIterator first, InputIterator last)
      : m_first(std::move(first)), m_last(std::move(last)) {}
};

template <class ExecutionSpace, class InputIterator>
void reverse_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first, InputIterator last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  if (last >= first + 2) {
    // only need half
    const auto num_elements = Kokkos::Experimental::distance(first, last) / 2;
    ::Kokkos::parallel_for(label,
                           RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                           StdReverseFunctor(first, last));
    ex.fence("Kokkos::reverse: fence after operation");
  }
}

template <class TeamHandleType, class InputIterator>
KOKKOS_FUNCTION void reverse_team_impl(const TeamHandleType& teamHandle,
                                       InputIterator first,
                                       InputIterator last) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  // run
  if (last >= first + 2) {
    // only need half
    const auto num_elements = Kokkos::Experimental::distance(first, last) / 2;
    ::Kokkos::parallel_for(TeamThreadRange(teamHandle, 0, num_elements),
                           StdReverseFunctor(first, last));
    teamHandle.team_barrier();
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
