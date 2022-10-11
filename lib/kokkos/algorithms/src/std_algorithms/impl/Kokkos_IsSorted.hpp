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

#ifndef KOKKOS_STD_ALGORITHMS_IS_SORTED_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_IS_SORTED_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ComparatorType>
struct StdIsSortedFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  ComparatorType m_comparator;

  KOKKOS_FUNCTION
  void operator()(const index_type i, std::size_t& update) const {
    const auto& val_i   = m_first[i];
    const auto& val_ip1 = m_first[i + 1];

    if (m_comparator(val_ip1, val_i)) {
      ++update;
    }
  }

  KOKKOS_FUNCTION
  StdIsSortedFunctor(IteratorType _first1, ComparatorType comparator)
      : m_first(std::move(_first1)), m_comparator(std::move(comparator)) {}
};

template <class ExecutionSpace, class IteratorType, class ComparatorType>
bool is_sorted_impl(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements <= 1) {
    return true;
  }

  // use num_elements-1 because each index handles i and i+1
  const auto num_elements_minus_one = num_elements - 1;
  using functor_type = StdIsSortedFunctor<IteratorType, ComparatorType>;

  // result is incremented by one if sorting breaks at index i
  std::size_t result = 0;
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_minus_one),
      functor_type(first, std::move(comp)), result);

  return result == 0;
}

template <class ExecutionSpace, class IteratorType>
bool is_sorted_impl(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last) {
  using value_type = typename IteratorType::value_type;
  using pred_t     = Impl::StdAlgoLessThanBinaryPredicate<value_type>;
  return is_sorted_impl(label, ex, first, last, pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
