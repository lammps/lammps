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

#ifndef KOKKOS_STD_ALGORITHMS_UNIQUE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_UNIQUE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Move.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_AdjacentFind.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class InputIt, class OutputIt,
          class BinaryPredicateType>
struct StdUniqueFunctor {
  InputIt m_first_from;
  InputIt m_last_from;
  OutputIt m_first_dest;
  BinaryPredicateType m_pred;

  KOKKOS_FUNCTION
  StdUniqueFunctor(InputIt first_from, InputIt last_from, OutputIt first_dest,
                   BinaryPredicateType pred)
      : m_first_from(std::move(first_from)),
        m_last_from(std::move(last_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, IndexType& update,
                  const bool final_pass) const {
    auto& val_i         = m_first_from[i];
    const auto& val_ip1 = m_first_from[i + 1];

    if (final_pass) {
      if (!m_pred(val_i, val_ip1)) {
        m_first_dest[update] = std::move(val_i);
      }
    }

    if (!m_pred(val_i, val_ip1)) {
      update += 1;
    }
  }
};

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType unique_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements == 0) {
    return first;
  } else if (num_elements == 1) {
    return last;
  } else {
    // ----------
    // step 1:
    // find first location of adjacent equal elements
    // ----------
    auto it_found =
        ::Kokkos::Experimental::adjacent_find(ex, first, last, pred);

    // if none, all elements are unique, so nothing to do
    if (it_found == last) {
      return last;
    } else {
      // if here, we found some equal adjacent elements,
      // so count all preceeding unique elements
      const auto num_unique_found_in_step_one = it_found - first;

      // ----------
      // step 2:
      // ----------
      // since we found some unique elements, we don't need to explore
      // the full range [first, last), but only need to focus on the
      // remaining range [it_found, last)
      const auto num_elements_to_explore = last - it_found;

      // create a tmp view to use to *move* all unique elements
      // using the same algorithm used for unique_copy but we now move things
      using value_type    = typename IteratorType::value_type;
      using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
      tmp_view_type tmp_view("std_unique_tmp_view", num_elements_to_explore);

      // scan extent is: num_elements_to_explore - 1
      // for same reason as the one explained in unique_copy
      const auto scan_size = num_elements_to_explore - 1;
      auto tmp_first       = ::Kokkos::Experimental::begin(tmp_view);
      using output_it      = decltype(tmp_first);

      using index_type = typename IteratorType::difference_type;
      using func_type =
          StdUniqueFunctor<index_type, IteratorType, output_it, PredicateType>;
      index_type count = 0;
      ::Kokkos::parallel_scan(
          label, RangePolicy<ExecutionSpace>(ex, 0, scan_size),
          func_type(it_found, last, tmp_first, pred), count);

      // move last element too, for the same reason as the unique_copy
      auto unused_r =
          Impl::move_impl("Kokkos::move_from_unique", ex, it_found + scan_size,
                          last, tmp_first + count);
      (void)unused_r;  // r1 not used

      // ----------
      // step 3
      // ----------
      // move back from tmp to original range,
      // ensuring we start overwriting after the original unique found
      using tmp_readwrite_iterator_type = decltype(begin(tmp_view));
      using step3_func_t =
          StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;

      ::Kokkos::parallel_for(
          "unique_step3_parfor",
          RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
          step3_func_t(begin(tmp_view),
                       (first + num_unique_found_in_step_one)));

      ex.fence("Kokkos::unique: fence after operation");

      // return iterator to one passed the last written
      // (the +1 is needed to account for the last element, see above)
      return (first + num_unique_found_in_step_one + count + 1);
    }
  }
}

template <class ExecutionSpace, class IteratorType>
IteratorType unique_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last) {
  using value_type    = typename IteratorType::value_type;
  using binary_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return unique_impl(label, ex, first, last, binary_pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
