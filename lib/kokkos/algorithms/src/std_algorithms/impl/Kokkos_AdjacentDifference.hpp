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

#ifndef KOKKOS_STD_ALGORITHMS_ADJACENT_DIFFERENCE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ADJACENT_DIFFERENCE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ValueType1, class ValueType2, class RetType = ValueType2>
struct StdAdjacentDifferenceDefaultBinaryOpFunctor {
  KOKKOS_FUNCTION
  constexpr RetType operator()(const ValueType1& a, const ValueType2& b) const {
    return a - b;
  }
};

template <class InputIteratorType, class OutputIteratorType,
          class BinaryOperator>
struct StdAdjacentDiffFunctor {
  using index_type = typename InputIteratorType::difference_type;

  const InputIteratorType m_first_from;
  const OutputIteratorType m_first_dest;
  BinaryOperator m_op;

  KOKKOS_FUNCTION
  void operator()(const index_type i) const {
    const auto& my_value = m_first_from[i];
    if (i == 0) {
      m_first_dest[i] = my_value;
    } else {
      const auto& left_value = m_first_from[i - 1];
      m_first_dest[i]        = m_op(my_value, left_value);
    }
  }

  KOKKOS_FUNCTION
  StdAdjacentDiffFunctor(InputIteratorType first_from,
                         OutputIteratorType first_dest, BinaryOperator op)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_op(std::move(op)) {}
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOp>
OutputIteratorType adjacent_difference_impl(const std::string& label,
                                            const ExecutionSpace& ex,
                                            InputIteratorType first_from,
                                            InputIteratorType last_from,
                                            OutputIteratorType first_dest,
                                            BinaryOp bin_op) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  if (first_from == last_from) {
    return first_dest;
  }

  // aliases
  using value_type    = typename OutputIteratorType::value_type;
  using aux_view_type = ::Kokkos::View<value_type*, ExecutionSpace>;
  using functor_t =
      StdAdjacentDiffFunctor<InputIteratorType, OutputIteratorType, BinaryOp>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  aux_view_type aux_view("aux_view", num_elements);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         functor_t(first_from, first_dest, bin_op));
  ex.fence("Kokkos::adjacent_difference: fence after operation");

  // return
  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
