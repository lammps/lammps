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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_COPY_IF_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_COPY_IF_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class InputIterator, class OutputIterator,
          class PredicateType, class ValueType>
struct StdReplaceIfCopyFunctor {
  InputIterator m_first_from;
  OutputIterator m_first_dest;
  PredicateType m_pred;
  ValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    const auto& myvalue_from = m_first_from[i];

    if (m_pred(myvalue_from)) {
      m_first_dest[i] = m_new_value;
    } else {
      m_first_dest[i] = myvalue_from;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceIfCopyFunctor(InputIterator first_from, OutputIterator first_dest,
                          PredicateType pred, ValueType new_value)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)),
        m_new_value(std::move(new_value)) {}
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class PredicateType, class ValueType>
OutputIteratorType replace_copy_if_impl(const std::string& label,
                                        const ExecutionSpace& ex,
                                        InputIteratorType first_from,
                                        InputIteratorType last_from,
                                        OutputIteratorType first_dest,
                                        PredicateType pred,
                                        const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_t =
      StdReplaceIfCopyFunctor<index_type, InputIteratorType, OutputIteratorType,
                              PredicateType, ValueType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_t(first_from, first_dest, std::move(pred), new_value));
  ex.fence("Kokkos::replace_copy_if: fence after operation");

  // return
  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
