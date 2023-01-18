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

#ifndef KOKKOS_STD_ALGORITHMS_TRANSFORM_HPP
#define KOKKOS_STD_ALGORITHMS_TRANSFORM_HPP

#include "impl/Kokkos_Transform.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator first1, InputIterator last1,
          OutputIterator d_first, UnaryOperation unary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, d_first, std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator first1, InputIterator last1, OutputIterator d_first,
          UnaryOperation unary_op) {
  return Impl::transform_impl(label, ex, first1, last1, d_first,
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source), end(source), begin(dest),
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source), end(source),
                              begin(dest), std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator1 first1, InputIterator1 last1,
          InputIterator2 first2, OutputIterator d_first,
          BinaryOperation binary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
          OutputIterator d_first, BinaryOperation binary_op) {
  return Impl::transform_impl(label, ex, first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source1), end(source1), begin(source2),
                              begin(dest), std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source1), end(source1),
                              begin(source2), begin(dest),
                              std::move(binary_op));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
