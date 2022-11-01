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

#ifndef KOKKOS_STD_ALGORITHMS_PARTITION_COPY_HPP
#define KOKKOS_STD_ALGORITHMS_PARTITION_COPY_HPP

#include "impl/Kokkos_PartitionCopy.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const ExecutionSpace& ex, InputIteratorType from_first,
    InputIteratorType from_last, OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_impl(
      "Kokkos::partition_copy_iterator_api_default", ex, from_first, from_last,
      to_first_true, to_first_false, std::move(p));
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorTrueType, class OutputIteratorFalseType,
          class PredicateType>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType from_first, InputIteratorType from_last,
    OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_impl(label, ex, from_first, from_last,
                                   to_first_true, to_first_false, std::move(p));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class PredicateType>
auto partition_copy(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_impl("Kokkos::partition_copy_view_api_default",
                                   ex, cbegin(view_from), cend(view_from),
                                   begin(view_dest_true),
                                   begin(view_dest_false), std::move(p));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class PredicateType>
auto partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_impl(label, ex, cbegin(view_from),
                                   cend(view_from), begin(view_dest_true),
                                   begin(view_dest_false), std::move(p));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
