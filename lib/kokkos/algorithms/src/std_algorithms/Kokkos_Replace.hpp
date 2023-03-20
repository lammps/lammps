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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_HPP

#include "impl/Kokkos_Replace.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class Iterator, class ValueType>
void replace(const ExecutionSpace& ex, Iterator first, Iterator last,
             const ValueType& old_value, const ValueType& new_value) {
  return Impl::replace_impl("Kokkos::replace_iterator_api", ex, first, last,
                            old_value, new_value);
}

template <class ExecutionSpace, class Iterator, class ValueType>
void replace(const std::string& label, const ExecutionSpace& ex, Iterator first,
             Iterator last, const ValueType& old_value,
             const ValueType& new_value) {
  return Impl::replace_impl(label, ex, first, last, old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class ValueType>
void replace(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_impl("Kokkos::replace_view_api", ex, KE::begin(view),
                            KE::end(view), old_value, new_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class ValueType>
void replace(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_impl(label, ex, KE::begin(view), KE::end(view),
                            old_value, new_value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
