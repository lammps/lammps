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

#ifndef KOKKOS_STD_ALGORITHMS_ROTATE_HPP
#define KOKKOS_STD_ALGORITHMS_ROTATE_HPP

#include "impl/Kokkos_Rotate.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType>
IteratorType rotate(const ExecutionSpace& ex, IteratorType first,
                    IteratorType n_first, IteratorType last) {
  return Impl::rotate_impl("Kokkos::rotate_iterator_api_default", ex, first,
                           n_first, last);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType n_first,
                    IteratorType last) {
  return Impl::rotate_impl(label, ex, first, n_first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto rotate(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_impl("Kokkos::rotate_view_api_default", ex, begin(view),
                           begin(view) + n_location, end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto rotate(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_impl(label, ex, begin(view), begin(view) + n_location,
                           end(view));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
