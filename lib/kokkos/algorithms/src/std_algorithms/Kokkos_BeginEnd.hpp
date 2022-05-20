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

#ifndef KOKKOS_BEGIN_END_HPP
#define KOKKOS_BEGIN_END_HPP

#include <Kokkos_View.hpp>
#include "Kokkos_RandomAccessIterator.hpp"
#include "Kokkos_Constraints.hpp"

/// \file Kokkos_BeginEnd.hpp
/// \brief Kokkos begin, end, cbegin, cend

namespace Kokkos {
namespace Experimental {

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto begin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto end(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v, v.extent(0));
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cbegin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cend(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv, cv.extent(0));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
