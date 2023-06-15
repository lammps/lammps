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

#ifndef KOKKOS_STD_ALGORITHMS_ITER_SWAP_HPP
#define KOKKOS_STD_ALGORITHMS_ITER_SWAP_HPP

#include <Kokkos_Core.hpp>
#include "impl/Kokkos_Constraints.hpp"
#include "Kokkos_Swap.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType1, class IteratorType2>
struct StdIterSwapFunctor {
  IteratorType1 m_a;
  IteratorType2 m_b;

  KOKKOS_FUNCTION
  void operator()(int i) const {
    (void)i;
    ::Kokkos::Experimental::swap(*m_a, *m_b);
  }

  KOKKOS_FUNCTION
  StdIterSwapFunctor(IteratorType1 _a, IteratorType2 _b)
      : m_a(std::move(_a)), m_b(std::move(_b)) {}
};

template <class IteratorType1, class IteratorType2>
void iter_swap_impl(IteratorType1 a, IteratorType2 b) {
  // is there a better way to do this maybe?
  ::Kokkos::parallel_for(
      1, StdIterSwapFunctor<IteratorType1, IteratorType2>(a, b));
  Kokkos::DefaultExecutionSpace().fence(
      "Kokkos::iter_swap: fence after operation");
}
}  // namespace Impl
//----------------------------------------------------------------------------

// iter_swap
template <class IteratorType1, class IteratorType2>
void iter_swap(IteratorType1 a, IteratorType2 b) {
  Impl::iter_swap_impl(a, b);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
