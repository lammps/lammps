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

#ifndef KOKKOS_STD_ALGORITHMS_REVERSE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REVERSE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_Swap.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator>
struct StdReverseFunctor {
  using index_type = typename InputIterator::difference_type;
  static_assert(std::is_signed<index_type>::value,
                "Kokkos: StdReverseFunctor requires signed index type");

  InputIterator m_first;
  InputIterator m_last;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    // the swap below is doing the same thing, but
    // for Intel 18.0.5 does not work.
    // But putting the impl directly here, it works.
#ifdef KOKKOS_COMPILER_INTEL
    typename InputIterator::value_type tmp = std::move(m_first[i]);
    m_first[i]                             = std::move(m_last[-i - 1]);
    m_last[-i - 1]                         = std::move(tmp);
#else
    ::Kokkos::Experimental::swap(m_first[i], m_last[-i - 1]);
#endif
  }

  StdReverseFunctor(InputIterator first, InputIterator last)
      : m_first(std::move(first)), m_last(std::move(last)) {}
};

template <class ExecutionSpace, class InputIterator>
void reverse_impl(const std::string& label, const ExecutionSpace& ex,
                  InputIterator first, InputIterator last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // aliases
  using func_t = StdReverseFunctor<InputIterator>;

  // run
  if (last >= first + 2) {
    // only need half
    const auto num_elements = Kokkos::Experimental::distance(first, last) / 2;
    ::Kokkos::parallel_for(label,
                           RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                           func_t(first, last));
    ex.fence("Kokkos::reverse: fence after operation");
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
