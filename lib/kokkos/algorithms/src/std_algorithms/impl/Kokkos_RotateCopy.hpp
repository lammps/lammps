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

#ifndef KOKKOS_STD_ALGORITHMS_ROTATE_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ROTATE_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class InputIterator, class OutputIterator>
struct StdRotateCopyFunctor {
  InputIterator m_first;
  InputIterator m_last;
  InputIterator m_first_n;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const {
    const IndexType shift = m_last - m_first_n;

    if (i < shift) {
      m_dest_first[i] = m_first_n[i];
    } else {
      m_dest_first[i] = m_first[i - shift];
    }
  }

  StdRotateCopyFunctor(InputIterator first, InputIterator last,
                       InputIterator first_n, OutputIterator dest_first)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_first_n(std::move(first_n)),
        m_dest_first(std::move(dest_first)) {}
};

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator rotate_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator n_first, InputIterator last,
                                OutputIterator d_first) {
  /*
    algorithm is implemented as follows:

    first 	   n_first		last
    |		      |                  |
    o  o  o  o  o  o  o  o  o  o  o  o

    dest+0 -> first_n
    dest+1 -> first_n+1
    dest+2 -> first_n+2
    dest+3 -> first
    dest+4 -> first+1
    dest+5 -> first+2
    dest+6 -> first+3
    dest+7 -> first+4
    dest+8 -> first+5
    ...
    let shift = last - first_n;

    then we have:
    if (i < shift){
      *(dest_first + i) = *(first_n + i);
    }
    else{
      *(dest_first + i) = *(from + i - shift);
    }
  */

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  if (first == last) {
    return d_first;
  }

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_type =
      StdRotateCopyFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_type(first, last, n_first, d_first));

  ex.fence("Kokkos::rotate_copy: fence after operation");

  // return
  return d_first + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
