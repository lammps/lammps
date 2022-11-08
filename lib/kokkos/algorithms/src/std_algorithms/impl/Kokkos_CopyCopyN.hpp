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

#ifndef KOKKOS_STD_ALGORITHMS_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class InputIterator, class OutputIterator>
struct StdCopyFunctor {
  InputIterator m_first;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(IndexType i) const { m_dest_first[i] = m_first[i]; }

  KOKKOS_FUNCTION
  StdCopyFunctor(InputIterator _first, OutputIterator _dest_first)
      : m_first(std::move(_first)), m_dest_first(std::move(_dest_first)) {}
};

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator copy_impl(const std::string& label, const ExecutionSpace& ex,
                         InputIterator first, InputIterator last,
                         OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using index_type = typename InputIterator::difference_type;
  using func_t     = StdCopyFunctor<index_type, InputIterator, OutputIterator>;

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         func_t(first, d_first));
  ex.fence("Kokkos::copy: fence after operation");

  // return
  return d_first + num_elements;
}

template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n_impl(const std::string& label, const ExecutionSpace& ex,
                           InputIterator first_from, Size count,
                           OutputIterator first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);

  if (count > 0) {
    return copy_impl(label, ex, first_from, first_from + count, first_dest);
  } else {
    return first_dest;
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
