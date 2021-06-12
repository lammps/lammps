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

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

void test_is_specialization_of() {
  using Kokkos::Impl::is_specialization_of;
  static_assert(is_specialization_of<Kokkos::pair<float, int>, Kokkos::pair>{},
                "");
  static_assert(!is_specialization_of<Kokkos::View<int*>, Kokkos::pair>{}, "");
  static_assert(is_specialization_of<Kokkos::View<int*>, Kokkos::View>{}, "");
  // NOTE Not removing cv-qualifiers
  static_assert(!is_specialization_of<Kokkos::View<int*> const, Kokkos::View>{},
                "");
  // NOTE Would not compile because Kokkos::Array takes a non-type template
  // parameter
  // static_assert(is_specialization_of<Kokkos::Array<int, 4>, Kokkos::Array>{},
  // "");
  // But this is fine of course
  static_assert(!is_specialization_of<Kokkos::Array<float, 2>, Kokkos::pair>{},
                "");
}

template <std::size_t... Idxs, class... Args>
std::size_t do_comma_emulation_test(std::integer_sequence<std::size_t, Idxs...>,
                                    Args... args) {
  // Count the bugs, since ASSERT_EQ is a statement and not an expression
  std::size_t bugs = 0;
  // Ensure in-order evaluation
  std::size_t i = 0;
  KOKKOS_IMPL_FOLD_COMMA_OPERATOR(bugs += std::size_t(Idxs != i++) /*, ...*/);
  // Ensure expansion of multiple packs works
  KOKKOS_IMPL_FOLD_COMMA_OPERATOR(bugs += std::size_t(Idxs != args) /*, ...*/);
  return bugs;
}

TEST(utilities, comma_operator_emulation) {
  ASSERT_EQ(
      0, do_comma_emulation_test(std::make_index_sequence<5>{}, 0, 1, 2, 3, 4));
}

}  // namespace Test
