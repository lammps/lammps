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

#include <Kokkos_DetectionIdiom.hpp>

#define STATIC_ASSERT(cond) static_assert(cond, "");

void test_nonesuch() {
  using Kokkos::nonesuch;
  STATIC_ASSERT(!std::is_constructible<nonesuch>::value);
  STATIC_ASSERT(!std::is_destructible<nonesuch>::value);
  STATIC_ASSERT(!std::is_copy_constructible<nonesuch>::value);
  STATIC_ASSERT(!std::is_move_constructible<nonesuch>::value);
#ifdef KOKKOS_ENABLE_CXX17
  STATIC_ASSERT(!std::is_aggregate<nonesuch>::value);
#endif
}

#undef STATIC_ASSERT

namespace Example {
// Example from https://en.cppreference.com/w/cpp/experimental/is_detected
template <class T>
using copy_assign_t = decltype(std::declval<T&>() = std::declval<const T&>());

struct Meow {};
struct Purr {
  void operator=(const Purr&) = delete;
};

static_assert(Kokkos::is_detected<copy_assign_t, Meow>::value,
              "Meow should be copy assignable!");
static_assert(!Kokkos::is_detected<copy_assign_t, Purr>::value,
              "Purr should not be copy assignable!");
static_assert(Kokkos::is_detected_exact<Meow&, copy_assign_t, Meow>::value,
              "Copy assignment of Meow should return Meow&!");

template <class T>
using diff_t = typename T::difference_type;

template <class Ptr>
using difference_type = Kokkos::detected_or_t<std::ptrdiff_t, diff_t, Ptr>;

struct Woof {
  using difference_type = int;
};
struct Bark {};

static_assert(std::is_same<difference_type<Woof>, int>::value,
              "Woof's difference_type should be int!");
static_assert(std::is_same<difference_type<Bark>, std::ptrdiff_t>::value,
              "Bark's difference_type should be ptrdiff_t!");
}  // namespace Example
