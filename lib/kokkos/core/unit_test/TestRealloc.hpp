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
#ifndef TESTREALLOC_HPP_
#define TESTREALLOC_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace TestViewRealloc {

struct Default {};
struct WithoutInitializing {};

template <typename View, typename... Args>
inline void realloc_dispatch(Default, View& v, Args&&... args) {
  Kokkos::realloc(v, std::forward<Args>(args)...);
}

template <typename View, typename... Args>
inline void realloc_dispatch(WithoutInitializing, View& v, Args&&... args) {
  Kokkos::realloc(Kokkos::WithoutInitializing, v, std::forward<Args>(args)...);
}

template <class DeviceType, class Tag = Default>
void impl_testRealloc() {
  const size_t sizes[8] = {2, 3, 4, 5, 6, 7, 8, 9};

  // Check #904 fix (no reallocation if dimensions didn't change).
  {
    using view_type = Kokkos::View<int*, DeviceType>;
    view_type view_1d("view_1d", sizes[0]);
    const int* oldPointer = view_1d.data();
    auto const& oldLabel  = view_1d.label();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_1d, sizes[0]);
    auto const& newLabel = view_1d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_1d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int**, DeviceType>;
    view_type view_2d("view_2d", sizes[0], sizes[1]);
    auto const& oldLabel  = view_2d.label();
    const int* oldPointer = view_2d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_2d, sizes[0], sizes[1]);
    auto const& newLabel = view_2d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_2d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int***, DeviceType>;
    view_type view_3d("view_3d", sizes[0], sizes[1], sizes[2]);
    auto const& oldLabel  = view_3d.label();
    const int* oldPointer = view_3d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_3d, sizes[0], sizes[1], sizes[2]);
    auto const& newLabel = view_3d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_3d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int****, DeviceType>;
    view_type view_4d("view_4d", sizes[0], sizes[1], sizes[2], sizes[3]);
    auto const& oldLabel  = view_4d.label();
    const int* oldPointer = view_4d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_4d, sizes[0], sizes[1], sizes[2], sizes[3]);
    auto const& newLabel = view_4d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_4d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*****, DeviceType>;
    view_type view_5d("view_5d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4]);
    auto const& oldLabel  = view_5d.label();
    const int* oldPointer = view_5d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_5d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4]);
    auto const& newLabel = view_5d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_5d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int******, DeviceType>;
    view_type view_6d("view_6d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5]);
    const int* oldPointer = view_6d.data();
    auto const& oldLabel  = view_6d.label();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_6d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5]);
    auto const& newLabel = view_6d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_6d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*******, DeviceType>;
    view_type view_7d("view_7d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6]);
    auto const& oldLabel  = view_7d.label();
    const int* oldPointer = view_7d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_7d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5], sizes[6]);
    auto const& newLabel = view_7d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_7d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int********, DeviceType>;
    view_type view_8d("view_8d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6], sizes[7]);
    auto const& oldLabel  = view_8d.label();
    const int* oldPointer = view_8d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_8d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5], sizes[6], sizes[7]);
    auto const& newLabel = view_8d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_8d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
}

template <class DeviceType>
void testRealloc() {
  {
    impl_testRealloc<DeviceType>();  // with data initialization
  }
  {
    impl_testRealloc<DeviceType,
                     WithoutInitializing>();  // without data initialization
  }
}

}  // namespace TestViewRealloc
#endif  // TESTREALLOC_HPP_
