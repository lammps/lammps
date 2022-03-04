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

#include <impl/Kokkos_HostSharedPtr.hpp>

#include <gtest/gtest.h>

using Kokkos::Impl::HostSharedPtr;

TEST(TEST_CATEGORY, host_shared_ptr_use_count) {
  using T = int;
  {
    HostSharedPtr<T> p1;
    EXPECT_EQ(p1.use_count(), 0);
  }
  {
    HostSharedPtr<T> p1(nullptr);
    EXPECT_EQ(p1.use_count(), 0);
  }
  {
    HostSharedPtr<T> p1(new T());
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T(), [](T* p) { delete p; });
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2(p1);  // copy construction
    EXPECT_EQ(p1.use_count(), 2);
    EXPECT_EQ(p2.use_count(), 2);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2(std::move(p1));  // move construction
    EXPECT_EQ(p2.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2;
    p2 = p1;  // copy assignment
    EXPECT_EQ(p1.use_count(), 2);
    EXPECT_EQ(p2.use_count(), 2);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2;
    p2 = std::move(p1);  // move assignment
    EXPECT_EQ(p2.use_count(), 1);
  }
}

TEST(TEST_CATEGORY, host_shared_ptr_get) {
  using T = int;
  {
    HostSharedPtr<T> p1;
    EXPECT_EQ(p1.get(), nullptr);
  }
  {
    HostSharedPtr<T> p1(nullptr);
    EXPECT_EQ(p1.get(), nullptr);
  }
  {
    T* p_i = new T();
    HostSharedPtr<T> p1(p_i);
    EXPECT_EQ(p1.get(), p_i);
  }
  {
    T* p_i = new T();
    HostSharedPtr<T> p1(p_i, [](T* p) { delete p; });
    EXPECT_EQ(p1.get(), p_i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    EXPECT_EQ(p1.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2(p1);  // copy construction
    EXPECT_EQ(p1.get(), &i);
    EXPECT_EQ(p1.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2(std::move(p1));  // move construction
    EXPECT_EQ(p1.get(), nullptr);
    EXPECT_EQ(p2.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2;
    p2 = p1;  // copy assignment
    EXPECT_EQ(p1.get(), &i);
    EXPECT_EQ(p1.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2;
    p2 = std::move(p1);  // move assignment
    EXPECT_EQ(p1.get(), nullptr);
    EXPECT_EQ(p2.get(), &i);
  }
}
