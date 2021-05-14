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
#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

using Kokkos::Impl::HostSharedPtr;

namespace {

class Data {
  Kokkos::Array<char, 64> d;

 public:
  KOKKOS_FUNCTION void write(char const* c) {
    for (int i = 0; i < 64 && c; ++i, ++c) {
      d[i] = *c;
    }
  }
};

template <class SmartPtr>
struct CheckAccessStoredPointerAndDereferenceOnDevice {
  SmartPtr m_device_ptr;
  using ElementType = typename SmartPtr::element_type;
  static_assert(std::is_same<ElementType, Data>::value, "");

  CheckAccessStoredPointerAndDereferenceOnDevice(SmartPtr device_ptr)
      : m_device_ptr(device_ptr) {
    int errors;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1), *this,
                            errors);
    EXPECT_EQ(errors, 0);
  }

  KOKKOS_FUNCTION void operator()(int, int& e) const {
    auto raw_ptr = m_device_ptr.get();  // get

    auto tmp = new (raw_ptr) ElementType();

    auto& obj = *m_device_ptr;  // operator*
    if (&obj != raw_ptr) ++e;

    m_device_ptr->write("hello world");  // operator->

    tmp->~ElementType();
  }
};

template <class Ptr>
CheckAccessStoredPointerAndDereferenceOnDevice<Ptr>
check_access_stored_pointer_and_dereference_on_device(Ptr p) {
  return {p};
}

template <class SmartPtr>
struct CheckSpecialMembersOnDevice {
  SmartPtr m_device_ptr;

  KOKKOS_FUNCTION void operator()(int, int& e) const {
    SmartPtr p1 = m_device_ptr;   // copy construction
    SmartPtr p2 = std::move(p1);  // move construction

    p1 = p2;             // copy assignment
    p2 = std::move(p1);  // move assignment

    SmartPtr p3;  // default constructor
    if (p3) ++e;
    SmartPtr p4{nullptr};
    if (p4) ++e;
  }

  CheckSpecialMembersOnDevice(SmartPtr device_ptr) : m_device_ptr(device_ptr) {
    int errors;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1), *this,
                            errors);
    EXPECT_EQ(errors, 0);
  }
};

template <class Ptr>
CheckSpecialMembersOnDevice<Ptr> check_special_members_on_device(Ptr p) {
  return {p};
}

}  // namespace

TEST(TEST_CATEGORY, host_shared_ptr_dereference_on_device) {
  using T = Data;

  using MemorySpace = TEST_EXECSPACE::memory_space;

  HostSharedPtr<T> device_ptr(
      static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(sizeof(T))),
      [](T* p) { Kokkos::kokkos_free<MemorySpace>(p); });

  check_access_stored_pointer_and_dereference_on_device(device_ptr);
}

// FIXME_OPENMPTARGET
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, host_shared_ptr_special_members_on_device) {
  using T = Data;

  using MemorySpace = TEST_EXECSPACE::memory_space;

  HostSharedPtr<T> device_ptr(
      static_cast<T*>(Kokkos::kokkos_malloc<MemorySpace>(sizeof(T))),
      [](T* p) { Kokkos::kokkos_free<MemorySpace>(p); });

  check_special_members_on_device(device_ptr);
}
#endif
