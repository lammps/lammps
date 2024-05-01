//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <impl/Kokkos_StringManipulation.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

using Kokkos::Impl::HostSharedPtr;

namespace {

class Data {
  char d[64];

 public:
  KOKKOS_FUNCTION void write(char const* s) {
    Kokkos::Impl::strncpy(d, s, sizeof(d));
  }
};

template <class SmartPtr>
struct CheckAccessStoredPointerAndDereferenceOnDevice {
  SmartPtr m_device_ptr;
  using ElementType = typename SmartPtr::element_type;
  static_assert(std::is_same<ElementType, Data>::value);

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

// FIXME_OPENMPTARGET
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA) && \
    !defined(KOKKOS_ENABLE_OPENMPTARGET)
namespace {

struct Bar {
  double val;
};

struct Foo {
  Foo(bool allocate = false) : ptr(allocate ? new Bar : nullptr) {}
  Kokkos::Impl::HostSharedPtr<Bar> ptr;
  int use_count() { return ptr.use_count(); }
};

template <class DevMemSpace, class HostMemSpace>
void host_shared_ptr_test_reference_counting() {
  using ExecSpace = typename DevMemSpace::execution_space;
  bool is_gpu =
      !Kokkos::SpaceAccessibility<ExecSpace, Kokkos::HostSpace>::accessible;

  // Create two tracked instances
  Foo f1(true), f2(true);
  // Scope Views
  {
    Foo* fp_d_ptr =
        static_cast<Foo*>(Kokkos::kokkos_malloc<DevMemSpace>(sizeof(Foo)));
    Kokkos::View<Foo, DevMemSpace> fp_d(fp_d_ptr);
    // If using UVM or on the CPU don't make an extra HostCopy
    Foo* fp_h_ptr = std::is_same<DevMemSpace, HostMemSpace>::value
                        ? fp_d_ptr
                        : static_cast<Foo*>(
                              Kokkos::kokkos_malloc<HostMemSpace>(sizeof(Foo)));
    Kokkos::View<Foo, HostMemSpace> fp_h(fp_h_ptr);
    ASSERT_EQ(1, f1.use_count());
    ASSERT_EQ(1, f2.use_count());

    // Just for the sake of it initialize the data of the host copy
    new (fp_h.data()) Foo();
    // placement new in kernel
    //  if on GPU: should not increase use_count, fp_d will not be tracked
    //  if on Host: refcount will increase fp_d is tracked
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, 1),
        KOKKOS_LAMBDA(int) { new (fp_d.data()) Foo(f1); });
    Kokkos::fence();
    Kokkos::deep_copy(fp_h, fp_d);

    if (is_gpu)
      ASSERT_EQ(1, f1.use_count());
    else
      ASSERT_EQ(2, f1.use_count());
    ASSERT_EQ(1, f2.use_count());

    // assignment operator on host, will increase f2 use_count
    //   if default device is GPU: fp_h was untracked
    //   if default device is CPU: fp_h was tracked and use_count was 2 for
    //   aliasing f1, in which case use_count will be decreased here
    fp_h() = f2;
    ASSERT_EQ(1, f1.use_count());
    ASSERT_EQ(2, f2.use_count());

    Kokkos::deep_copy(fp_d, fp_h);
    ASSERT_EQ(1, f1.use_count());
    ASSERT_EQ(2, f2.use_count());

    // assignment in kernel:
    //  If on GPU: should not increase use_count of f1 and fp_d will not be
    //  tracked.
    //  If on Host: use_count will increase of f1, fp_d is tracked,
    //  use_count of f2 goes down.
    //  Since we are messing with the use count on the device: make host copy
    //  untracked first. Note if fp_d and fp_h alias each other (e.g. compiling
    //  for CPU only) that means fp_d() will be untracked too during assignemnt
    fp_h() = Foo();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, 1),
        KOKKOS_LAMBDA(int) { fp_d() = f1; });
    Kokkos::fence();
    Kokkos::deep_copy(fp_h, fp_d);

    if (is_gpu)
      ASSERT_EQ(1, f1.use_count());
    else
      ASSERT_EQ(2, f1.use_count());
    ASSERT_EQ(1, f2.use_count());

    // Assign non-tracked ptr
    //   if  if_gpu will not change use_count
    //   if !is_gpu will decrease use_count of f1
    fp_h() = Foo();
    ASSERT_EQ(1, f1.use_count());
    ASSERT_EQ(1, f2.use_count());
    fp_h() = f2;
    ASSERT_EQ(1, f1.use_count());
    ASSERT_EQ(2, f2.use_count());

    // before deleting host version make sure its not tracked
    fp_h() = Foo();
    if (fp_h_ptr != fp_d_ptr) Kokkos::kokkos_free<HostMemSpace>(fp_h_ptr);
    Kokkos::kokkos_free<DevMemSpace>(fp_d_ptr);
  }

  ASSERT_EQ(1, f1.use_count());
  ASSERT_EQ(1, f2.use_count());
}
}  // namespace

TEST(TEST_CATEGORY, host_shared_ptr_tracking) {
  host_shared_ptr_test_reference_counting<typename TEST_EXECSPACE::memory_space,
                                          Kokkos::HostSpace>();
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value)
    host_shared_ptr_test_reference_counting<Kokkos::CudaUVMSpace,
                                            Kokkos::CudaUVMSpace>();
#endif
#ifdef KOKKOS_ENABLE_SYCL
  if (std::is_same<TEST_EXECSPACE, Kokkos::Experimental::SYCL>::value)
    host_shared_ptr_test_reference_counting<
        Kokkos::Experimental::SYCLSharedUSMSpace,
        Kokkos::Experimental::SYCLSharedUSMSpace>();
#endif
#ifdef KOKKOS_ENABLE_HIP
  if (std::is_same<TEST_EXECSPACE, Kokkos::HIP>::value) {
    host_shared_ptr_test_reference_counting<Kokkos::HIPHostPinnedSpace,
                                            Kokkos::HIPHostPinnedSpace>();
    host_shared_ptr_test_reference_counting<Kokkos::HIPManagedSpace,
                                            Kokkos::HIPManagedSpace>();
  }
#endif
}

#endif  // KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
