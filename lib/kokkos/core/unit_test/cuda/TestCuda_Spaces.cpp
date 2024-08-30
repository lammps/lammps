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

#include <Kokkos_Core.hpp>
#include <TestCuda_Category.hpp>

namespace Test {

__global__ void test_abort() { Kokkos::abort("test_abort"); }

__global__ void test_cuda_spaces_int_value(int *ptr) {
  if (*ptr == 42) {
    *ptr = 2 * 42;
  }
}

TEST(cuda, space_access) {
  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                                Kokkos::HostSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                      Kokkos::CudaHostPinnedSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::CudaSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::CudaSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::CudaUVMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                      Kokkos::CudaUVMSpace>::accessible);

  //--------------------------------------

  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                                Kokkos::CudaSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                      Kokkos::CudaUVMSpace>::assignable);

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                Kokkos::CudaSpace, Kokkos::CudaHostPinnedSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                      Kokkos::CudaHostPinnedSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                       Kokkos::HostSpace>::accessible);

  //--------------------------------------

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                      Kokkos::CudaUVMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                       Kokkos::CudaSpace>::assignable);

  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                                Kokkos::CudaSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                       Kokkos::HostSpace>::accessible);

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                Kokkos::CudaUVMSpace, Kokkos::CudaHostPinnedSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaUVMSpace,
                                      Kokkos::CudaHostPinnedSpace>::accessible);

  //--------------------------------------

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                      Kokkos::CudaHostPinnedSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                                Kokkos::HostSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                       Kokkos::CudaSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                       Kokkos::CudaSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                       Kokkos::CudaUVMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                      Kokkos::CudaUVMSpace>::accessible);

  //--------------------------------------

  static_assert(
      !Kokkos::SpaceAccessibility<Kokkos::Cuda, Kokkos::HostSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::Cuda, Kokkos::CudaSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<Kokkos::Cuda,
                                           Kokkos::CudaUVMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::Cuda,
                                 Kokkos::CudaHostPinnedSpace>::accessible);

  static_assert(!Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                            Kokkos::CudaSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                           Kokkos::CudaUVMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 Kokkos::CudaHostPinnedSpace>::accessible);

  static_assert(std::is_same<Kokkos::Impl::HostMirror<Kokkos::CudaSpace>::Space,
                             Kokkos::HostSpace>::value);

  static_assert(
      std::is_same<Kokkos::Impl::HostMirror<Kokkos::CudaUVMSpace>::Space,
                   Kokkos::Device<Kokkos::HostSpace::execution_space,
                                  Kokkos::CudaUVMSpace>>::value);

  static_assert(
      std::is_same<Kokkos::Impl::HostMirror<Kokkos::CudaHostPinnedSpace>::Space,
                   Kokkos::CudaHostPinnedSpace>::value);

  static_assert(std::is_same<Kokkos::Device<Kokkos::HostSpace::execution_space,
                                            Kokkos::CudaUVMSpace>,
                             Kokkos::Device<Kokkos::HostSpace::execution_space,
                                            Kokkos::CudaUVMSpace>>::value);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::Impl::HostMirror<Kokkos::Cuda>::Space,
                                 Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::CudaSpace>::Space,
                Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::CudaUVMSpace>::Space,
                Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::CudaHostPinnedSpace>::Space,
                Kokkos::HostSpace>::accessible);
#ifdef KOKKOS_ENABLE_CUDA_UVM
  using uvm_view = Kokkos::View<double *, Kokkos::CudaUVMSpace>;
  static_assert(std::is_same<uvm_view::HostMirror::execution_space,
                             Kokkos::DefaultHostExecutionSpace>::value,
                "Verify HostMirror execution space is really a host space");
#endif
}

TEST(cuda, uvm) {
  int *uvm_ptr = static_cast<int *>(
      Kokkos::kokkos_malloc<Kokkos::CudaUVMSpace>("uvm_ptr", sizeof(int)));

  *uvm_ptr = 42;

  Kokkos::fence();
  test_cuda_spaces_int_value<<<1, 1>>>(uvm_ptr);
  Kokkos::fence();

  EXPECT_EQ(*uvm_ptr, int(2 * 42));

  Kokkos::kokkos_free<Kokkos::CudaUVMSpace>(uvm_ptr);
}

template <class MemSpace, class ExecSpace>
struct TestViewCudaAccessible {
  enum { N = 1000 };

  using V = Kokkos::View<double *, MemSpace>;

  V m_base;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagInit &, const int i) const { m_base[i] = i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagTest &, const int i, long &error_count) const {
    if (m_base[i] != i + 1) ++error_count;
  }

  TestViewCudaAccessible() : m_base("base", N) {}

  static void run() {
    TestViewCudaAccessible self;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename MemSpace::execution_space, TagInit>(0, N),
        self);
    typename MemSpace::execution_space().fence();

    // Next access is a different execution space, must complete prior kernel.
    long error_count = -1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, TagTest>(0, N), self,
                            error_count);
    EXPECT_EQ(error_count, 0);
  }
};

TEST(cuda, impl_view_accessible) {
  TestViewCudaAccessible<Kokkos::CudaSpace, Kokkos::Cuda>::run();

  TestViewCudaAccessible<Kokkos::CudaUVMSpace, Kokkos::Cuda>::run();
  TestViewCudaAccessible<Kokkos::CudaUVMSpace,
                         Kokkos::HostSpace::execution_space>::run();

  TestViewCudaAccessible<Kokkos::CudaHostPinnedSpace, Kokkos::Cuda>::run();
  TestViewCudaAccessible<Kokkos::CudaHostPinnedSpace,
                         Kokkos::HostSpace::execution_space>::run();
}

template <class MemSpace>
struct TestViewCudaTexture {
  enum { N = 1000 };

  using V = Kokkos::View<double *, MemSpace>;
  using T = Kokkos::View<const double *, MemSpace, Kokkos::MemoryRandomAccess>;

  V m_base;
  T m_tex;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagInit &, const int i) const { m_base[i] = i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagTest &, const int i, long &error_count) const {
    if (m_tex[i] != i + 1) ++error_count;
  }

  TestViewCudaTexture() : m_base("base", N), m_tex(m_base) {}

  static void run() {
    EXPECT_TRUE((std::is_same<typename V::reference_type, double &>::value));
    EXPECT_TRUE(
        (std::is_same<typename T::reference_type, const double>::value));

    EXPECT_TRUE(V::reference_type_is_lvalue_reference);   // An ordinary view.
    EXPECT_FALSE(T::reference_type_is_lvalue_reference);  // Texture fetch
                                                          // returns by value.

    TestViewCudaTexture self;
    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::Cuda, TagInit>(0, N),
                         self);

    long error_count = -1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Kokkos::Cuda, TagTest>(0, N),
                            self, error_count);
    EXPECT_EQ(error_count, 0);
  }
};

TEST(cuda, impl_view_texture) {
  TestViewCudaTexture<Kokkos::CudaSpace>::run();
  TestViewCudaTexture<Kokkos::CudaUVMSpace>::run();
}

// couldn't create a random-access subview of a view of const T in Kokkos::Cuda
namespace issue_5594 {

template <typename View>
struct InitFunctor {
  InitFunctor(const View &view) : view_(view) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { view_(i) = i; }
  View view_;
};

template <typename V1, typename V2>
struct Issue5594Functor {
  Issue5594Functor(const V1 &v1) : v1_(v1) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, int &lerr) const {
    V2 v2(&v1_(0),
          v1_.size());  // failure here -- create subview in execution space
    lerr += v1_(i) != v2(i);  // check that subview is correct
  }
  V1 v1_;
};

template <typename View>
View create_view() {
  using execution_space = typename View::execution_space;
  View view("", 10);
  // MSVC+CUDA errors on CTAD here
  InitFunctor<View> iota(view);
  Kokkos::parallel_for("test_view_subview_const_randomaccess",
                       Kokkos::RangePolicy<execution_space>(0, view.extent(0)),
                       iota);
  return view;
}

// creating a RandomAccess subview of a view of const T in Kokkos::Cuda
template <typename Exec, typename Mem>
void test_view_subview_const_randomaccess() {
  using view_t         = Kokkos::View<int *, Mem>;
  using view_const_t   = Kokkos::View<const int *, Mem>;
  using u_view_const_t = Kokkos::View<
      const int *, Mem,
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  // create non-const view with known values
  view_t nonConst = create_view<view_t>();
  // get a const version of the values
  view_const_t view(nonConst);

  // create a subview in the execution space and check that it worked
  Issue5594Functor<view_const_t, u_view_const_t> checker(view);
  int errCount;
  Kokkos::parallel_reduce("test_view_subview_const_randomaccess",
                          Kokkos::RangePolicy<Exec>(0, view.extent(0)), checker,
                          errCount);
  EXPECT_TRUE(0 == errCount);
}
}  // namespace issue_5594

TEST(cuda, view_subview_const_randomaccess) {
  issue_5594::test_view_subview_const_randomaccess<Kokkos::Cuda,
                                                   Kokkos::CudaSpace>();
  issue_5594::test_view_subview_const_randomaccess<Kokkos::Cuda,
                                                   Kokkos::CudaUVMSpace>();
}

}  // namespace Test
