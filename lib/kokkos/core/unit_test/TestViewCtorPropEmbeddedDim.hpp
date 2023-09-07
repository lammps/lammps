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

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <type_traits>
#include <typeinfo>

namespace Test {

namespace {

template <typename ExecSpace>
struct TestViewCtorProp_EmbeddedDim {
  using ViewIntType    = typename Kokkos::View<int**, ExecSpace>;
  using ViewDoubleType = typename Kokkos::View<double*, ExecSpace>;

  // Cuda 7.0 has issues with using a lambda in parallel_for to initialize the
  // view - replace with this functor
  template <class ViewType>
  struct Functor {
    ViewType v;

    Functor(const ViewType& v_) : v(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const { v(i) = i; }
  };

  static void test_vcpt(const int N0, const int N1) {
    // Create views to test
    {
      using VIT = typename TestViewCtorProp_EmbeddedDim::ViewIntType;
      using VDT = typename TestViewCtorProp_EmbeddedDim::ViewDoubleType;

      VIT vi1("vi1", N0, N1);
      VDT vd1("vd1", N0);

      // TEST: Test for common type between two views, one with type double,
      // other with type int Deduce common value_type and construct a view with
      // that type
      {
        // Two views
        auto view_alloc_arg = Kokkos::common_view_alloc_prop(vi1, vd1);
        using CommonViewValueType =
            typename decltype(view_alloc_arg)::value_type;
        using CVT     = typename Kokkos::View<CommonViewValueType*, ExecSpace>;
        using HostCVT = typename CVT::HostMirror;

        // Construct View using the common type; for case of specialization, an
        // 'embedded_dim' would be stored by view_alloc_arg
        CVT cv1(Kokkos::view_alloc("cv1", view_alloc_arg), N0 * N1);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, N0 * N1),
                             Functor<CVT>(cv1));

        HostCVT hcv1 = Kokkos::create_mirror_view(cv1);
        Kokkos::deep_copy(hcv1, cv1);

        ASSERT_EQ((std::is_same<CommonViewValueType, double>::value), true);
        ASSERT_EQ(
            (std::is_same<typename decltype(view_alloc_arg)::scalar_array_type,
                          CommonViewValueType>::value),
            true);
#if 0
      // debug output
      for ( int i = 0; i < N0*N1; ++i ) {
        printf(" Output check: hcv1(%d) = %lf\n ", i, hcv1(i) );
      }

      printf( " Common value type view: %s \n", typeid( CVT() ).name() );
      printf( " Common value type: %s \n", typeid( CommonViewValueType() ).name() );
      if ( std::is_same< CommonViewValueType, double >::value == true ) {
        printf("Proper common value_type\n");
      }
      else {
        printf("WRONG common value_type\n");
      }
      // end debug output
#endif
      }

      {
        // Single view
        auto view_alloc_arg = Kokkos::common_view_alloc_prop(vi1);
        using CommonViewValueType =
            typename decltype(view_alloc_arg)::value_type;
        using CVT     = typename Kokkos::View<CommonViewValueType*, ExecSpace>;
        using HostCVT = typename CVT::HostMirror;

        // Construct View using the common type; for case of specialization, an
        // 'embedded_dim' would be stored by view_alloc_arg
        CVT cv1(Kokkos::view_alloc("cv1", view_alloc_arg), N0 * N1);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, N0 * N1),
                             Functor<CVT>(cv1));

        HostCVT hcv1 = Kokkos::create_mirror_view(cv1);
        Kokkos::deep_copy(hcv1, cv1);

        ASSERT_EQ((std::is_same<CommonViewValueType, int>::value), true);
      }
    }

  }  // end test_vcpt

};  // end struct

}  // namespace

TEST(TEST_CATEGORY, viewctorprop_embedded_dim) {
  TestViewCtorProp_EmbeddedDim<TEST_EXECSPACE>::test_vcpt(2, 3);
}

TEST(TEST_CATEGORY,
     viewctorpop_view_allocate_without_initializing_backward_compatility) {
  using deprecated_view_alloc = Kokkos::ViewAllocateWithoutInitializing;
  Kokkos::View<int**, TEST_EXECSPACE> v(deprecated_view_alloc("v"), 5, 7);
}

}  // namespace Test
