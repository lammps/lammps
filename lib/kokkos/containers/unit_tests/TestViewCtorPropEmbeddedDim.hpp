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

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>

#include <type_traits>
#include <typeinfo>

namespace Test {

namespace {

template <typename ExecSpace>
struct TestViewCtorProp_EmbeddedDim {
  using ViewIntType    = typename Kokkos::View<int**, ExecSpace>;
  using ViewDoubleType = typename Kokkos::View<double*, ExecSpace>;

  using DynRankViewIntType    = typename Kokkos::DynRankView<int, ExecSpace>;
  using DynRankViewDoubleType = typename Kokkos::DynRankView<double, ExecSpace>;

  // Cuda 7.0 has issues with using a lamda in parallel_for to initialize the
  // view - replace with this functor
  template <class ViewType>
  struct Functor {
    ViewType v;

    Functor(const ViewType& v_) : v(v_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const { v(i) = i; }
  };

  static void test_vcpt(const int N0, const int N1) {
    // Create two views to test
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
        typedef
            typename decltype(view_alloc_arg)::value_type CommonViewValueType;
        typedef typename Kokkos::View<CommonViewValueType*, ExecSpace> CVT;
        typedef typename CVT::HostMirror HostCVT;

        // Construct View using the common type; for case of specialization, an
        // 'embedded_dim' would be stored by view_alloc_arg
        CVT cv1(Kokkos::view_alloc("cv1", view_alloc_arg), N0 * N1);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, N0 * N1),
                             Functor<CVT>(cv1));

        HostCVT hcv1 = Kokkos::create_mirror_view(cv1);
        Kokkos::deep_copy(hcv1, cv1);

        ASSERT_EQ((std::is_same<CommonViewValueType, double>::value), true);
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
        typedef
            typename decltype(view_alloc_arg)::value_type CommonViewValueType;
        typedef typename Kokkos::View<CommonViewValueType*, ExecSpace> CVT;
        typedef typename CVT::HostMirror HostCVT;

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

    // Create two dynamic rank views to test
    {
      using VIT = typename TestViewCtorProp_EmbeddedDim::DynRankViewIntType;
      using VDT = typename TestViewCtorProp_EmbeddedDim::DynRankViewDoubleType;

      VIT vi1("vi1", N0, N1);
      VDT vd1("vd1", N0);

      // TEST: Test for common type between two views, one with type double,
      // other with type int Deduce common value_type and construct a view with
      // that type
      {
        // Two views
        auto view_alloc_arg = Kokkos::common_view_alloc_prop(vi1, vd1);
        typedef
            typename decltype(view_alloc_arg)::value_type CommonViewValueType;
        typedef typename Kokkos::View<CommonViewValueType*, ExecSpace> CVT;
        typedef typename CVT::HostMirror HostCVT;

        // Construct View using the common type; for case of specialization, an
        // 'embedded_dim' would be stored by view_alloc_arg
        CVT cv1(Kokkos::view_alloc("cv1", view_alloc_arg), N0 * N1);

        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, N0 * N1),
                             Functor<CVT>(cv1));

        HostCVT hcv1 = Kokkos::create_mirror_view(cv1);
        Kokkos::deep_copy(hcv1, cv1);

        ASSERT_EQ((std::is_same<CommonViewValueType, double>::value), true);
      }

      {
        // Single views
        auto view_alloc_arg = Kokkos::common_view_alloc_prop(vi1);
        typedef
            typename decltype(view_alloc_arg)::value_type CommonViewValueType;
        typedef typename Kokkos::View<CommonViewValueType*, ExecSpace> CVT;
        typedef typename CVT::HostMirror HostCVT;

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
}  // namespace Test
