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

#include <Kokkos_Core.hpp>

#define KOKKOS_PRAGMA_UNROLL(a)

namespace {

template <class Scalar, class ExecutionSpace>
struct SumPlain {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<Scalar*, execution_space> type;

  type view;

  SumPlain(type view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int /*i*/, Scalar& val) { val += Scalar(); }
};

template <class Scalar, class ExecutionSpace>
struct SumInitJoinFinalValueType {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<Scalar*, execution_space> type;
  typedef Scalar value_type;

  type view;

  SumInitJoinFinalValueType(type view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const { val = value_type(); }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& val, volatile value_type& src) const {
    val += src;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int /*i*/, value_type& val) const { val += value_type(); }
};

template <class Scalar, class ExecutionSpace>
struct SumInitJoinFinalValueType2 {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<Scalar*, execution_space> type;
  typedef Scalar value_type;

  type view;

  SumInitJoinFinalValueType2(type view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void init(volatile value_type& val) const { val = value_type(); }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& val, const volatile value_type& src) const {
    val += src;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int /*i*/, value_type& val) const { val += value_type(); }
};

template <class Scalar, class ExecutionSpace>
struct SumInitJoinFinalValueTypeArray {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<Scalar*, execution_space> type;
  typedef Scalar value_type[];

  type view;
  int n;

  SumInitJoinFinalValueTypeArray(type view_, int n_) : view(view_), n(n_) {}

  KOKKOS_INLINE_FUNCTION
  void init(value_type val) const {
    for (int k = 0; k < n; k++) {
      val[k] = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type val, const volatile value_type src) const {
    for (int k = 0; k < n; k++) {
      val[k] += src[k];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, value_type val) const {
    for (int k = 0; k < n; k++) {
      val[k] += k * i;
    }
  }
};

template <class Scalar, class ExecutionSpace>
struct SumWrongInitJoinFinalValueType {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<Scalar*, execution_space> type;
  typedef Scalar value_type;

  type view;

  SumWrongInitJoinFinalValueType(type view_) : view(view_) {}

  KOKKOS_INLINE_FUNCTION
  void init(double& val) const { val = double(); }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& val, const value_type& src) const {
    val += src;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int /*i*/, value_type& val) const { val += value_type(); }
};

template <class Scalar, class ExecutionSpace>
void TestTemplateMetaFunctions() {
  typedef typename Kokkos::View<Scalar*, ExecutionSpace> type;
  type a("A", 100);
  /*
    int sum_plain_has_init_arg = Kokkos::Impl::FunctorHasInit< SumPlain<Scalar,
    ExecutionSpace>, Scalar & >::value; ASSERT_EQ( sum_plain_has_init_arg, 0 );
    int sum_initjoinfinalvaluetype_has_init_arg = Kokkos::Impl::FunctorHasInit<
    SumInitJoinFinalValueType<Scalar, ExecutionSpace>, Scalar >::value;
    ASSERT_EQ( sum_initjoinfinalvaluetype_has_init_arg, 1 );
    int sum_initjoinfinalvaluetype_has_init_arg2 = Kokkos::Impl::FunctorHasInit<
    SumInitJoinFinalValueType2<Scalar,ExecutionSpace>, Scalar >::value;
    ASSERT_EQ( sum_initjoinfinalvaluetype_has_init_arg2, 1 );
    int sum_wronginitjoinfinalvaluetype_has_init_arg =
    Kokkos::Impl::FunctorHasInit< SumWrongInitJoinFinalValueType<Scalar,
    ExecutionSpace>, Scalar >::value; ASSERT_EQ(
    sum_wronginitjoinfinalvaluetype_has_init_arg, 0 );

    //int sum_initjoinfinalvaluetypearray_has_init_arg =
    Kokkos::Impl::FunctorHasInit< SumInitJoinFinalValueTypeArray<Scalar,
    ExecutionSpace>, Scalar[] >::value;
    //ASSERT_EQ( sum_initjoinfinalvaluetypearray_has_init_arg, 1 );

    //printf( "Values Init: %i %i %i\n", sum_plain_has_init_arg,
    sum_initjoinfinalvaluetype_has_init_arg,
    sum_wronginitjoinfinalvaluetype_has_init_arg );

    int sum_plain_has_join_arg = Kokkos::Impl::FunctorHasJoin< SumPlain<Scalar,
    ExecutionSpace>, Scalar >::value; ASSERT_EQ( sum_plain_has_join_arg, 0 );
    int sum_initjoinfinalvaluetype_has_join_arg = Kokkos::Impl::FunctorHasJoin<
    SumInitJoinFinalValueType<Scalar, ExecutionSpace>, Scalar >::value;
    ASSERT_EQ( sum_initjoinfinalvaluetype_has_join_arg, 1 );
    int sum_initjoinfinalvaluetype_has_join_arg2 = Kokkos::Impl::FunctorHasJoin<
    SumInitJoinFinalValueType2<Scalar, ExecutionSpace>, Scalar >::value;
    ASSERT_EQ( sum_initjoinfinalvaluetype_has_join_arg2, 1 );
    int sum_wronginitjoinfinalvaluetype_has_join_arg =
    Kokkos::Impl::FunctorHasJoin< SumWrongInitJoinFinalValueType<Scalar,
    ExecutionSpace>, Scalar >::value; ASSERT_EQ(
    sum_wronginitjoinfinalvaluetype_has_join_arg, 0 );

    //printf( "Values Join: %i %i %i\n", sum_plain_has_join_arg,
    sum_initjoinfinalvaluetype_has_join_arg,
    sum_wronginitjoinfinalvaluetype_has_join_arg );
  */
}

}  // namespace

namespace Test {
TEST(TEST_CATEGORY, template_meta_functions) {
  TestTemplateMetaFunctions<int, TEST_EXECSPACE>();
}
}  // namespace Test
