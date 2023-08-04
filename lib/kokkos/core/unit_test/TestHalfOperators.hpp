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

#ifndef TESTHALFOPERATOR_HPP_
#define TESTHALFOPERATOR_HPP_
namespace Test {
using namespace Kokkos::Experimental;
using ExecutionSpace = TEST_EXECSPACE;
using ScalarType     = double;
using ViewType       = Kokkos::View<ScalarType*, ExecutionSpace>;
using ViewTypeHost   = Kokkos::View<ScalarType*, Kokkos::HostSpace>;
KOKKOS_FUNCTION
const half_t& accept_ref(const half_t& a) { return a; }
KOKKOS_FUNCTION
double accept_ref_expected(const half_t& a) {
  double tmp = static_cast<double>(a);
  return tmp;
}
#if !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_FUNCTION
const bhalf_t& accept_ref(const bhalf_t& a) { return a; }
KOKKOS_FUNCTION
double accept_ref_expected(const bhalf_t& a) {
  double tmp = static_cast<double>(a);
  return tmp;
}
#endif  // !KOKKOS_BHALF_T_IS_FLOAT

enum OP_TESTS {
  ASSIGN,
  ASSIGN_CHAINED,
  UNA,
  UNS,
  PREFIX_INC,
  PREFIX_DEC,
  POSTFIX_INC,
  POSTFIX_DEC,
  CADD_H_H,
  CADD_H_S,
  CADD_S_H,
  CADD_H_D,
  CADD_D_H,
  CSUB_H_H,
  CSUB_H_S,
  CSUB_S_H,
  CSUB_H_D,
  CSUB_D_H,
  CMUL_H_H,
  CMUL_H_S,
  CMUL_S_H,
  CMUL_H_D,
  CMUL_D_H,
  CDIV_H_H,
  CDIV_H_S,
  CDIV_S_H,
  CDIV_H_D,
  CDIV_D_H,
  ADD_H_H,
  ADD_H_S,
  ADD_S_H,
  ADD_H_D,
  ADD_D_H,
  ADD_H_H_SZ,
  ADD_H_S_SZ,
  ADD_S_H_SZ,
  ADD_H_D_SZ,
  ADD_D_H_SZ,
  ADD_SI_H,
  ADD_SI_H_SZ,
  ADD_I_H,
  ADD_I_H_SZ,
  ADD_LI_H,
  ADD_LI_H_SZ,
  ADD_LLI_H,
  ADD_LLI_H_SZ,
  ADD_USI_H,
  ADD_USI_H_SZ,
  ADD_UI_H,
  ADD_UI_H_SZ,
  ADD_ULI_H,
  ADD_ULI_H_SZ,
  ADD_ULLI_H,
  ADD_ULLI_H_SZ,
  ADD_H_SI,
  ADD_H_SI_SZ,
  ADD_H_I,
  ADD_H_I_SZ,
  ADD_H_LI,
  ADD_H_LI_SZ,
  ADD_H_LLI,
  ADD_H_LLI_SZ,
  ADD_H_USI,
  ADD_H_USI_SZ,
  ADD_H_UI,
  ADD_H_UI_SZ,
  ADD_H_ULI,
  ADD_H_ULI_SZ,
  ADD_H_ULLI,
  ADD_H_ULLI_SZ,
  SUB_H_H,
  SUB_H_S,
  SUB_S_H,
  SUB_H_D,
  SUB_D_H,
  SUB_H_H_SZ,
  SUB_H_S_SZ,
  SUB_S_H_SZ,
  SUB_H_D_SZ,
  SUB_D_H_SZ,
  SUB_SI_H,
  SUB_SI_H_SZ,
  SUB_I_H,
  SUB_I_H_SZ,
  SUB_LI_H,
  SUB_LI_H_SZ,
  SUB_LLI_H,
  SUB_LLI_H_SZ,
  SUB_USI_H,
  SUB_USI_H_SZ,
  SUB_UI_H,
  SUB_UI_H_SZ,
  SUB_ULI_H,
  SUB_ULI_H_SZ,
  SUB_ULLI_H,
  SUB_ULLI_H_SZ,
  SUB_H_SI,
  SUB_H_SI_SZ,
  SUB_H_I,
  SUB_H_I_SZ,
  SUB_H_LI,
  SUB_H_LI_SZ,
  SUB_H_LLI,
  SUB_H_LLI_SZ,
  SUB_H_USI,
  SUB_H_USI_SZ,
  SUB_H_UI,
  SUB_H_UI_SZ,
  SUB_H_ULI,
  SUB_H_ULI_SZ,
  SUB_H_ULLI,
  SUB_H_ULLI_SZ,
  MUL_H_H,
  MUL_H_S,
  MUL_S_H,
  MUL_H_D,
  MUL_D_H,
  MUL_H_H_SZ,
  MUL_H_S_SZ,
  MUL_S_H_SZ,
  MUL_H_D_SZ,
  MUL_D_H_SZ,
  MUL_SI_H,
  MUL_SI_H_SZ,
  MUL_I_H,
  MUL_I_H_SZ,
  MUL_LI_H,
  MUL_LI_H_SZ,
  MUL_LLI_H,
  MUL_LLI_H_SZ,
  MUL_USI_H,
  MUL_USI_H_SZ,
  MUL_UI_H,
  MUL_UI_H_SZ,
  MUL_ULI_H,
  MUL_ULI_H_SZ,
  MUL_ULLI_H,
  MUL_ULLI_H_SZ,
  MUL_H_SI,
  MUL_H_SI_SZ,
  MUL_H_I,
  MUL_H_I_SZ,
  MUL_H_LI,
  MUL_H_LI_SZ,
  MUL_H_LLI,
  MUL_H_LLI_SZ,
  MUL_H_USI,
  MUL_H_USI_SZ,
  MUL_H_UI,
  MUL_H_UI_SZ,
  MUL_H_ULI,
  MUL_H_ULI_SZ,
  MUL_H_ULLI,
  MUL_H_ULLI_SZ,
  DIV_H_H,
  DIV_H_S,
  DIV_S_H,
  DIV_H_D,
  DIV_D_H,
  DIV_H_H_SZ,
  DIV_H_S_SZ,
  DIV_S_H_SZ,
  DIV_H_D_SZ,
  DIV_D_H_SZ,
  DIV_SI_H,
  DIV_SI_H_SZ,
  DIV_I_H,
  DIV_I_H_SZ,
  DIV_LI_H,
  DIV_LI_H_SZ,
  DIV_LLI_H,
  DIV_LLI_H_SZ,
  DIV_USI_H,
  DIV_USI_H_SZ,
  DIV_UI_H,
  DIV_UI_H_SZ,
  DIV_ULI_H,
  DIV_ULI_H_SZ,
  DIV_ULLI_H,
  DIV_ULLI_H_SZ,
  DIV_H_SI,
  DIV_H_SI_SZ,
  DIV_H_I,
  DIV_H_I_SZ,
  DIV_H_LI,
  DIV_H_LI_SZ,
  DIV_H_LLI,
  DIV_H_LLI_SZ,
  DIV_H_USI,
  DIV_H_USI_SZ,
  DIV_H_UI,
  DIV_H_UI_SZ,
  DIV_H_ULI,
  DIV_H_ULI_SZ,
  DIV_H_ULLI,
  DIV_H_ULLI_SZ,
  NEG,
  AND,
  OR,
  EQ,
  NEQ,
  LT,
  GT,
  LE,
  GE,  // TODO: TW,
  PASS_BY_REF,
  AO_IMPL_HALF,
  AO_HALF_T,
  N_OP_TESTS
};

// volatile-qualified parameter type 'volatile half_type' is deprecated
#if !defined(KOKKOS_ENABLE_CXX20) && !defined(KOKKOS_ENABLE_CXX23)
template <class view_type, class half_type>
struct Functor_TestHalfVolatileOperators {
  volatile half_type h_lhs, h_rhs;
  view_type actual_lhs, expected_lhs;
  double d_lhs, d_rhs;
  Functor_TestHalfVolatileOperators(volatile half_type lhs = half_type(0),
                                    volatile half_type rhs = half_type(0))
      : h_lhs(lhs), h_rhs(rhs) {
    actual_lhs   = view_type("actual_lhs", N_OP_TESTS);
    expected_lhs = view_type("expected_lhs", N_OP_TESTS);
    half_type nv_tmp;
    nv_tmp = h_lhs;
    d_lhs  = static_cast<double>(nv_tmp);
    nv_tmp = h_rhs;
    d_rhs  = static_cast<double>(nv_tmp);
    if (std::is_same<view_type, ViewTypeHost>::value) {
      auto run_on_host = *this;
      run_on_host(0);
    } else {
      Kokkos::parallel_for("Test::Functor_TestHalfVolatileOperators",
                           Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this);
    }
  }

  KOKKOS_FUNCTION
  void operator()(int) const {
    volatile half_type tmp_lhs;
    half_type nv_tmp;

    // Initialze output views to catch missing test invocations
    for (int i = 0; i < N_OP_TESTS; ++i) {
      actual_lhs(i)   = 1;
      expected_lhs(i) = -1;
    }

    nv_tmp               = h_lhs;
    actual_lhs(ASSIGN)   = static_cast<double>(nv_tmp);
    expected_lhs(ASSIGN) = d_lhs;

    actual_lhs(LT)   = h_lhs < h_rhs;
    expected_lhs(LT) = d_lhs < d_rhs;

    actual_lhs(LE)   = h_lhs <= h_rhs;
    expected_lhs(LE) = d_lhs <= d_rhs;

    actual_lhs(NEQ)   = h_lhs != h_rhs;
    expected_lhs(NEQ) = d_lhs != d_rhs;

    actual_lhs(GT)   = h_lhs > h_rhs;
    expected_lhs(GT) = d_lhs > d_rhs;

    actual_lhs(GE)   = h_lhs >= h_rhs;
    expected_lhs(GE) = d_lhs >= d_rhs;

    actual_lhs(EQ)   = h_lhs == h_rhs;
    expected_lhs(EQ) = d_lhs == d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs += h_rhs;
    nv_tmp                 = tmp_lhs;
    actual_lhs(CADD_H_H)   = static_cast<double>(nv_tmp);
    expected_lhs(CADD_H_H) = d_lhs;
    expected_lhs(CADD_H_H) += d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs -= h_rhs;
    nv_tmp                 = tmp_lhs;
    actual_lhs(CSUB_H_H)   = static_cast<double>(nv_tmp);
    expected_lhs(CSUB_H_H) = d_lhs;
    expected_lhs(CSUB_H_H) -= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs *= h_rhs;
    nv_tmp                 = tmp_lhs;
    actual_lhs(CMUL_H_H)   = static_cast<double>(nv_tmp);
    expected_lhs(CMUL_H_H) = d_lhs;
    expected_lhs(CMUL_H_H) *= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs /= h_rhs;
    nv_tmp                 = tmp_lhs;
    actual_lhs(CDIV_H_H)   = static_cast<double>(nv_tmp);
    expected_lhs(CDIV_H_H) = d_lhs;
    expected_lhs(CDIV_H_H) /= d_rhs;
  }
};
#endif

template <class view_type, class half_type>
struct Functor_TestHalfOperators {
  half_type h_lhs, h_rhs;
  double d_lhs, d_rhs;
  view_type actual_lhs, expected_lhs;

  Functor_TestHalfOperators(half_type lhs = half_type(0),
                            half_type rhs = half_type(0))
      : h_lhs(lhs), h_rhs(rhs) {
    actual_lhs   = view_type("actual_lhs", N_OP_TESTS);
    expected_lhs = view_type("expected_lhs", N_OP_TESTS);
    d_lhs        = static_cast<double>(h_lhs);
    d_rhs        = static_cast<double>(h_rhs);

    if (std::is_same<view_type, ViewTypeHost>::value) {
      auto run_on_host = *this;
      run_on_host(0);
    } else {
      Kokkos::parallel_for("Test::Functor_TestHalfOperators",
                           Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this);
    }
  }

  // BEGIN: Binary Arithmetic test helpers
  template <class LhsType, class RhsType, class ExpectedResultType>
  KOKKOS_INLINE_FUNCTION void test_add(int op_test_idx,
                                       int op_test_sz_idx) const {
    auto sum = static_cast<LhsType>(h_lhs) + static_cast<RhsType>(h_rhs);
    actual_lhs(op_test_idx) = static_cast<double>(sum);

    if (std::is_same<RhsType, half_type>::value &&
        std::is_same<LhsType, half_type>::value) {
      expected_lhs(op_test_idx) = d_lhs + d_rhs;
    } else {
      if (std::is_same<LhsType, half_type>::value)
        expected_lhs(op_test_idx) = d_lhs + static_cast<RhsType>(d_rhs);
      if (std::is_same<RhsType, half_type>::value)
        expected_lhs(op_test_idx) = static_cast<LhsType>(d_lhs) + d_rhs;
    }

    actual_lhs(op_test_sz_idx)   = sizeof(sum);
    expected_lhs(op_test_sz_idx) = sizeof(ExpectedResultType);
  }

  template <class LhsType, class RhsType, class ExpectedResultType>
  KOKKOS_INLINE_FUNCTION void test_sub(int op_test_idx,
                                       int op_test_sz_idx) const {
    auto result = static_cast<LhsType>(h_lhs) - static_cast<RhsType>(h_rhs);
    actual_lhs(op_test_idx) = static_cast<double>(result);

    if (std::is_same<RhsType, half_type>::value &&
        std::is_same<LhsType, half_type>::value) {
      expected_lhs(op_test_idx) = d_lhs - d_rhs;
    } else {
      if (std::is_same<LhsType, half_type>::value)
        expected_lhs(op_test_idx) = d_lhs - static_cast<RhsType>(d_rhs);
      if (std::is_same<RhsType, half_type>::value)
        expected_lhs(op_test_idx) = static_cast<LhsType>(d_lhs) - d_rhs;
    }

    actual_lhs(op_test_sz_idx)   = sizeof(result);
    expected_lhs(op_test_sz_idx) = sizeof(ExpectedResultType);
  }

  template <class LhsType, class RhsType, class ExpectedResultType>
  KOKKOS_INLINE_FUNCTION void test_mul(int op_test_idx,
                                       int op_test_sz_idx) const {
    auto result = static_cast<LhsType>(h_lhs) * static_cast<RhsType>(h_rhs);
    actual_lhs(op_test_idx) = static_cast<double>(result);

    if (std::is_same<RhsType, half_type>::value &&
        std::is_same<LhsType, half_type>::value) {
      expected_lhs(op_test_idx) = d_lhs * d_rhs;
    } else {
      if (std::is_same<LhsType, half_type>::value)
        expected_lhs(op_test_idx) = d_lhs * static_cast<RhsType>(d_rhs);
      if (std::is_same<RhsType, half_type>::value)
        expected_lhs(op_test_idx) = static_cast<LhsType>(d_lhs) * d_rhs;
    }

    actual_lhs(op_test_sz_idx)   = sizeof(result);
    expected_lhs(op_test_sz_idx) = sizeof(ExpectedResultType);
  }

  template <class LhsType, class RhsType, class ExpectedResultType>
  KOKKOS_INLINE_FUNCTION void test_div(int op_test_idx,
                                       int op_test_sz_idx) const {
    auto result = static_cast<LhsType>(h_lhs) / static_cast<RhsType>(h_rhs);
    actual_lhs(op_test_idx) = static_cast<double>(result);

    if (std::is_same<RhsType, half_type>::value &&
        std::is_same<LhsType, half_type>::value) {
      expected_lhs(op_test_idx) = d_lhs / d_rhs;
    } else {
      if (std::is_same<LhsType, half_type>::value)
        expected_lhs(op_test_idx) = d_lhs / static_cast<RhsType>(d_rhs);
      if (std::is_same<RhsType, half_type>::value)
        expected_lhs(op_test_idx) = static_cast<LhsType>(d_lhs) / d_rhs;
    }

    actual_lhs(op_test_sz_idx)   = sizeof(result);
    expected_lhs(op_test_sz_idx) = sizeof(ExpectedResultType);
  }
  // END: Binary Arithmetic test helpers

  KOKKOS_FUNCTION
  void operator()(int) const {
    half_type tmp_lhs, tmp2_lhs, *tmp_ptr;
    double tmp_d_lhs;
    float tmp_s_lhs;
#if !defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
    using half_impl_type = typename half_type::impl_type;
#else
    using half_impl_type = half_type;
#endif  // !defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
    half_impl_type half_tmp;

    // Initialze output views to catch missing test invocations
    for (int i = 0; i < N_OP_TESTS; ++i) {
      actual_lhs(i)   = 1;
      expected_lhs(i) = -1;
    }

    tmp_lhs              = h_lhs;
    actual_lhs(ASSIGN)   = static_cast<double>(tmp_lhs);
    expected_lhs(ASSIGN) = d_lhs;

    tmp_lhs  = 0;
    tmp2_lhs = tmp_lhs           = h_lhs;
    actual_lhs(ASSIGN_CHAINED)   = static_cast<double>(tmp2_lhs);
    expected_lhs(ASSIGN_CHAINED) = d_lhs;

    actual_lhs(UNA)   = static_cast<double>(+h_lhs);
    expected_lhs(UNA) = +d_lhs;

    actual_lhs(UNS)   = static_cast<double>(-h_lhs);
    expected_lhs(UNS) = -d_lhs;

    tmp_lhs                  = h_lhs;
    tmp_d_lhs                = d_lhs;
    actual_lhs(PREFIX_INC)   = static_cast<double>(++tmp_lhs);
    expected_lhs(PREFIX_INC) = ++tmp_d_lhs;

    actual_lhs(PREFIX_DEC)   = static_cast<double>(--tmp_lhs);
    expected_lhs(PREFIX_DEC) = --tmp_d_lhs;

    // if (h_lhs != tmp_lhs) {
    //  printf("tmp_lhs = %f, h_lhs = %f\n", __half2float(tmp_lhs),
    //  __half2float(h_lhs)); Kokkos::abort("Error in half_type prefix
    //  operators");
    //}

    actual_lhs(POSTFIX_INC)   = static_cast<double>(tmp_lhs++);
    expected_lhs(POSTFIX_INC) = tmp_d_lhs++;

    actual_lhs(POSTFIX_DEC)   = static_cast<double>(tmp_lhs--);
    expected_lhs(POSTFIX_DEC) = tmp_d_lhs--;

    // if (h_lhs != tmp_lhs) {
    //  printf("tmp_lhs = %f, h_lhs = %f\n", __half2float(tmp_lhs),
    //  __half2float(h_lhs)); Kokkos::abort("Error in half_type postfix
    //  operators");
    //}

    tmp_lhs = h_lhs;
    tmp_lhs += h_rhs;
    actual_lhs(CADD_H_H)   = static_cast<double>(tmp_lhs);
    expected_lhs(CADD_H_H) = d_lhs;
    expected_lhs(CADD_H_H) += d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs += static_cast<float>(d_rhs);
    actual_lhs(CADD_H_S)   = static_cast<double>(tmp_lhs);
    expected_lhs(CADD_H_S) = d_lhs;
    expected_lhs(CADD_H_S) += d_rhs;

    tmp_s_lhs = static_cast<float>(h_lhs);
    tmp_s_lhs += h_rhs;
    actual_lhs(CADD_S_H)   = static_cast<double>(tmp_s_lhs);
    expected_lhs(CADD_S_H) = d_lhs;
    expected_lhs(CADD_S_H) += d_rhs;

    tmp_lhs = static_cast<double>(h_lhs);
    tmp_lhs += static_cast<double>(d_rhs);
    actual_lhs(CADD_H_D)   = static_cast<double>(tmp_lhs);
    expected_lhs(CADD_H_D) = d_lhs;
    expected_lhs(CADD_H_D) += d_rhs;

    tmp_d_lhs = static_cast<double>(h_lhs);
    tmp_d_lhs += h_rhs;
    actual_lhs(CADD_D_H)   = static_cast<double>(tmp_d_lhs);
    expected_lhs(CADD_D_H) = d_lhs;
    expected_lhs(CADD_D_H) += d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs -= h_rhs;
    actual_lhs(CSUB_H_H)   = static_cast<double>(tmp_lhs);
    expected_lhs(CSUB_H_H) = d_lhs;
    expected_lhs(CSUB_H_H) -= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs -= static_cast<float>(d_rhs);
    actual_lhs(CSUB_H_S)   = static_cast<double>(tmp_lhs);
    expected_lhs(CSUB_H_S) = d_lhs;
    expected_lhs(CSUB_H_S) -= d_rhs;

    tmp_s_lhs = static_cast<float>(h_lhs);
    tmp_s_lhs -= h_rhs;
    actual_lhs(CSUB_S_H)   = static_cast<double>(tmp_s_lhs);
    expected_lhs(CSUB_S_H) = d_lhs;
    expected_lhs(CSUB_S_H) -= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs -= d_rhs;
    actual_lhs(CSUB_H_D)   = static_cast<double>(tmp_lhs);
    expected_lhs(CSUB_H_D) = d_lhs;
    expected_lhs(CSUB_H_D) -= d_rhs;

    tmp_d_lhs = static_cast<double>(h_lhs);
    tmp_d_lhs -= h_rhs;
    actual_lhs(CSUB_D_H)   = tmp_d_lhs;
    expected_lhs(CSUB_D_H) = d_lhs;
    expected_lhs(CSUB_D_H) -= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs *= h_rhs;
    actual_lhs(CMUL_H_H)   = static_cast<double>(tmp_lhs);
    expected_lhs(CMUL_H_H) = d_lhs;
    expected_lhs(CMUL_H_H) *= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs *= static_cast<float>(d_rhs);
    actual_lhs(CMUL_H_S)   = static_cast<double>(tmp_lhs);
    expected_lhs(CMUL_H_S) = d_lhs;
    expected_lhs(CMUL_H_S) *= d_rhs;

    tmp_s_lhs = static_cast<float>(h_lhs);
    tmp_s_lhs *= h_rhs;
    actual_lhs(CMUL_S_H)   = static_cast<double>(tmp_s_lhs);
    expected_lhs(CMUL_S_H) = d_lhs;
    expected_lhs(CMUL_S_H) *= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs *= d_rhs;
    actual_lhs(CMUL_H_D)   = static_cast<double>(tmp_lhs);
    expected_lhs(CMUL_H_D) = d_lhs;
    expected_lhs(CMUL_H_D) *= d_rhs;

    tmp_d_lhs = static_cast<double>(h_lhs);
    tmp_d_lhs *= h_rhs;
    actual_lhs(CMUL_D_H)   = tmp_d_lhs;
    expected_lhs(CMUL_D_H) = d_lhs;
    expected_lhs(CMUL_D_H) *= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs /= h_rhs;
    actual_lhs(CDIV_H_H)   = static_cast<double>(tmp_lhs);
    expected_lhs(CDIV_H_H) = d_lhs;
    expected_lhs(CDIV_H_H) /= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs /= static_cast<float>(d_rhs);
    actual_lhs(CDIV_H_S)   = static_cast<double>(tmp_lhs);
    expected_lhs(CDIV_H_S) = d_lhs;
    expected_lhs(CDIV_H_S) /= d_rhs;

    tmp_s_lhs = static_cast<float>(h_lhs);
    tmp_s_lhs /= h_rhs;
    actual_lhs(CDIV_S_H)   = static_cast<double>(tmp_s_lhs);
    expected_lhs(CDIV_S_H) = d_lhs;
    expected_lhs(CDIV_S_H) /= d_rhs;

    tmp_lhs = h_lhs;
    tmp_lhs /= d_rhs;
    actual_lhs(CDIV_H_D)   = static_cast<double>(tmp_lhs);
    expected_lhs(CDIV_H_D) = d_lhs;
    expected_lhs(CDIV_H_D) /= d_rhs;

    tmp_d_lhs = static_cast<double>(h_lhs);
    tmp_d_lhs /= h_rhs;
    actual_lhs(CDIV_D_H)   = tmp_d_lhs;
    expected_lhs(CDIV_D_H) = d_lhs;
    expected_lhs(CDIV_D_H) /= d_rhs;

    test_add<half_type, half_type, half_type>(ADD_H_H, ADD_H_H_SZ);
    test_add<float, half_type, float>(ADD_S_H, ADD_S_H_SZ);
    test_add<double, half_type, double>(ADD_D_H, ADD_D_H_SZ);
    test_add<short int, half_type, half_type>(ADD_SI_H, ADD_SI_H_SZ);
    test_add<int, half_type, half_type>(ADD_I_H, ADD_I_H_SZ);
    test_add<long int, half_type, half_type>(ADD_LI_H, ADD_LI_H_SZ);
    test_add<long long int, half_type, half_type>(ADD_LLI_H, ADD_LLI_H_SZ);
    test_add<half_type, float, float>(ADD_H_S, ADD_H_S_SZ);
    test_add<half_type, double, double>(ADD_H_D, ADD_H_D_SZ);
    test_add<half_type, short int, half_type>(ADD_H_SI, ADD_H_SI_SZ);
    test_add<half_type, int, half_type>(ADD_H_I, ADD_H_I_SZ);
    test_add<half_type, long int, half_type>(ADD_H_LI, ADD_H_LI_SZ);
    test_add<half_type, long long int, half_type>(ADD_H_LLI, ADD_H_LLI_SZ);

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_lhs >= 0) {
      test_add<unsigned short int, half_type, half_type>(ADD_USI_H,
                                                         ADD_USI_H_SZ);
      test_add<unsigned int, half_type, half_type>(ADD_UI_H, ADD_UI_H_SZ);
      test_add<unsigned long int, half_type, half_type>(ADD_ULI_H,
                                                        ADD_ULI_H_SZ);
      test_add<unsigned long long int, half_type, half_type>(ADD_ULLI_H,
                                                             ADD_ULLI_H_SZ);
    } else {
      actual_lhs(ADD_USI_H)     = expected_lhs(ADD_USI_H);
      actual_lhs(ADD_USI_H_SZ)  = expected_lhs(ADD_USI_H_SZ);
      actual_lhs(ADD_UI_H)      = expected_lhs(ADD_UI_H);
      actual_lhs(ADD_UI_H_SZ)   = expected_lhs(ADD_UI_H_SZ);
      actual_lhs(ADD_ULI_H)     = expected_lhs(ADD_ULI_H);
      actual_lhs(ADD_ULI_H_SZ)  = expected_lhs(ADD_ULI_H_SZ);
      actual_lhs(ADD_ULLI_H)    = expected_lhs(ADD_ULLI_H);
      actual_lhs(ADD_ULLI_H_SZ) = expected_lhs(ADD_ULLI_H_SZ);
    }

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_rhs >= 0) {
      test_add<half_type, unsigned short int, half_type>(ADD_H_USI,
                                                         ADD_H_USI_SZ);
      test_add<half_type, unsigned int, half_type>(ADD_H_UI, ADD_H_UI_SZ);
      test_add<half_type, unsigned long int, half_type>(ADD_H_ULI,
                                                        ADD_H_ULI_SZ);
      test_add<half_type, unsigned long long int, half_type>(ADD_H_ULLI,
                                                             ADD_H_ULLI_SZ);
    } else {
      actual_lhs(ADD_H_USI)     = expected_lhs(ADD_H_USI);
      actual_lhs(ADD_H_USI_SZ)  = expected_lhs(ADD_H_USI_SZ);
      actual_lhs(ADD_H_UI)      = expected_lhs(ADD_H_UI);
      actual_lhs(ADD_H_UI_SZ)   = expected_lhs(ADD_H_UI_SZ);
      actual_lhs(ADD_H_ULI)     = expected_lhs(ADD_H_ULI);
      actual_lhs(ADD_H_ULI_SZ)  = expected_lhs(ADD_H_ULI_SZ);
      actual_lhs(ADD_H_ULLI)    = expected_lhs(ADD_H_ULLI);
      actual_lhs(ADD_H_ULLI_SZ) = expected_lhs(ADD_H_ULLI_SZ);
    }

    test_sub<half_type, half_type, half_type>(SUB_H_H, SUB_H_H_SZ);
    test_sub<float, half_type, float>(SUB_S_H, SUB_S_H_SZ);
    test_sub<double, half_type, double>(SUB_D_H, SUB_D_H_SZ);
    test_sub<short int, half_type, half_type>(SUB_SI_H, SUB_SI_H_SZ);
    test_sub<int, half_type, half_type>(SUB_I_H, SUB_I_H_SZ);
    test_sub<long int, half_type, half_type>(SUB_LI_H, SUB_LI_H_SZ);
    test_sub<long long int, half_type, half_type>(SUB_LLI_H, SUB_LLI_H_SZ);
    test_sub<half_type, float, float>(SUB_H_S, SUB_H_S_SZ);
    test_sub<half_type, double, double>(SUB_H_D, SUB_H_D_SZ);
    test_sub<half_type, short int, half_type>(SUB_H_SI, SUB_H_SI_SZ);
    test_sub<half_type, int, half_type>(SUB_H_I, SUB_H_I_SZ);
    test_sub<half_type, long int, half_type>(SUB_H_LI, SUB_H_LI_SZ);
    test_sub<half_type, long long int, half_type>(SUB_H_LLI, SUB_H_LLI_SZ);

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_lhs >= half_type(0)) {
      test_sub<unsigned short int, half_type, half_type>(SUB_USI_H,
                                                         SUB_USI_H_SZ);
      test_sub<unsigned int, half_type, half_type>(SUB_UI_H, SUB_UI_H_SZ);
      test_sub<unsigned long int, half_type, half_type>(SUB_ULI_H,
                                                        SUB_ULI_H_SZ);
      test_sub<unsigned long long int, half_type, half_type>(SUB_ULLI_H,
                                                             SUB_ULLI_H_SZ);
    } else {
      actual_lhs(SUB_USI_H)     = expected_lhs(SUB_USI_H);
      actual_lhs(SUB_USI_H_SZ)  = expected_lhs(SUB_USI_H_SZ);
      actual_lhs(SUB_UI_H)      = expected_lhs(SUB_UI_H);
      actual_lhs(SUB_UI_H_SZ)   = expected_lhs(SUB_UI_H_SZ);
      actual_lhs(SUB_ULI_H)     = expected_lhs(SUB_ULI_H);
      actual_lhs(SUB_ULI_H_SZ)  = expected_lhs(SUB_ULI_H_SZ);
      actual_lhs(SUB_ULLI_H)    = expected_lhs(SUB_ULLI_H);
      actual_lhs(SUB_ULLI_H_SZ) = expected_lhs(SUB_ULLI_H_SZ);
    }

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_rhs >= half_type(0)) {
      test_sub<half_type, unsigned short int, half_type>(SUB_H_USI,
                                                         SUB_H_USI_SZ);
      test_sub<half_type, unsigned int, half_type>(SUB_H_UI, SUB_H_UI_SZ);
      test_sub<half_type, unsigned long int, half_type>(SUB_H_ULI,
                                                        SUB_H_ULI_SZ);
      test_sub<half_type, unsigned long long int, half_type>(SUB_H_ULLI,
                                                             SUB_H_ULLI_SZ);
    } else {
      actual_lhs(SUB_H_USI)     = expected_lhs(SUB_H_USI);
      actual_lhs(SUB_H_USI_SZ)  = expected_lhs(SUB_H_USI_SZ);
      actual_lhs(SUB_H_UI)      = expected_lhs(SUB_H_UI);
      actual_lhs(SUB_H_UI_SZ)   = expected_lhs(SUB_H_UI_SZ);
      actual_lhs(SUB_H_ULI)     = expected_lhs(SUB_H_ULI);
      actual_lhs(SUB_H_ULI_SZ)  = expected_lhs(SUB_H_ULI_SZ);
      actual_lhs(SUB_H_ULLI)    = expected_lhs(SUB_H_ULLI);
      actual_lhs(SUB_H_ULLI_SZ) = expected_lhs(SUB_H_ULLI_SZ);
    }

    test_mul<half_type, half_type, half_type>(MUL_H_H, MUL_H_H_SZ);
    test_mul<float, half_type, float>(MUL_S_H, MUL_S_H_SZ);
    test_mul<double, half_type, double>(MUL_D_H, MUL_D_H_SZ);
    test_mul<short int, half_type, half_type>(MUL_SI_H, MUL_SI_H_SZ);
    test_mul<int, half_type, half_type>(MUL_I_H, MUL_I_H_SZ);
    test_mul<long int, half_type, half_type>(MUL_LI_H, MUL_LI_H_SZ);
    test_mul<long long int, half_type, half_type>(MUL_LLI_H, MUL_LLI_H_SZ);
    test_mul<half_type, float, float>(MUL_H_S, MUL_H_S_SZ);
    test_mul<half_type, double, double>(MUL_H_D, MUL_H_D_SZ);
    test_mul<half_type, short int, half_type>(MUL_H_SI, MUL_H_SI_SZ);
    test_mul<half_type, int, half_type>(MUL_H_I, MUL_H_I_SZ);
    test_mul<half_type, long int, half_type>(MUL_H_LI, MUL_H_LI_SZ);
    test_mul<half_type, long long int, half_type>(MUL_H_LLI, MUL_H_LLI_SZ);

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_lhs >= half_type(0)) {
      test_mul<unsigned short int, half_type, half_type>(MUL_USI_H,
                                                         MUL_USI_H_SZ);
      test_mul<unsigned int, half_type, half_type>(MUL_UI_H, MUL_UI_H_SZ);
      test_mul<unsigned long int, half_type, half_type>(MUL_ULI_H,
                                                        MUL_ULI_H_SZ);
      test_mul<unsigned long long int, half_type, half_type>(MUL_ULLI_H,
                                                             MUL_ULLI_H_SZ);
    } else {
      actual_lhs(MUL_USI_H)     = expected_lhs(MUL_USI_H);
      actual_lhs(MUL_UI_H)      = expected_lhs(MUL_UI_H);
      actual_lhs(MUL_ULI_H)     = expected_lhs(MUL_ULI_H);
      actual_lhs(MUL_ULLI_H)    = expected_lhs(MUL_ULLI_H);
      actual_lhs(MUL_USI_H_SZ)  = expected_lhs(MUL_USI_H_SZ);
      actual_lhs(MUL_UI_H_SZ)   = expected_lhs(MUL_UI_H_SZ);
      actual_lhs(MUL_ULI_H_SZ)  = expected_lhs(MUL_ULI_H_SZ);
      actual_lhs(MUL_ULLI_H_SZ) = expected_lhs(MUL_ULLI_H_SZ);
    }

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_rhs >= half_type(0)) {
      test_mul<half_type, unsigned short int, half_type>(MUL_H_USI,
                                                         MUL_H_USI_SZ);
      test_mul<half_type, unsigned int, half_type>(MUL_H_UI, MUL_H_UI_SZ);
      test_mul<half_type, unsigned long int, half_type>(MUL_H_ULI,
                                                        MUL_H_ULI_SZ);
      test_mul<half_type, unsigned long long int, half_type>(MUL_H_ULLI,
                                                             MUL_H_ULLI_SZ);
    } else {
      actual_lhs(MUL_H_USI)     = expected_lhs(MUL_H_USI);
      actual_lhs(MUL_H_UI)      = expected_lhs(MUL_H_UI);
      actual_lhs(MUL_H_ULI)     = expected_lhs(MUL_H_ULI);
      actual_lhs(MUL_H_ULLI)    = expected_lhs(MUL_H_ULLI);
      actual_lhs(MUL_H_USI_SZ)  = expected_lhs(MUL_H_USI_SZ);
      actual_lhs(MUL_H_UI_SZ)   = expected_lhs(MUL_H_UI_SZ);
      actual_lhs(MUL_H_ULI_SZ)  = expected_lhs(MUL_H_ULI_SZ);
      actual_lhs(MUL_H_ULLI_SZ) = expected_lhs(MUL_H_ULLI_SZ);
    }

    test_div<half_type, half_type, half_type>(DIV_H_H, DIV_H_H_SZ);
    test_div<float, half_type, float>(DIV_S_H, DIV_S_H_SZ);
    test_div<double, half_type, double>(DIV_D_H, DIV_D_H_SZ);
    test_div<short int, half_type, half_type>(DIV_SI_H, DIV_SI_H_SZ);
    test_div<int, half_type, half_type>(DIV_I_H, DIV_I_H_SZ);
    test_div<long int, half_type, half_type>(DIV_LI_H, DIV_LI_H_SZ);
    test_div<long long int, half_type, half_type>(DIV_LLI_H, DIV_LLI_H_SZ);
    test_div<half_type, float, float>(DIV_H_S, DIV_H_S_SZ);
    test_div<half_type, double, double>(DIV_H_D, DIV_H_D_SZ);

    // Check for division by zero due to truncation by half_type -> integral
    // cast
    if (h_rhs >= half_type(1) || h_rhs <= half_type(-1)) {
      test_div<half_type, short int, half_type>(DIV_H_SI, DIV_H_SI_SZ);
      test_div<half_type, int, half_type>(DIV_H_I, DIV_H_I_SZ);
      test_div<half_type, long int, half_type>(DIV_H_LI, DIV_H_LI_SZ);
      test_div<half_type, long long int, half_type>(DIV_H_LLI, DIV_H_LLI_SZ);
    } else {
      actual_lhs(DIV_H_SI)     = expected_lhs(DIV_H_SI);
      actual_lhs(DIV_H_I)      = expected_lhs(DIV_H_I);
      actual_lhs(DIV_H_LI)     = expected_lhs(DIV_H_LI);
      actual_lhs(DIV_H_LLI)    = expected_lhs(DIV_H_LLI);
      actual_lhs(DIV_H_SI_SZ)  = expected_lhs(DIV_H_SI_SZ);
      actual_lhs(DIV_H_I_SZ)   = expected_lhs(DIV_H_I_SZ);
      actual_lhs(DIV_H_LI_SZ)  = expected_lhs(DIV_H_LI_SZ);
      actual_lhs(DIV_H_LLI_SZ) = expected_lhs(DIV_H_LLI_SZ);
    }

    // Check for potential overflow due to negative half_type -> unsigned
    // integral cast
    if (h_lhs >= half_type(0)) {
      test_div<unsigned short int, half_type, half_type>(DIV_USI_H,
                                                         DIV_USI_H_SZ);
      test_div<unsigned int, half_type, half_type>(DIV_UI_H, DIV_UI_H_SZ);
      test_div<unsigned long int, half_type, half_type>(DIV_ULI_H,
                                                        DIV_ULI_H_SZ);
      test_div<unsigned long long int, half_type, half_type>(DIV_ULLI_H,
                                                             DIV_ULLI_H_SZ);
    } else {
      actual_lhs(DIV_USI_H)     = expected_lhs(DIV_USI_H);
      actual_lhs(DIV_UI_H)      = expected_lhs(DIV_UI_H);
      actual_lhs(DIV_ULI_H)     = expected_lhs(DIV_ULI_H);
      actual_lhs(DIV_ULLI_H)    = expected_lhs(DIV_ULLI_H);
      actual_lhs(DIV_USI_H_SZ)  = expected_lhs(DIV_USI_H_SZ);
      actual_lhs(DIV_UI_H_SZ)   = expected_lhs(DIV_UI_H_SZ);
      actual_lhs(DIV_ULI_H_SZ)  = expected_lhs(DIV_ULI_H_SZ);
      actual_lhs(DIV_ULLI_H_SZ) = expected_lhs(DIV_ULLI_H_SZ);
    }

    // Check for division by zero due to truncation by half_type -> integral
    // cast
    if (h_rhs >= half_type(1)) {
      test_div<half_type, unsigned short int, half_type>(DIV_H_USI,
                                                         DIV_H_USI_SZ);
      test_div<half_type, unsigned int, half_type>(DIV_H_UI, DIV_H_UI_SZ);
      test_div<half_type, unsigned long int, half_type>(DIV_H_ULI,
                                                        DIV_H_ULI_SZ);
      test_div<half_type, unsigned long long int, half_type>(DIV_H_ULLI,
                                                             DIV_H_ULLI_SZ);
    } else {
      actual_lhs(DIV_H_USI)     = expected_lhs(DIV_H_USI);
      actual_lhs(DIV_H_USI_SZ)  = expected_lhs(DIV_H_USI_SZ);
      actual_lhs(DIV_H_UI)      = expected_lhs(DIV_H_UI);
      actual_lhs(DIV_H_UI_SZ)   = expected_lhs(DIV_H_UI_SZ);
      actual_lhs(DIV_H_ULI)     = expected_lhs(DIV_H_ULI);
      actual_lhs(DIV_H_ULI_SZ)  = expected_lhs(DIV_H_ULI_SZ);
      actual_lhs(DIV_H_ULLI)    = expected_lhs(DIV_H_ULLI);
      actual_lhs(DIV_H_ULLI_SZ) = expected_lhs(DIV_H_ULLI_SZ);
    }

    // TODO: figure out why operator{!,&&,||} are returning __nv_bool
    actual_lhs(NEG)   = static_cast<double>(!h_lhs);
    expected_lhs(NEG) = !d_lhs;

    actual_lhs(AND)   = static_cast<double>(half_type(0) && h_lhs);
    expected_lhs(AND) = double(0) && d_lhs;

    actual_lhs(OR)   = static_cast<double>(h_lhs || half_type(1));
    expected_lhs(OR) = d_lhs || double(1);

    actual_lhs(EQ)   = h_lhs == h_rhs;
    expected_lhs(EQ) = d_lhs == d_rhs;

    actual_lhs(NEQ)   = h_lhs != h_rhs;
    expected_lhs(NEQ) = d_lhs != d_rhs;

    actual_lhs(LT)   = h_lhs < h_rhs;
    expected_lhs(LT) = d_lhs < d_rhs;

    actual_lhs(GT)   = h_lhs > h_rhs;
    expected_lhs(GT) = d_lhs > d_rhs;

    actual_lhs(LE)   = h_lhs <= h_rhs;
    expected_lhs(LE) = d_lhs <= d_rhs;

    actual_lhs(GE)   = h_lhs >= h_rhs;
    expected_lhs(GE) = d_lhs >= d_rhs;

    // actual_lhs(TW)   = h_lhs <=> h_rhs;  // Need C++20?
    // expected_lhs(TW) = d_lhs <=> d_rhs;  // Need C++20?

    actual_lhs(PASS_BY_REF) = static_cast<double>(accept_ref(h_lhs));

    // Use accept_ref and accept_ref_expected to ensure the compiler
    // does not optimize out the casts half_type -> double -> half_type.
    // Note that these casts are accompanied by rounding. For the bhalf_t
    // epsilon, these rounding policies used for casting is enough to cause
    // the unit tests to fail.
    // In short, one cannot simply assign static_cast<double>(h_lhs) to
    // expected_lhs(PASS_BY_REF).
    expected_lhs(PASS_BY_REF) = accept_ref_expected(h_lhs);

    half_tmp = static_cast<float>(h_lhs);
    tmp_ptr  = &(tmp_lhs = half_tmp);
    if (tmp_ptr != &tmp_lhs)
      Kokkos::abort("Error in half_type address-of operator");
    actual_lhs(AO_IMPL_HALF)   = static_cast<double>(*tmp_ptr);
    expected_lhs(AO_IMPL_HALF) = d_lhs;

    tmp2_lhs = h_lhs;
    tmp_ptr  = &(tmp_lhs = tmp2_lhs);
    if (tmp_ptr != &tmp_lhs)
      Kokkos::abort("Error in half_type address-of operator");
    actual_lhs(AO_HALF_T)   = static_cast<double>(tmp_ptr[0]);
    expected_lhs(AO_HALF_T) = d_lhs;

    // TODO: Check upcasting and downcasting in large expressions involving
    // integral and floating point types
  }
};

template <class half_type>
void __test_half_operators(half_type h_lhs, half_type h_rhs) {
  double epsilon = Kokkos::Experimental::epsilon<half_type>::value;

  Functor_TestHalfOperators<ViewType, half_type> f_device(h_lhs, h_rhs);
  Functor_TestHalfOperators<ViewTypeHost, half_type> f_host(h_lhs, h_rhs);
  typename ViewType::HostMirror f_device_actual_lhs =
      Kokkos::create_mirror_view(f_device.actual_lhs);
  typename ViewType::HostMirror f_device_expected_lhs =
      Kokkos::create_mirror_view(f_device.expected_lhs);

  ExecutionSpace().fence();
  Kokkos::deep_copy(f_device_actual_lhs, f_device.actual_lhs);
  Kokkos::deep_copy(f_device_expected_lhs, f_device.expected_lhs);
  for (int op_test = 0; op_test < N_OP_TESTS; op_test++) {
    // printf("op_test = %d\n", op_test);
    ASSERT_NEAR(f_device_actual_lhs(op_test), f_device_expected_lhs(op_test),
                epsilon);
    ASSERT_NEAR(f_host.actual_lhs(op_test), f_host.expected_lhs(op_test),
                epsilon);
  }

// volatile-qualified parameter type 'volatile half_type' is deprecated
#if !defined(KOKKOS_ENABLE_CXX20) && !defined(KOKKOS_ENABLE_CXX23)
  // Test partial volatile support
  volatile half_type _h_lhs = h_lhs;
  volatile half_type _h_rhs = h_rhs;
  Functor_TestHalfVolatileOperators<ViewType, half_type> f_volatile_device(
      _h_lhs, _h_rhs);
  Functor_TestHalfVolatileOperators<ViewTypeHost, half_type> f_volatile_host(
      _h_lhs, _h_rhs);

  ExecutionSpace().fence();
  Kokkos::deep_copy(f_device_actual_lhs, f_device.actual_lhs);
  Kokkos::deep_copy(f_device_expected_lhs, f_device.expected_lhs);
  for (int op_test = 0; op_test < N_OP_TESTS; op_test++) {
    // printf("op_test = %d\n", op_test);
    if (op_test == ASSIGN || op_test == LT || op_test == LE || op_test == NEQ ||
        op_test == EQ || op_test == GT || op_test == GE ||
        op_test == CADD_H_H || op_test == CSUB_H_H || op_test == CMUL_H_H ||
        op_test == CDIV_H_H) {
      ASSERT_NEAR(f_device_actual_lhs(op_test), f_device_expected_lhs(op_test),
                  epsilon);
      ASSERT_NEAR(f_host.actual_lhs(op_test), f_host.expected_lhs(op_test),
                  epsilon);
    }
  }
#endif

  // is_trivially_copyable is false with the addition of explicit
  // copy constructors that are required for supporting reductions
  // ASSERT_TRUE(std::is_trivially_copyable<half_type>::value);

  constexpr size_t n       = 2;
  constexpr size_t n_bytes = sizeof(half_type) * n;
  const half_type h_arr0 = half_type(0x89ab), h_arr1 = half_type(0xcdef);
  half_type h_arr[n];
  char c_arr[n_bytes], *h_arr_ptr = nullptr;
  size_t i;

  h_arr[0]  = h_arr0;
  h_arr[1]  = h_arr1;
  h_arr_ptr = reinterpret_cast<char*>(h_arr);

  std::memcpy(c_arr, h_arr, n_bytes);
  for (i = 0; i < n_bytes; i++) ASSERT_EQ(c_arr[i], h_arr_ptr[i]);

  ASSERT_EQ(h_arr[0], h_arr0);
  ASSERT_EQ(h_arr[1], h_arr1);
}

void test_half_operators() {
  half_t h_lhs = half_t(0.23458), h_rhs = half_t(0.67898);
  for (int i = -3; i < 2; i++) {
    // printf("%f OP %f\n", float(h_lhs + cast_to_half(i + 1)), float(h_rhs +
    // cast_to_half(i)));
    __test_half_operators<half_t>(h_lhs + cast_to_half(i + 1),
                                  h_rhs + cast_to_half(i));
    // TODO: __test_half_operators(h_lhs + cast_to_half(i + 1), half_t(0));
    // TODO: __test_half_operators(half_t(0), h_rhs + cast_to_half(i));
  }
}

void test_bhalf_operators() {
  bhalf_t h_lhs = bhalf_t(0.23458), h_rhs = bhalf_t(0.67898);
  for (int i = -2; i < 2; i++) {
    // printf("%f OP %f\n", float(h_lhs + cast_to_bhalf(i + 1)), float(h_rhs +
    // cast_to_bhalf(i)));
    __test_half_operators<bhalf_t>(h_lhs + cast_to_bhalf(i + 1),
                                   h_rhs + cast_to_bhalf(i));
  }
}

TEST(TEST_CATEGORY, half_operators) { test_half_operators(); }
TEST(TEST_CATEGORY, bhalf_operators) { test_bhalf_operators(); }
}  // namespace Test
#endif  // TESTHALFOPERATOR_HPP_
