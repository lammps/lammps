#pragma once
#ifndef __GEMM_HPP__
#define __GEMM_HPP__

/// \file gemm.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho {

  using namespace std;

  template<int ArgTransA, int ArgTransB, int ArgAlgo,
           int ArgVariant = Variant::One,
           template<int,int> class ControlType = Control>
  struct Gemm {

    // data-parallel interface
    // =======================
    template<typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename ExecViewTypeA::policy_type &policy,
                      const typename ExecViewTypeA::policy_type::member_type &member,
                      const ScalarType alpha,
                      typename ExecViewTypeA::matrix_type &A,
                      typename ExecViewTypeB::matrix_type &B,
                      const ScalarType beta,
                      typename ExecViewTypeC::matrix_type &C);

    // task-data parallel interface
    // ============================
    template<typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeB,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename ExecViewTypeA::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      typename ExecViewTypeA::matrix_type _A;
      typename ExecViewTypeB::matrix_type _B;
      typename ExecViewTypeC::matrix_type _C;

      policy_type _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const policy_type & P,
                  const ScalarType alpha,
                  const typename ExecViewTypeA::matrix_type & A,
                  const typename ExecViewTypeB::matrix_type & B,
                  const ScalarType beta,
                  const typename ExecViewTypeC::matrix_type & C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _B(B),
          _C(C),
          _policy(P)
      { }

      string Label() const { return "Gemm"; }

      // task execution
      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Gemm::invoke<ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>(_policy, _policy.member_single(),
                             _alpha, _A, _B, _beta, _C);
      }

      // task-data execution
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        r_val = Gemm::invoke<ScalarType,ExecViewTypeA,ExecViewTypeB,ExecViewTypeC>(_policy, member,
                             _alpha, _A, _B, _beta, _C);
      }

    };

  };

}


// #include "gemm_nt_nt.hpp"
#include "gemm_ct_nt.hpp"

#endif
