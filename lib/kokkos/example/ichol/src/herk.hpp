#pragma once
#ifndef __HERK_HPP__
#define __HERK_HPP__

/// \file herk.hpp
/// \brief Sparse hermitian rank one update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho {

  using namespace std;

  template<int ArgUplo, int ArgTrans, int ArgAlgo,
           int ArgVariant = Variant::One,
           template<int,int> class ControlType = Control>
  struct Herk {

    // data-parallel interface
    // =======================
    template<typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename ExecViewTypeA::policy_type &policy,
                      const typename ExecViewTypeA::policy_type::member_type &member,
                      const ScalarType alpha,
                      typename ExecViewTypeA::matrix_type &A,
                      const ScalarType beta,
                      typename ExecViewTypeC::matrix_type &C);

    // task-data parallel interface
    // ============================
    template<typename ScalarType,
             typename ExecViewTypeA,
             typename ExecViewTypeC>
    class TaskFunctor {
    public:
      typedef typename ExecViewTypeA::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;

    private:
      ScalarType _alpha, _beta;
      typename ExecViewTypeA::matrix_type _A;
      typename ExecViewTypeC::matrix_type _C;

      policy_type _policy;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const policy_type & P,
                  const ScalarType alpha,
                  const typename ExecViewTypeA::matrix_type & A,
                  const ScalarType beta,
                  const typename ExecViewTypeC::matrix_type & C)
        : _alpha(alpha),
          _beta(beta),
          _A(A),
          _C(C),
          _policy(P)
      { }

      string Label() const { return "Herk"; }

      // task execution
      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Herk::invoke<ScalarType,ExecViewTypeA,ExecViewTypeC>(_policy, _policy.member_single(), 
                             _alpha, _A, _beta, _C);
      }

      // task-data execution
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {
        r_val = Herk::invoke<ScalarType,ExecViewTypeA,ExecViewTypeC>(_policy, member, 
                             _alpha, _A, _beta, _C);
      }

    };

  };

}

#include "herk_u_ct.hpp"

#endif
