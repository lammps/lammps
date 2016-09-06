#pragma once
#ifndef __CHOL_HPP__
#define __CHOL_HPP__

/// \file chol.hpp
/// \brief Incomplete Cholesky factorization front interface.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "control.hpp"
#include "partition.hpp"

namespace Tacho { 

  using namespace std;

  // tasking interface
  // * default behavior is for non-by-blocks tasks
  // * control is only used for by-blocks algorithms
  // ===============================================
  template<int ArgUplo, int ArgAlgo, 
           int ArgVariant = Variant::One,                  
           template<int,int> class ControlType = Control>  
  class Chol {
  public:
    
    // function interface
    // ==================
    template<typename ExecViewType>
    KOKKOS_INLINE_FUNCTION
    static int invoke(typename ExecViewType::policy_type &policy, 
                      const typename ExecViewType::policy_type::member_type &member, 
                      typename ExecViewType::matrix_type &A);

    // task-data parallel interface
    // ============================
    template<typename ExecViewType>
    class TaskFunctor {
    public:
      typedef typename ExecViewType::policy_type policy_type;
      typedef typename policy_type::member_type member_type;
      typedef int value_type;
      
    private:
      typename ExecViewType::matrix_type _A;
      
      policy_type _policy;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor(const policy_type & P ,
                  const typename ExecViewType::matrix_type & A)
        : _A(A),
          _policy(P)
      { } 
      
      string Label() const { return "Chol"; }
      
      // task execution
      KOKKOS_INLINE_FUNCTION
      void apply(value_type &r_val) {
        r_val = Chol::invoke<ExecViewType>(_policy, _policy.member_single(), _A);
      }

      // task-data execution
      KOKKOS_INLINE_FUNCTION
      void apply(const member_type &member, value_type &r_val) {

        const int result = Chol::invoke<ExecViewType>(_policy, member, _A);

        if ( 0 == member.team_rank() ) { r_val = result ; }

      }

    };

  };
}


// unblocked version blas operations
#include "scale.hpp"

// blocked version blas operations
#include "gemm.hpp"
#include "trsm.hpp"
#include "herk.hpp"

// cholesky
#include "chol_u.hpp"

#endif
