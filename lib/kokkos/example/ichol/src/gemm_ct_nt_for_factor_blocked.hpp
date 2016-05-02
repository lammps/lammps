#pragma once
#ifndef __GEMM_CT_NT_FOR_FACTOR_BLOCKED_HPP__
#define __GEMM_CT_NT_FOR_FACTOR_BLOCKED_HPP__

/// \file gemm_ct_nt_for_factor_blocked.hpp
/// \brief Sparse matrix-matrix multiplication on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  using namespace std;

  // Gemm used in the factorization phase
  // ====================================
  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB,
           typename CrsExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::ConjTranspose,Trans::NoTranspose,
       AlgoGemm::ForFactorBlocked>
  ::invoke(typename CrsExecViewTypeA::policy_type &policy,
           const typename CrsExecViewTypeA::policy_type::member_type &member,
           const ScalarType alpha,
           typename CrsExecViewTypeA::matrix_type &A,
           typename CrsExecViewTypeB::matrix_type &B,
           const ScalarType beta,
           typename CrsExecViewTypeC::matrix_type &C) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;


if ( false && member.team_rank() == 0 ) {
 printf("Gemm [%d +%d)x[%d +%d)\n"
       , C.OffsetRows()
       , C.NumRows()
       , C.OffsetCols()
       , C.NumCols()
       );
}

    // scale the matrix C with beta
    scaleCrsMatrix<ScalarType,CrsExecViewTypeC>(member, beta, C);

    // Sparse matrix-matrix multiply:
    // C(i,j) += alpha*A'(i,k)*B(k,j)

    const ordinal_type mA = A.NumRows();
    for (ordinal_type k=0;k<mA;++k) {
      row_view_type &a = A.RowView(k);
      const ordinal_type nnz_a = a.NumNonZeros();

      row_view_type &b = B.RowView(k);
      const ordinal_type nnz_b = b.NumNonZeros();

      if (nnz_a > 0 && nnz_b > 0 ) {
#if 0
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, nnz_a),
          [&](const ordinal_type i) {
             const ordinal_type row_at_i  = a.Col(i);
             const value_type   val_at_ik = a.Value(i);
             // const value_type   val_at_ik = conj(a.Value(i));

             row_view_type &c = C.RowView(row_at_i);

             ordinal_type idx = 0;
             for (ordinal_type j=0;j<nnz_b && (idx > -2);++j) {
                const ordinal_type col_at_j  = b.Col(j);
                const value_type   val_at_kj = b.Value(j);

                idx = c.Index(col_at_j, idx);
                if (idx >= 0)
                  c.Value(idx) += alpha*val_at_ik*val_at_kj;
                }
          });
#else
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, nnz_a * nnz_b ),
          [&](const ordinal_type ii) {
             const ordinal_type i = ii / nnz_a ;
             const ordinal_type j = ii % nnz_a ;

             row_view_type &c = C.RowView( a.Col(i) );

             // Binary search for c's index of b.Col(j)
             const ordinal_type idx = c.Index( b.Col(j) );

             if (idx >= 0) {
               // const value_type   val_at_ik = conj(a.Value(i));
               c.Value(idx) += alpha * a.Value(i) * b.Value(j);
             }
          });
#endif

        member.team_barrier();
      }
    }

    return 0;
  }

}

#endif
