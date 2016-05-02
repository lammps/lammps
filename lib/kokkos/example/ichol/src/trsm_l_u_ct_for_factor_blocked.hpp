#pragma once
#ifndef __TRSM_L_U_CT_FOR_FACTOR_BLOCKED_HPP__
#define __TRSM_L_U_CT_FOR_FACTOR_BLOCKED_HPP__

/// \file trsm_l_u_ct_for_factor_blocked.hpp
/// \brief Sparse triangular solve on given sparse patterns and multiple rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///

namespace Tacho {

  using namespace std;

  // Trsm used in the factorization phase: data parallel on b1t
  // ==========================================================
  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForFactorBlocked,Variant::One>
  ::invoke(typename CrsExecViewTypeA::policy_type &policy,
           const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           typename CrsExecViewTypeA::matrix_type &A,
           typename CrsExecViewTypeB::matrix_type &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;


if ( false && member.team_rank() == 0 ) {
 printf("Trsm [%d +%d)x[%d +%d)\n"
       , B.OffsetRows()
       , B.NumRows()
       , B.OffsetCols()
       , B.NumCols()
       );
}

    // scale the matrix B with alpha
    scaleCrsMatrix<ScalarType,CrsExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        // const value_type cdiag = std::conj(a.Value(0)); // for complex<T>
        const value_type cdiag = a.Value(0);

        // invert
        row_view_type &b1 = B.RowView(k);
        const ordinal_type nnz_b1 = b1.NumNonZeros();

        if (diagA != Diag::Unit && nnz_b1 > 0) {
          // b1t = b1t / conj(diag)
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nnz_b1),
                               [&](const ordinal_type j) {
                                 b1.Value(j) /= cdiag;
                               });
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nnz_b1),
                               [&](const ordinal_type j) {
                                 // grab b1
                                 const ordinal_type col_at_j = b1.Col(j);
                                 const value_type   val_at_j = b1.Value(j);

                                 for (ordinal_type i=1;i<nnz_a;++i) {
                                   // grab a12t
                                   const ordinal_type row_at_i = a.Col(i);
                                   // const value_type   val_at_i = conj(a.Value(i));
                                   const value_type   val_at_i = a.Value(i);

                                   // grab b2
                                   row_view_type &b2 = B.RowView(row_at_i);

                                   // check and update
                                   ordinal_type idx = 0;
                                   idx = b2.Index(col_at_j, idx);
                                   if (idx >= 0)
                                     b2.Value(idx) -= val_at_i*val_at_j;
                                 }
                               });
        }
        member.team_barrier();
      }
    }

    return 0;
  }

  // Trsm used in the factorization phase: data parallel on a1t
  // ==========================================================
  template<>
  template<typename ScalarType,
           typename CrsExecViewTypeA,
           typename CrsExecViewTypeB>
  KOKKOS_INLINE_FUNCTION
  int
  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,
       AlgoTrsm::ForFactorBlocked,Variant::Two>
  ::invoke(typename CrsExecViewTypeA::policy_type &policy,
           const typename CrsExecViewTypeA::policy_type::member_type &member,
           const int diagA,
           const ScalarType alpha,
           typename CrsExecViewTypeA::matrix_type &A,
           typename CrsExecViewTypeB::matrix_type &B) {
    typedef typename CrsExecViewTypeA::ordinal_type      ordinal_type;
    typedef typename CrsExecViewTypeA::value_type        value_type;
    typedef typename CrsExecViewTypeA::row_view_type     row_view_type;

    // scale the matrix B with alpha
    scaleCrsMatrix<ScalarType,CrsExecViewTypeB>(member, alpha, B);

    // Solve a system: AX = B -> B := inv(A) B
    const ordinal_type mA = A.NumRows();
    const ordinal_type nB = B.NumCols();

    if (nB > 0) {
      for (ordinal_type k=0;k<mA;++k) {
        row_view_type &a = A.RowView(k);
        // const value_type cdiag = conj(a.Value(0));
        const value_type cdiag = a.Value(0);

        // invert
        row_view_type &b1 = B.RowView(k);
        const ordinal_type nnz_b1 = b1.NumNonZeros();

        if (diagA != Diag::Unit && nnz_b1 > 0) {
          // b1t = b1t / conj(diag)
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nnz_b1),
                               [&](const ordinal_type j) {
                                 b1.Value(j) /= cdiag;
                               });
          member.team_barrier();
        }

        // update
        const ordinal_type nnz_a = a.NumNonZeros();
        if (nnz_a > 0) {
          // B2 = B2 - trans(conj(a12t)) b1t
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 1, nnz_a),
                               [&](const ordinal_type i) {
                                 // grab a12t
                                 const ordinal_type row_at_i = a.Col(i);
                                 // const value_type   val_at_i = conj(a.Value(i));
                                 const value_type   val_at_i = a.Value(i);

                                 // grab b2
                                 row_view_type &b2 = B.RowView(row_at_i);

                                 ordinal_type idx = 0;
                                 for (ordinal_type j=0;j<nnz_b1 && (idx > -2);++j) {
                                   // grab b1
                                   const ordinal_type col_at_j = b1.Col(j);
                                   const value_type   val_at_j = b1.Value(j);

                                   // check and update
                                   idx = b2.Index(col_at_j, idx);
                                   if (idx >= 0)
                                     b2.Value(idx) -= val_at_i*val_at_j;
                                 }
                               });
          member.team_barrier();
        }
      }
    }

    return 0;
  }

}

#endif
