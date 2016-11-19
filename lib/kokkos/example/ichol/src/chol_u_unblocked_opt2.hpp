#pragma once
#ifndef __CHOL_U_UNBLOCKED_OPT2_HPP__
#define __CHOL_U_UNBLOCKED_OPT2_HPP__

/// \file chol_u_unblocked_opt2.hpp
/// \brief Unblocked incomplete Chloesky factorization; version for data parallel sharing L1 cache.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "partition.hpp"

namespace Tacho {

  using namespace std;

  template<>
  template<typename CrsExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  Chol<Uplo::Upper,AlgoChol::UnblockedOpt,Variant::Two>
  ::invoke(typename CrsExecViewType::policy_type &policy,
           const typename CrsExecViewType::policy_type::member_type &member,
           typename CrsExecViewType::matrix_type &A) {

    typedef typename CrsExecViewType::value_type        value_type;
    typedef typename CrsExecViewType::ordinal_type      ordinal_type;
    typedef typename CrsExecViewType::row_view_type     row_view_type;

if ( false && member.team_rank() == 0 ) {
 printf("Chol [%d +%d)x[%d +%d) begin\n"
       , A.OffsetRows()
       , A.NumRows()
       , A.OffsetCols()
       , A.NumCols()
       );
}

    // row_view_type r1t, r2t;

    for (ordinal_type k=0;k<A.NumRows();++k) {
      //r1t.setView(A, k);
      row_view_type &r1t = A.RowView(k);

      // extract diagonal from alpha11
      value_type &alpha = r1t.Value(0);

      if (member.team_rank() == 0) {
        // if encounter null diag or wrong index, return -(row + 1)
        if (abs(alpha) == 0.0 || r1t.Col(0) != k)
          return -(k + 1);

        // error handling should be more carefully designed

        // sqrt on diag
        // alpha = sqrt(real(alpha));
        alpha = sqrt(alpha);
      }
      member.team_barrier();


if ( false && member.team_rank() == 0 ) {
 printf("Chol [%d +%d)x[%d +%d) local row %d\n"
       , A.OffsetRows()
       , A.NumRows()
       , A.OffsetCols()
       , A.NumCols()
       , int(k)
       );
}


      const ordinal_type nnz_r1t = r1t.NumNonZeros();

      if (nnz_r1t) {
        // inverse scale
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 1, nnz_r1t),
                             [&](const ordinal_type j) {
                               r1t.Value(j) /= alpha;
                             });

        member.team_barrier();


if ( false && member.team_rank() == 0 ) {
 printf("Chol [%d +%d)x[%d +%d) local row %d nnz_r1t\n"
       , A.OffsetRows()
       , A.NumRows()
       , A.OffsetCols()
       , A.NumCols()
       , int(k)
       );
}

        // hermitian rank update
        for (ordinal_type i=1;i<nnz_r1t;++i) {
          const ordinal_type row_at_i = r1t.Col(i);
          // const value_type   val_at_i = conj(r1t.Value(i));
          const value_type   val_at_i = r1t.Value(i);

          //r2t.setView(A, row_at_i);
          row_view_type &r2t = A.RowView(row_at_i);

          ordinal_type member_idx = 0 ;

          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, i, nnz_r1t),
                               [&](const ordinal_type j) {
                                 if (member_idx > -2) {
                                   const ordinal_type col_at_j = r1t.Col(j);
                                   member_idx = r2t.Index(col_at_j, member_idx);
                                   if (member_idx >= 0) {
                                     const value_type   val_at_j = r1t.Value(j);
                                     r2t.Value(member_idx) -= val_at_i*val_at_j;
                                   }
                                 }
                               });
        }
      }


if ( false ) {
member.team_barrier();
if ( member.team_rank() == 0 ) {
 printf("Chol [%d +%d)x[%d +%d) local row %d end\n"
       , A.OffsetRows()
       , A.NumRows()
       , A.OffsetCols()
       , A.NumCols()
       , int(k)
       );
}
}

    }


if ( false ) {
member.team_barrier();
if ( member.team_rank() == 0 ) {
 printf("Chol [%d +%d)x[%d +%d) end\n"
       , A.OffsetRows()
       , A.NumRows()
       , A.OffsetCols()
       , A.NumCols()
       );
}
}


    return 0;
  }

}

#endif
