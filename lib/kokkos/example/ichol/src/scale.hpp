#pragma once
#ifndef __SCALE_HPP__
#define __SCALE_HPP__

/// \file scale.hpp
/// \brief Scaling sparse matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  using namespace std;

  template<typename T> struct ScaleTraits {
    typedef T scale_type;
    // assume built-in types have appropriate type conversion
    static constexpr T one = 1 ;
    static constexpr T zero = 0 ;
  };


  template<typename ScalarType,
           typename CrsExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  scaleCrsMatrix(const typename CrsExecViewType::policy_type::member_type &member,
                 const ScalarType alpha,
                 typename CrsExecViewType::matrix_type &A) {
    typedef typename CrsExecViewType::ordinal_type  ordinal_type;
    typedef typename CrsExecViewType::value_type    value_type;
    typedef typename CrsExecViewType::row_view_type row_view_type;

    if (alpha == ScaleTraits<value_type>::one) {
      // do nothing
    } else {
      const ordinal_type mA = A.NumRows();
      if (mA > 0) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, mA),
                             [&](const ordinal_type i) {
                               row_view_type &row = A.RowView(i);
                               for (ordinal_type j=0;j<row.NumNonZeros();++j)
                                 row.Value(j) *= alpha;
                             });
        member.team_barrier();
      }
    }

    return 0;
  }

  template<typename ScalarType,
           typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  int
  scaleDenseMatrix(const typename DenseExecViewType::policy_type::member_type &member,
                   const ScalarType alpha,
                   DenseExecViewType &A) {
    typedef typename DenseExecViewType::ordinal_type  ordinal_type;
    typedef typename DenseExecViewType::value_type    value_type;

    if (alpha == ScaleTraits<value_type>::one) {
      // do nothing
    } else {
      if (A.BaseObject().ColStride() > A.BaseObject().RowStride()) {
        const ordinal_type nA = A.NumCols();
        if (nA > 0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, nA),
                               [&](const ordinal_type j) {
                                 for (ordinal_type i=0;i<A.NumRows();++i)
                                   A.Value(i, j) *= alpha;
                               });
          member.team_barrier();
        }
      } else {
        const ordinal_type mA = A.NumRows();
        if (mA > 0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, mA),
                               [&](const ordinal_type i) {
                                 for (ordinal_type j=0;j<A.NumCols();++j)
                                   A.Value(i, j) *= alpha;
                               });
          member.team_barrier();
        }
      }
    }

    return 0;
  }

}

#endif

