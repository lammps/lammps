#pragma once
#ifndef __NORM_HPP__
#define __NORM_HPP__

/// \file norm.hpp
/// \brief Compute norm of sparse or dense matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  using namespace std;

  template<typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  auto
  normOneDenseMatrix(DenseExecViewType &A) -> decltype(real(typename DenseExecViewType::value_type())) {
    typedef typename DenseExecViewType::ordinal_type  ordinal_type;
    typedef typename DenseExecViewType::value_type    value_type;
    typedef decltype(real(value_type())) norm_type;

    const ordinal_type mA = A.NumRows();
    const ordinal_type nA = A.NumCols();

    norm_type r_val = 0.0;

    for (ordinal_type j=0;j<nA;++j) {
      norm_type col_sum_at_j = 0.0;
      for (ordinal_type i=0;i<mA;++i)
        col_sum_at_j += abs(A.Value(i,j));
      r_val = max(r_val, col_sum_at_j);
    }
    return r_val;
  }

  template<typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  auto
  normInfDenseMatrix(DenseExecViewType &A) -> decltype(real(typename DenseExecViewType::value_type())) {
    typedef typename DenseExecViewType::ordinal_type  ordinal_type;
    typedef typename DenseExecViewType::value_type    value_type;
    typedef decltype(real(value_type())) norm_type;

    const ordinal_type mA = A.NumRows();
    const ordinal_type nA = A.NumCols();

    norm_type r_val = 0.0;

    for (ordinal_type i=0;i<mA;++i) {
      norm_type row_sum_at_i = 0.0;
      for (ordinal_type j=0;j<nA;++j) 
        row_sum_at_i += abs(A.Value(i,j));
      r_val = max(r_val, row_sum_at_i);
    }
    return r_val;
  }
  
  template<typename DenseExecViewType>
  KOKKOS_INLINE_FUNCTION
  auto
  normFrobeniusDenseMatrix(DenseExecViewType &A) -> decltype(real(typename DenseExecViewType::value_type())) {
    typedef typename DenseExecViewType::ordinal_type  ordinal_type;
    typedef typename DenseExecViewType::value_type    value_type;
    typedef decltype(real(value_type())) norm_type;

    const ordinal_type mA = A.NumRows();
    const ordinal_type nA = A.NumCols();

    norm_type r_val = 0.0;

    for (ordinal_type i=0;i<mA;++i) 
      for (ordinal_type j=0;j<nA;++j) {
        value_type val = A.Value(i,j);
        // r_val += conj(val)*val;
        r_val += val*val;
      }
    return sqrt(r_val);
  }

}

#endif

