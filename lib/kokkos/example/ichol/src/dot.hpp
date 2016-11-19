#pragma once
#ifndef __DOT_HPP__
#define __DOT_HPP__

/// \file dot.hpp
/// \brief Sparse dot product.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

/// dot_type result = x^H y

namespace Tacho { 

  using namespace std;

  template<typename T> struct DotTraits {
    typedef T dot_type;

    static KOKKOS_FORCEINLINE_FUNCTION 
    dot_type 
    // dot(const T &x, const T &y) { return conj<T>(x)*y; }
    dot(const T &x, const T &y) { return x*y; }
  }; 

  template<typename CrsRowViewType>
  KOKKOS_INLINE_FUNCTION 
  typename CrsRowViewType::value_type
  dot(const CrsRowViewType x, const CrsRowViewType y) {
    typedef typename CrsRowViewType::ordinal_type ordinal_type;
    typedef typename CrsRowViewType::value_type   value_type;

    typedef DotTraits<value_type> dot_traits;

    value_type r_val(0);

    const ordinal_type nnz_x = x.NumNonZeros();
    const ordinal_type nnz_y = y.NumNonZeros();

    for (ordinal_type jx=0, jy=0;jx<nnz_x && jy<nnz_y;) {
      const ordinal_type diff = x.Col(jx) - y.Col(jy);
      const ordinal_type sign = (0 < diff) - (diff < 0);
      switch (sign) {
      case  0:
        r_val += dot_traits::dot(x.Value(jx++), y.Value(jy++));
        break;
      case -1: ++jx; break;
      case  1: ++jy; break;
      }
    }
    
    return r_val;
  }

  template<typename CrsRowViewType>
  KOKKOS_INLINE_FUNCTION 
  typename CrsRowViewType::value_type
  dot(const CrsRowViewType x) {
    typedef typename CrsRowViewType::ordinal_type ordinal_type;
    typedef typename CrsRowViewType::value_type   value_type;

    typedef DotTraits<value_type> dot_traits;

    value_type r_val(0);

    const ordinal_type nnz = x.NumNonZeros();

    for (ordinal_type j=0;j<nnz;++j) 
      r_val += dot_traits::dot(x.Value(j), x.Value(j));
    
    return r_val;
  }

}

#endif
