#pragma once
#ifndef __CRS_ROW_VIEW_HPP__
#define __CRS_ROW_VIEW_HPP__

/// \file crs_row_view.hpp
/// \brief A view to a row extracted from CrsMatrixView.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  using namespace std;

  /// \class CrsRowView
  template<typename CrsMatBaseType>
  class CrsRowView {
  public:
    typedef typename CrsMatBaseType::ordinal_type           ordinal_type;
    typedef typename CrsMatBaseType::value_type             value_type;
    typedef typename CrsMatBaseType::ordinal_type_array_ptr ordinal_type_array_ptr;
    typedef typename CrsMatBaseType::value_type_array_ptr   value_type_array_ptr;
    
  private:
    // row info
    ordinal_type _offn, _n;    

    // this assumes a contiguous memory buffer
    ordinal_type_array_ptr _aj, _ajn; // column index compressed format in row
    value_type_array_ptr   _ax;                // values 

    static KOKKOS_INLINE_FUNCTION
    typename CrsMatBaseType::ordinal_type_array_ptr
    lower_bound( typename CrsMatBaseType::ordinal_type_array_ptr begin ,
                 typename CrsMatBaseType::ordinal_type_array_ptr const end ,
                 typename CrsMatBaseType::ordinal_type           const val )
      {
         typename CrsMatBaseType::ordinal_type_array_ptr it = begin ;
         int count = end - begin ;
         int step = 0 ;
         while (count>0) {
           it = begin ;
           it += ( step = (count >> 1) );
           if (*it<val) {
             begin=++it;
             count-=step+1;
           }
           else { count=step; }
         }
         return begin;
      }

  public:
    KOKKOS_INLINE_FUNCTION
    ordinal_type OffsetCols() const { return _offn; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZeros() const { return _ajn - _aj; } 

    KOKKOS_INLINE_FUNCTION
    ordinal_type Col(const ordinal_type j) const { return _aj[j] - _offn; }

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type j) { return _ax[j]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type j) const { return _ax[j]; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type Index(const ordinal_type col ) const {
      const ordinal_type loc = _offn + col ;
      // binary search
      ordinal_type_array_ptr aj = CrsRowView::lower_bound(_aj, _ajn, loc);

      // if found, return index for the location, 
      // otherwise return -1 (not found), -2 (end of array)
      return (aj < _ajn ? (*aj == loc ? aj - _aj : -1) : -2);
    }

    KOKKOS_INLINE_FUNCTION
    ordinal_type Index(const ordinal_type col,
                       const ordinal_type prev ) const {
      const ordinal_type loc = _offn + col;
      ordinal_type_array_ptr aj = _aj + prev;

      // binary search
      // aj = lower_bound(aj, _ajn, loc);

      // linear search from prev: this gains about 45 % faster
      for ( ;aj < _ajn && *aj<loc; ++aj); 

      // if found, return index for the location, 
      // otherwise return -1 (not found), -2 (end of array)
      return (aj < _ajn ? (*aj == loc ? aj - _aj : -1) : -2);
    }

    KOKKOS_INLINE_FUNCTION
    value_type ValueAtColumn(const ordinal_type col) const {
      const ordinal_type j = Index(col);
      return (j < 0 ? value_type(0) : _ax[j]);
    }

    KOKKOS_INLINE_FUNCTION
    CrsRowView()
      : _offn(0),
        _n(0),
        _aj(),
        _ajn(),
        _ax() 
    { }


    KOKKOS_INLINE_FUNCTION
    CrsRowView(const ordinal_type           offn,
               const ordinal_type           n,
               const ordinal_type_array_ptr aj,
               const ordinal_type_array_ptr ajn,
               const value_type_array_ptr   ax) 
      : _offn(offn),
        _n(n),
        _aj(aj),
        _ajn(ajn),
        _ax(ax) 
    { }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(const CrsMatrixView<CrsMatBaseType> &A, 
               const ordinal_type i) {
      this->setView(A, i);
    }

    KOKKOS_INLINE_FUNCTION
    CrsRowView(const CrsMatBaseType &A, 
               const ordinal_type i) {
      this->setView(A, i);
    }

    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatrixView<CrsMatBaseType> &A, 
                 const ordinal_type i) {
      _offn = A.OffsetCols();
      _n    = A.NumCols();

      const ordinal_type ii = A.OffsetRows() + i;

      const typename CrsMatBaseType::ordinal_type_array_ptr cols = A.BaseObject().ColsInRow(ii);
      const typename CrsMatBaseType::ordinal_type_array_ptr next = A.BaseObject().ColsInRow(ii+1);
      const typename CrsMatBaseType::value_type_array_ptr   vals = A.BaseObject().ValuesInRow(ii);

      // [cols..next) is sorted so a log(N) search could performed
      _aj  = CrsRowView::lower_bound(cols, next, _offn);
      _ajn = CrsRowView::lower_bound(_aj,  next, _offn+_n);

      _ax  = &vals[_aj - cols];
    }

    KOKKOS_INLINE_FUNCTION
    void setView(const CrsMatBaseType &A, 
                 const ordinal_type i) {
      _offn = 0;
      _n    = A.NumCols();
      _aj   = A.ColsInRow(i);
      _ajn  = A.ColsInRow(i+1);
      _ax   = A.ValuesInRow(i);
    }

    ostream& showMe(ostream &os) const {                                                
      const ordinal_type nnz = NumNonZeros();
      const ordinal_type offset = OffsetCols();
      os << "  offset = " << offset
         << ", nnz = " << nnz
         << endl; 
      for (ordinal_type j=0;j<nnz;++j) {
        const value_type val = _ax[j];
        os << "(" << Col(j) << ", "
           << val << ")"
           << endl;
      }
      return os;
    }
  };
}

#endif
