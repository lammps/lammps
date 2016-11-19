#pragma once
#ifndef __CRS_MATRIX_BASE_HPP__
#define __CRS_MATRIX_BASE_HPP__

/// \file crs_matrix_base.hpp
/// \brief CRS matrix base object interfaces to user provided input matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"
#include "coo.hpp"

namespace Tacho { 

  using namespace std;

  template< typename , typename > class TaskView ;

  template < typename CrsMatrixType >
  struct GetCrsMatrixRowViewType {
    typedef int type ;
  };


  template < typename CrsMatrixViewType , typename TaskFactoryType >
  struct GetCrsMatrixRowViewType
    < TaskView<CrsMatrixViewType,TaskFactoryType> >
  {
    typedef typename CrsMatrixViewType::row_view_type type ;
  };

  /// \class CrsMatrixBase
  /// \breif CRS matrix base object using Kokkos view and subview
  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType = OrdinalType,
           typename SpaceType = void,
           typename MemoryTraits = void>
  class CrsMatrixBase {
  public:
    typedef ValueType    value_type;
    typedef OrdinalType  ordinal_type;
    typedef SpaceType    space_type;
    typedef SizeType     size_type;
    typedef MemoryTraits memory_traits;

    // 1D view, layout does not matter; no template parameters for that
    typedef Kokkos::View<size_type*,   space_type,memory_traits> size_type_array;
    typedef Kokkos::View<ordinal_type*,space_type,memory_traits> ordinal_type_array;
    typedef Kokkos::View<value_type*,  space_type,memory_traits> value_type_array;

    typedef typename size_type_array::value_type*    size_type_array_ptr;
    typedef typename ordinal_type_array::value_type* ordinal_type_array_ptr;
    typedef typename value_type_array::value_type*   value_type_array_ptr;

    // range type
    template<typename T> using range_type = pair<T,T>;

    // external interface
    typedef Coo<CrsMatrixBase> ijv_type;
    
    friend class CrsMatrixHelper;

  private:

    ordinal_type       _m;       //!< # of rows
    ordinal_type       _n;       //!< # of cols
    size_type          _nnz;     //!< # of nonzeros
    size_type_array    _ap;      //!< pointers to column index and values
    ordinal_type_array _aj;      //!< column index compressed format
    value_type_array   _ax;      //!< values

  public:

    typedef typename GetCrsMatrixRowViewType< ValueType >::type row_view_type ;
    typedef Kokkos::View<row_view_type*,space_type> row_view_type_array;

    row_view_type_array _all_row_views ;

  protected:

    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n,
                              const size_type nnz) {
      _m = m;
      _n = n;
      _nnz = nnz;

      if (static_cast<ordinal_type>(_ap.dimension_0()) < m+1)
        _ap = size_type_array("CrsMatrixBase::RowPtrArray", m+1);
      
      if (static_cast<size_type>(_aj.dimension_0()) < nnz)
        _aj = ordinal_type_array("CrsMatrixBase::ColsArray", nnz);

      if (static_cast<size_type>(_ax.dimension_0()) < nnz)
        _ax = value_type_array("CrsMatrixBase::ValuesArray", nnz);
    }

    // Copy sparse matrix structure from coordinate format in 'mm'
    // to CRS format in Views _ap, _aj, a_x.
    void ijv2crs(const vector<ijv_type> &mm) {

      ordinal_type ii = 0;
      size_type jj = 0;
      
      ijv_type prev = mm[0];
      _ap[ii++] = 0;
      _aj[jj] = prev.Col();
      _ax[jj] = prev.Val();
      ++jj;
      
      for (typename vector<ijv_type>::const_iterator it=(mm.begin()+1);it<mm.end();++it) {
        ijv_type aij = (*it);
        
        // row index
        if (aij.Row() != prev.Row()) {
          _ap[ii++] = jj; 
        }
        
        if (aij == prev) {
          --jj;
          _aj[jj]  = aij.Col();
          _ax[jj] += aij.Val();
        } else {
          _aj[jj] = aij.Col();
          _ax[jj] = aij.Val();
        }
        ++jj;
        
        prev = aij;
      }
      
      // add the last index to terminate the storage
      _ap[ii++] = jj;
      _nnz = jj;
    }
    
  public:

    KOKKOS_INLINE_FUNCTION
    void setNumNonZeros() { 
      if (_m) 
        _nnz = _ap[_m];
    }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumRows() const { return _m; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }

    KOKKOS_INLINE_FUNCTION
    size_type NumNonZeros() const { return _nnz; }

    KOKKOS_INLINE_FUNCTION
    size_type_array_ptr RowPtr() const { return &_ap[0]; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type_array_ptr ColPtr() const { return &_aj[0]; }

    KOKKOS_INLINE_FUNCTION
    value_type_array_ptr ValuePtr() const { return &_ax[0];}

    KOKKOS_INLINE_FUNCTION
    size_type RowPtr(const ordinal_type i) const { return _ap[i]; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type_array_ptr ColsInRow(const ordinal_type i) const { return _aj.data() + _ap[i] ; }
    
    KOKKOS_INLINE_FUNCTION
    value_type_array_ptr ValuesInRow(const ordinal_type i) const { return _ax.data() + _ap[i] ; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type NumNonZerosInRow(const ordinal_type i) const { return (_ap[i+1] - _ap[i]); } 

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type k) { return _ax[k]; }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type k) const { return _ax[k]; }

    /// \brief Default constructor.
    KOKKOS_INLINE_FUNCTION
    CrsMatrixBase() 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap(),
        _aj(),
        _ax()
    { }

    /// \brief Constructor with label
    CrsMatrixBase(const string & ) 
      : _m(0),
        _n(0),
        _nnz(0),
        _ap(),
        _aj(),
        _ax()
    { }

    /// \brief Copy constructor (shallow copy), for deep-copy use a method copy
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    CrsMatrixBase(const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) 
      : _m(b._m),
        _n(b._n),
        _nnz(b._nnz),
        _ap(b._ap), 
        _aj(b._aj),
        _ax(b._ax) 
    { }

    /// \brief Constructor to allocate internal data structures.
    CrsMatrixBase(const string & ,
                  const ordinal_type m, 
                  const ordinal_type n, 
                  const ordinal_type nnz) 
      : _m(m),
        _n(n),
        _nnz(nnz),
        _ap("CrsMatrixBase::RowPtrArray", m+1),
        _aj("CrsMatrixBase::ColsArray", nnz),
        _ax("CrsMatrixBase::ValuesArray", nnz)
    { }

    /// \brief Constructor to attach external arrays to the matrix.
    CrsMatrixBase(const string &,
                  const ordinal_type m, 
                  const ordinal_type n, 
                  const ordinal_type nnz,
                  const size_type_array &ap,
                  const ordinal_type_array &aj,
                  const value_type_array &ax) 
      : _m(m),
        _n(n),
        _nnz(nnz),
        _ap(ap), 
        _aj(aj),
        _ax(ax) 
    { }
    
  // Allow the copy function access to the input CrsMatrixBase
  // private data.
  template<typename, typename, typename, typename, typename>
  friend class CrsMatrixBase ;

  public:
    /// \brief deep copy of matrix b, potentially different spaces
    template< typename SpT >
    int 
    copy(const CrsMatrixBase<ValueType,OrdinalType,SizeType,SpT,MemoryTraits> &b) {

      space_type::execution_space::fence();

      createInternalArrays(b._m, b._n, b._nnz);

      space_type::execution_space::fence();

      const auto ap_range = range_type<ordinal_type>(0, min(_ap.dimension_0(), b._ap.dimension_0()));
      const auto aj_range = range_type<size_type>   (0, min(_aj.dimension_0(), b._aj.dimension_0()));
      const auto ax_range = range_type<size_type>   (0, min(_ax.dimension_0(), b._ax.dimension_0()));

      Kokkos::deep_copy(Kokkos::subview(  _ap, ap_range), 
                        Kokkos::subview(b._ap, ap_range));
      Kokkos::deep_copy(Kokkos::subview(  _aj, aj_range),
                        Kokkos::subview(b._aj, aj_range));

      Kokkos::deep_copy(Kokkos::subview(  _ax, ax_range),
                        Kokkos::subview(b._ax, ax_range));

      space_type::execution_space::fence();

      return 0;
    }

    /// \brief deep copy of lower/upper triangular of matrix b
    int 
    copy(const int uplo, 
         const CrsMatrixBase &b) { 

      createInternalArrays(b._m, b._n, b._nnz);

      // assume that matrix b is sorted.
      switch (uplo) {
      case Uplo::Lower: {
        _nnz = 0;
        for (ordinal_type i=0;i<_m;++i) {
          size_type jbegin = b._ap[i];
          size_type jend   = b._ap[i+1];
          _ap[i] = _nnz;
          for (size_type j=jbegin;j<jend && (i >= b._aj[j]);++j,++_nnz) {
            _aj[_nnz] = b._aj[j];
            _ax[_nnz] = b._ax[j]; 
          }
        }
        _ap[_m] = _nnz;
        break;
      }
      case Uplo::Upper: {
        _nnz = 0;
        for (ordinal_type i=0;i<_m;++i) {
          size_type j = b._ap[i];
          size_type jend = b._ap[i+1];
          _ap[i] = _nnz;
          for ( ;j<jend && (i > b._aj[j]);++j) ;
          for ( ;j<jend;++j,++_nnz) {
            _aj[_nnz] = b._aj[j];
            _ax[_nnz] = b._ax[j]; 
          }
        }
        _ap[_m] = _nnz;
        break;
      }
      }

      return 0;
    }

    /// \brief deep copy of matrix b with given permutation vectors
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    int
    copy(const typename CrsMatrixBase<VT,OT,ST,SpT,MT>::ordinal_type_array &p,
         const typename CrsMatrixBase<VT,OT,ST,SpT,MT>::ordinal_type_array &ip,
         const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) {

      createInternalArrays(b._m, b._n, b._nnz);

      // Question:: do I need to use Kokkos::vector ? 
      //            in other words, where do we permute matrix in factoriztion ?
      //            permuting a matrix is a kernel ? 
      vector<ijv_type> tmp;

      // any chance to use parallel_for ?
      _nnz = 0;
      for (ordinal_type i=0;i<_m;++i) {
        ordinal_type ii = ip[i];

        size_type jbegin = b._ap[ii];
        size_type jend   = b._ap[ii+1];

        _ap[i] = _nnz;
        for (size_type j=jbegin;j<jend;++j) {
          ordinal_type jj = p[b._aj[j]];
          ijv_type aij(i, jj, b._ax[j]);
          tmp.push_back(aij);
        }

        sort(tmp.begin(), tmp.end(), less<ijv_type>());
        for (auto it=tmp.begin();it<tmp.end();++it) {
          ijv_type aij = (*it);

          _aj[_nnz] = aij.Col();
          _ax[_nnz] = aij.Val();
          ++_nnz;
        }
        tmp.clear();
      }
      _ap[_m] = _nnz;

      return 0;
    }

    /// \brief add the matrix b into this non-zero entires
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT,
             typename MT>
    int 
    add(const CrsMatrixBase<VT,OT,ST,SpT,MT> &b) { 

      const ordinal_type m = min(b._m, _m);
      for (ordinal_type i=0;i<m;++i) {
        const size_type jaend = _ap[i+1];
        const size_type jbend = b._ap[i+1];

        size_type ja = _ap[i];
        size_type jb = b._ap[i];
        
        for ( ;jb<jbend;++jb) {
          for ( ;(_aj[ja]<b._aj[jb] && ja<jaend);++ja);
          _ax[ja] += (_aj[ja] == b._aj[jb])*b._ax[jb];
        }
      }

      return 0;
    }

    int symmetrize(const int uplo, 
                   const bool conjugate = false) {
      vector<ijv_type> mm;
      mm.reserve(_nnz*2);

      for (ordinal_type i=0;i<_m;++i) {
        const size_type jbegin = _ap[i];
        const size_type jend   = _ap[i+1];
        for (size_type jj=jbegin;jj<jend;++jj) {
          const ordinal_type j = _aj[jj];
          const value_type val = (conjugate ? conj(_ax[j]) : _ax[j]);
          if        (uplo == Uplo::Lower && i > j) {
            mm.push_back(ijv_type(i, j, val));
            mm.push_back(ijv_type(j, i, val));
          } else if (uplo == Uplo::Upper && i < j) {
            mm.push_back(ijv_type(i, j, val));
            mm.push_back(ijv_type(j, i, val));
          } else if (i == j) {
            mm.push_back(ijv_type(i, i, val));
          }
        }
      }
      sort(mm.begin(), mm.end(), less<ijv_type>());

      createInternalArrays(_m, _n, mm.size());
      
      ijv2crs(mm);
      
      return 0;
    }

    int hermitianize(int uplo) {
      return symmetrize(uplo, true);
    }

    ostream& showMe(ostream &os) const {
      streamsize prec = os.precision();
      os.precision(8);
      os << scientific;

      os << " -- CrsMatrixBase -- " << endl
         << "    # of Rows          = " << _m << endl
         << "    # of Cols          = " << _n << endl
         << "    # of NonZeros      = " << _nnz << endl
         << endl
         << "    RowPtrArray length = " << _ap.dimension_0() << endl
         << "    ColArray    length = " << _aj.dimension_0() << endl 
         << "    ValueArray  length = " << _ax.dimension_0() << endl
         << endl;
      
      const int w = 10;
      if (_ap.size() && _aj.size() && _ax.size()) {
        os << setw(w) <<  "Row" << "  " 
           << setw(w) <<  "Col" << "  " 
           << setw(w) <<  "Val" << endl;
        for (ordinal_type i=0;i<_m;++i) {
          size_type jbegin = _ap[i], jend = _ap[i+1];
          for (size_type j=jbegin;j<jend;++j) {
            value_type val = _ax[j];
            os << setw(w) <<      i << "  " 
               << setw(w) << _aj[j] << "  " 
               << setw(w) <<    val << endl;
          }
        }
      }

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

    int importMatrixMarket(ifstream &file) {

      vector<ijv_type> mm; 
      const ordinal_type mm_base = 1; 

      {
        string header;
        if (file.is_open()) {
          getline(file, header);
          while (file.good()) {
            char c = file.peek();
            if (c == '%' || c == '\n') {
              file.ignore(256, '\n');
            continue;
            }
            break;
          }
        } else {
          ERROR(MSG_INVALID_INPUT(file));
        }

        // check the header
        bool symmetry = (header.find("symmetric") != string::npos);

        // read matrix specification
        ordinal_type m, n;
        size_type nnz;
        
        file >> m >> n >> nnz;
        
        mm.reserve(nnz*(symmetry ? 2 : 1));
        for (size_type i=0;i<nnz;++i) {
          ordinal_type row, col;
          value_type val;
          file >> row >> col >> val;
          
          row -= mm_base;
          col -= mm_base;
          
          mm.push_back(ijv_type(row, col, val));
          if (symmetry && row != col)
            mm.push_back(ijv_type(col, row, val));
        }
        sort(mm.begin(), mm.end(), less<ijv_type>());
      
        // construct workspace and set variables
        createInternalArrays(m, n, mm.size());
      }
      
      // change mm to crs
      ijv2crs(mm);
      
      return 0;
    }
    
    int exportMatrixMarket(ofstream &file,
                           const string comment,
                           const int uplo = 0) {
      streamsize prec = file.precision();
      file.precision(8);
      file << scientific;

      file << "%%MatrixMarket matrix coordinate "
           << (is_fundamental<value_type>::value ? "real " : "complex ")
           << ((uplo == Uplo::Upper || uplo == Uplo::Lower) ? "symmetric " : "general ")
           << endl;

      file << comment << endl;
      
      // cnt nnz
      size_type nnz = 0;
      for (ordinal_type i=0;i<_m;++i) {
        const size_type jbegin = _ap[i], jend = _ap[i+1];
        for (size_type j=jbegin;j<jend;++j) {
          if (uplo == Uplo::Upper && i <= _aj[j]) ++nnz;
          if (uplo == Uplo::Lower && i >= _aj[j]) ++nnz;
          if (!uplo) ++nnz;
        }
      }
      file << _m << " " << _n << " " << nnz << endl;

      const int w = 10;
      for (ordinal_type i=0;i<_m;++i) {
        const size_type jbegin = _ap[i], jend = _ap[i+1];
        for (size_type j=jbegin;j<jend;++j) {
          bool flag = false;
          if (uplo == Uplo::Upper && i <= _aj[j]) flag = true;
          if (uplo == Uplo::Lower && i >= _aj[j]) flag = true;
          if (!uplo) flag = true;
          if (flag) {
            value_type val = _ax[j];
            file << setw(w) << (     i+1) << "  " 
                 << setw(w) << (_aj[j]+1) << "  " 
                 << setw(w) <<    val << endl;
          }
        }
      }

      file.unsetf(ios::scientific);
      file.precision(prec);

      return 0;
    }

    //----------------------------------------------------------------------

    int convertGraph(size_type_array rptr,
                     ordinal_type_array cidx) const {
      ordinal_type ii = 0;
      size_type jj = 0;

      for (ordinal_type i=0;i<_m;++i) {
        size_type jbegin = _ap[i], jend = _ap[i+1];
        rptr[ii++] = jj;
        for (size_type j=jbegin;j<jend;++j)
          if (i != _aj[j])
            cidx[jj++] = _aj[j];
      }
      rptr[ii] = jj;

      return 0;
    }

    //----------------------------------------------------------------------

  };

}

#endif
