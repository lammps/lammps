
#ifndef __CRS_MATRIX_HELPER_IMPL_HPP__
#define __CRS_MATRIX_HELPER_IMPL_HPP__

/// \file crs_matrix_helper_impl.hpp
/// \brief This file includes utility functions to convert between flat and hierarchical matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"

namespace Tacho {

  using namespace std;

  template< typename CrsHierBase >
  struct FunctorFillRowViewArray {

    typedef typename CrsHierBase::ordinal_type         ordinal_type ;
    typedef typename CrsHierBase::row_view_type_array  row_view_type_array ;
    typedef typename CrsHierBase::value_type_array     ax_type ;

    typedef ordinal_type value_type ;

    row_view_type_array _all_row_views ;
    ax_type             _ax ;

    FunctorFillRowViewArray( const row_view_type_array & arg_all_row_views
                           , const ax_type             & arg_ax )
      : _all_row_views( arg_all_row_views )
      , _ax( arg_ax )
      {}

    KOKKOS_INLINE_FUNCTION
    void operator()( ordinal_type k , ordinal_type & value ) const
      { value += _ax(k).NumRows(); }

    KOKKOS_INLINE_FUNCTION
    void operator()( ordinal_type k , ordinal_type & value , bool final ) const
      {
        if ( final ) {
          const int begin = value ;
          const int end   = begin + _ax(k).NumRows();

          auto sub = Kokkos::subview( _all_row_views, Kokkos::pair<int,int>(begin,end) );

          _ax(k).setRowViewArray( sub );
        }

        value += _ax(k).NumRows();
      }
  };

  template< typename CrsHierBase >
  int CrsMatrixHelper::fillRowViewArray( CrsHierBase & device_HU )
  {
    typedef typename CrsHierBase::row_view_type_array row_view_type_array ;
    typedef typename CrsHierBase::space_type          space_type ;

    ordinal_type total_row_view_count = 0 ;

    Kokkos::RangePolicy< space_type >
      range_policy( 0 , device_HU.NumNonZeros() );

    space_type::fence();

    {
      FunctorFillRowViewArray< CrsHierBase >
         functor( row_view_type_array() , device_HU._ax );


      Kokkos::parallel_reduce( range_policy , functor , total_row_view_count );
    }

    device_HU._all_row_views =
      row_view_type_array("RowViews",total_row_view_count);

    space_type::fence();

    {
      FunctorFillRowViewArray< CrsHierBase >
         functor( device_HU._all_row_views , device_HU._ax );

      Kokkos::parallel_scan( range_policy , functor );
    }

    space_type::fence();

    return 0 ;
  }
  
  template<typename CrsFlatBase>
  int
  CrsMatrixHelper::filterZeros(CrsFlatBase &flat) {
    typedef typename CrsFlatBase::ordinal_type           ordinal_type;
    typedef typename CrsFlatBase::size_type              size_type;
    typedef typename CrsFlatBase::value_type             value_type;
    
    typedef typename CrsFlatBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    typedef typename CrsFlatBase::value_type_array_ptr   value_type_array_ptr;
    
    size_type nz = 0;
    const value_type zero(0);
    
    for (ordinal_type k=0;k<flat.NumNonZeros();++k) 
      nz += (flat.Value(k) == zero) ;
    
    if (nz) {
      CrsFlatBase resized(flat.Label() + "::ZeroFiltered", 
                          flat.NumRows(),
                          flat.NumCols(),
                          flat.NumNonZeros() - nz);
      
      ordinal_type_array_ptr rows = resized.RowPtr(); rows[0] = 0;
      ordinal_type_array_ptr cols = resized.ColPtr();
      value_type_array_ptr vals = resized.ValuePtr();    
      
      size_type nnz = 0;
      for (ordinal_type i=0;i<flat.NumRows();++i) {
        const ordinal_type nnz_in_row = flat.NumNonZerosInRow(i);
        const ordinal_type_array_ptr cols_in_row = flat.ColsInRow(i);
        const value_type_array_ptr vals_in_row = flat.ValuesInRow(i);
        
        for (ordinal_type j=0;j<nnz_in_row;++j) {
          if (vals_in_row[j] != zero) {
            cols[nnz] = cols_in_row[j];
            vals[nnz] = vals_in_row[j];
            ++nnz;
          }
        }
        rows[i+1] = nnz;
      }
      flat = resized;
      resized.setNumNonZeros();
    }

    return 0;
  }


  template<typename CrsFlatBase,
           typename CrsHierBase>
  int
  CrsMatrixHelper::flat2hier(CrsFlatBase &flat,
                             CrsHierBase &hier) {
    typedef typename CrsHierBase::ordinal_type           ordinal_type;
    typedef typename CrsHierBase::size_type              size_type;
    typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;

    size_type nnz = 0;

    hier.createInternalArrays(flat.NumRows(), flat.NumCols(), flat.NumNonZeros());

    for (ordinal_type i=0;i<flat.NumRows();++i) {
      ordinal_type jsize = flat.NumNonZerosInRow(i);

      hier._ap[i] = nnz;
      ordinal_type_array_ptr ci = flat.ColsInRow(i);
      for (ordinal_type j=0;j<jsize;++j,++nnz) {
        hier._aj[nnz] = ci[j];
        hier._ax[nnz].setView( flat,     i, 1,
                              /**/   ci[j], 1);
      }
    }

    hier._ap[flat.NumRows()] = nnz;
    hier._nnz = nnz;

    return 0;
  }

  template< typename CrsFlatBase ,
            typename CrsHierBase ,
            typename HostOrdinalTypeArray >
  int
  CrsMatrixHelper::flat2hier(int uplo,
                             CrsFlatBase &flat,
                             CrsHierBase &hier,
                             const typename CrsHierBase::ordinal_type       nblks,
                             const HostOrdinalTypeArray range ,
                             const HostOrdinalTypeArray tree) {
    switch(uplo) {
    case Uplo::Upper: return flat2hier_upper(flat, hier, nblks, range, tree);
    case Uplo::Lower: return flat2hier_lower(flat, hier, nblks, range, tree);
    }
    return -1;
  }

  template<typename CrsFlatBase,
           typename CrsHierBase,
           typename HostOrdinalTypeArray >
  int
  CrsMatrixHelper::flat2hier_upper(CrsFlatBase & device_flat, 
                                   CrsHierBase & device_hier,
                                   const typename CrsHierBase::ordinal_type       nblks,
                                   const HostOrdinalTypeArray range,
                                   const HostOrdinalTypeArray tree) {
    typedef typename CrsHierBase::ordinal_type            ordinal_type;
    typedef typename CrsHierBase::size_type               size_type;
    
    //typedef typename CrsHierBase::ordinal_type_array     ordinal_type_array;
    //typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    //typedef typename CrsHierBase::value_type_array_ptr   value_type_array_ptr;
    
    size_type nnz = 0;
    
    // count nnz and nnz in rows for the upper triangular hier matrix
    for (ordinal_type i=0;i<nblks;++i) 
      for (ordinal_type j=i;j != -1;++nnz,j=tree[j]) ;
    
    // create upper triangular block matrix
    device_hier.createInternalArrays(nblks, nblks, nnz);    

    typename CrsHierBase::size_type_array::HostMirror
      host_ap = Kokkos::create_mirror_view( device_hier._ap );

    typename CrsHierBase::ordinal_type_array::HostMirror
      host_aj = Kokkos::create_mirror_view( device_hier._aj );

    typename CrsHierBase::value_type_array::HostMirror
      host_ax = Kokkos::create_mirror_view( device_hier._ax );

    nnz = 0;
    for (ordinal_type i=0;i<nblks;++i) {
      host_ap[i] = nnz;
      for (ordinal_type j=i;j != -1;++nnz,j=tree[j]) {
        host_aj[nnz] = j;
        host_ax[nnz].setView( device_flat, range[i], (range[i+1] - range[i]),
                             /**/          range[j], (range[j+1] - range[j]));

        // this checking might more expensive 
        // and attempts to access device memory from the host
        // if (!host_ax[nnz].countNumNonZeros())
        //  --nnz;
      }
    }
    
    host_ap[nblks] = nnz;

    Kokkos::deep_copy( device_hier._ap , host_ap );
    Kokkos::deep_copy( device_hier._aj , host_aj );
    Kokkos::deep_copy( device_hier._ax , host_ax );

    device_hier._nnz = nnz;

    return 0;
  }

  // template<typename CrsFlatBase,
  //          typename CrsHierBase>
  // int
  // CrsMatrixHelper::flat2hier_upper(CrsFlatBase &flat,
  //                                  CrsHierBase &hier,
  //                                  const typename CrsHierBase::ordinal_type       nblks,
  //                                  const typename CrsHierBase::ordinal_type_array range,
  //                                  const typename CrsHierBase::ordinal_type_array tree) {
  //   typedef typename CrsHierBase::ordinal_type            ordinal_type;
  //   typedef typename CrsHierBase::size_type               size_type;

  //   typedef typename CrsHierBase::ordinal_type_array     ordinal_type_array;
  //   //typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
  //   //typedef typename CrsHierBase::value_type_array_ptr   value_type_array_ptr;

  //   ordinal_type_array sibling("CrsMatrixHelper::flat2hier_upper::sibling", nblks);

  //   // check the end of adjacent siblings (if not adjacent, they are separators)
  //   ordinal_type p = tree[0];
  //   for (ordinal_type i=1;i<nblks;++i) {
  //     const ordinal_type j = tree[i];
  //     if (p != j) {
  //       p = j;
  //       sibling[i-1] = -1;
  //     }
  //   }
  //   sibling[nblks-1] = -1;

  //   size_type nnz = 0;

  //   // count nnz and nnz in rows for the upper triangular hier matrix
  //   for (ordinal_type i=0;i<nblks;++i) {                  // search for all rows
  //     for (ordinal_type j=i;j != -1;j=tree[j]) {          // move up
  //       ordinal_type k=j;
  //       do {
  //         ++nnz;
  //       } while (sibling[k++] != -1);
  //     }
  //   }

  //   // create upper triangular block matrix
  //   hier.createInternalArrays(nblks, nblks, nnz);

  //   nnz = 0;
  //   for (ordinal_type i=0;i<nblks;++i) {
  //     hier._ap[i] = nnz;
  //     for (ordinal_type j=i;j != -1;j=tree[j]) {
  //       ordinal_type k=j;
  //       do {
  //         hier._aj[nnz] = k;
  //         hier._ax[nnz].setView( flat, range[i], (range[i+1] - range[i]),
  //                               /**/   range[k], (range[k+1] - range[k]));

  //         // this checking might more expensive
  //         if (hier._ax[nnz].hasNumNonZeros())
  //           ++nnz;
  //       } while (sibling[k++] != -1);
  //     }
  //   }
  //   hier._ap[nblks] = nnz;
  //   hier._nnz = nnz;

  //   return 0;
  // }

  template<typename CrsFlatBase,
           typename CrsHierBase,
           typename HostOrdinalTypeArray >
  int
  CrsMatrixHelper::flat2hier_lower(CrsFlatBase &flat,
                                   CrsHierBase &hier,
                                   const typename CrsHierBase::ordinal_type       nblks,
                                   const HostOrdinalTypeArray range,
                                   const HostOrdinalTypeArray tree) {
    ERROR(MSG_NOT_YET_IMPLEMENTED);

    // typedef typename CrsHierBase::ordinal_type           ordinal_type;
    // typedef typename CrsHierBase::size_type              size_type;

    // typedef typename CrsHierBase::ordinal_type_array     ordinal_type_array;
    // //typedef typename CrsHierBase::ordinal_type_array_ptr ordinal_type_array_ptr;
    // //typedef typename CrsHierBase::value_type_array_ptr   value_type_array_ptr;

    // ordinal_type_array tmp = ordinal_type_array("flat2hier:tmp", nblks+1);
    // size_type nnz = 0;

    // // count nnz and nnz in rows for lower triangular matrix
    // for (ordinal_type i=0;i<nblks;++i)
    //   for (ordinal_type j=i;j != -1;++nnz) {
    //     ++tmp[j];
    //     j = tree[j];
    //   }

    // // count nnz and nnz in rows for lower triangular matrix
    // hier.createInternalArrays(nblks, nblks, nnz);
    // for (ordinal_type i=1;i<(nblks+1);++i)
    //   hier._ap[i] = hier._ap[i-1] + tmp[i-1];

    // for (ordinal_type i=0;i<(nblks+1);++i)
    //   tmp[i] = hier._ap[i];

    // for (ordinal_type i=0;i<nblks;++i)
    //   for (ordinal_type j=i;j != -1;j=tree[j]) {
    //     hier._aj[tmp[j]] = i;
    //     hier._ax[tmp[j]].setView( flat, range[j], (range[j+1] - range[j]),
    //                              /**/   range[i], (range[i+1] - range[i]));
    //     ++tmp[j];
    //   }

    return 0;
  }

}


#endif

