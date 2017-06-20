/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef SPARSELINEARSYSTEM_HPP
#define SPARSELINEARSYSTEM_HPP

#include <cmath>
#include <impl/Kokkos_Timer.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

#include <LinAlgBLAS.hpp>

namespace Kokkos {

template< typename ScalarType , class Device >
struct CrsMatrix {
  typedef Device      execution_space ;
  typedef ScalarType  value_type ;

  typedef StaticCrsGraph< int , execution_space , void , int >  graph_type ;
  typedef View< value_type* , execution_space >   coefficients_type ;

  graph_type         graph ;
  coefficients_type  coefficients ;
};

//----------------------------------------------------------------------------

namespace Impl {

template< class Matrix , class OutputVector , class InputVector >
struct Multiply ;

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename AScalarType ,
          typename VScalarType ,
          class DeviceType >
struct Multiply< CrsMatrix<AScalarType,DeviceType> ,
                 View<VScalarType*,DeviceType > ,
                 View<VScalarType*,DeviceType > >
{
  typedef DeviceType                       execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef View<       VScalarType*, execution_space, MemoryUnmanaged >  vector_type ;
  typedef View< const VScalarType*, execution_space, MemoryUnmanaged >  vector_const_type ;

  typedef CrsMatrix< AScalarType , execution_space >    matrix_type ;

private:

  matrix_type        m_A ;
  vector_const_type  m_x ;
  vector_type        m_y ;

public:

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    double sum = 0 ;

#if defined( __INTEL_COMPILER )
#pragma simd reduction(+:sum)
#pragma ivdep
    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      sum += m_A.coefficients(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }
#else
    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      sum += m_A.coefficients(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }
#endif

    m_y(iRow) = sum ;
  }

  Multiply( const matrix_type & A ,
            const size_type nrow ,
            const size_type , // ncol ,
            const vector_type & x ,
            const vector_type & y )
    : m_A( A ), m_x( x ), m_y( y )
  {
    parallel_for( nrow , *this );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename AScalarType ,
          typename VScalarType ,
          class Device >
class Operator {
  typedef CrsMatrix<AScalarType,Device>  matrix_type ;
  typedef View<VScalarType*,Device>     vector_type ;

private:
  const CrsMatrix<AScalarType,Device> A ;

  ParallelDataMap                                         data_map ;
  AsyncExchange< VScalarType , Device , ParallelDataMap > exchange ;

public:

  Operator( const ParallelDataMap                  & arg_data_map ,
            const CrsMatrix<AScalarType,Device>    & arg_A )
    : A( arg_A )
    , data_map( arg_data_map )
    , exchange( arg_data_map , 1 )
    {}

  void apply( const View<VScalarType*,Device>  & x ,
              const View<VScalarType*,Device>  & y )
  {
    // Gather off-processor data for 'x'

    PackArray< vector_type >::pack( exchange.buffer() ,
                                    data_map.count_interior ,
                                    data_map.count_send , x );

    exchange.setup();

    // If interior & boundary matrices then could launch interior multiply

    exchange.send_receive();

    UnpackArray< vector_type >::unpack( x , exchange.buffer() ,
                                        data_map.count_owned ,
                                        data_map.count_receive );

    const typename Device::size_type nrow = data_map.count_owned ;
    const typename Device::size_type ncol = data_map.count_owned +
                                            data_map.count_receive ;

    Impl::Multiply<matrix_type,vector_type,vector_type>( A, nrow, ncol, x, y);
  }
};

//----------------------------------------------------------------------------

template< typename AScalarType , typename VScalarType , class Device >
void cgsolve(
  const ParallelDataMap                 data_map ,
  const CrsMatrix<AScalarType,Device>   A ,
  const View<VScalarType*,Device> b ,
  const View<VScalarType*,Device> x ,
  size_t & iteration ,
  double & normr ,
  double & iter_time ,
  const size_t maximum_iteration = 200 ,
  const double tolerance = std::numeric_limits<VScalarType>::epsilon() )
{
  typedef View<VScalarType*,Device> vector_type ;
  //typedef View<VScalarType,  Device> value_type ; // unused

  const size_t count_owned = data_map.count_owned ;
  const size_t count_total = data_map.count_owned + data_map.count_receive ;

  Operator<AScalarType,VScalarType,Device> matrix_operator( data_map , A );

  // Need input vector to matvec to be owned + received
  vector_type pAll ( "cg::p" , count_total );

  vector_type p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_owned) );
  vector_type r ( "cg::r" , count_owned );
  vector_type Ap( "cg::Ap", count_owned );

  /* r = b - A * x ; */

  /* p  = x      */ deep_copy( p , x );
  /* Ap = A * p  */ matrix_operator.apply( pAll , Ap );
  /* r  = b - Ap */ waxpby( count_owned , 1.0 , b , -1.0 , Ap , r );
  /* p  = r      */ deep_copy( p , r );

  double old_rdot = dot( count_owned , r , data_map.machine );

  normr     = std::sqrt( old_rdot );
  iteration = 0 ;

  Kokkos::Timer wall_clock ;

  while ( tolerance < normr && iteration < maximum_iteration ) {

    /* pAp_dot = dot( p , Ap = A * p ) */

    /* Ap = A * p  */ matrix_operator.apply( pAll , Ap );

    const double pAp_dot = dot( count_owned , p , Ap , data_map.machine );
    const double alpha   = old_rdot / pAp_dot ;

    /* x += alpha * p ;  */ axpy( count_owned,  alpha, p , x );
    /* r -= alpha * Ap ; */ axpy( count_owned, -alpha, Ap, r );

    const double r_dot = dot( count_owned , r , data_map.machine );
    const double beta  = r_dot / old_rdot ;

    /* p = r + beta * p ; */ xpby( count_owned , r , beta , p );

    normr = std::sqrt( old_rdot = r_dot );
    ++iteration ;
  }

  iter_time = wall_clock.seconds();
}

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_CUDA )

#if ( CUDA_VERSION < 6000 )
#pragma message "cusparse_v2.h"
#include <cusparse_v2.h>
#else
#pragma message "cusparse.h"
#include <cusparse.h>
#endif

namespace Kokkos {
namespace Impl {

struct CudaSparseSingleton {
  cusparseHandle_t   handle;
  cusparseMatDescr_t descra;

  CudaSparseSingleton()
  {
    cusparseCreate( & handle );
    cusparseCreateMatDescr( & descra );
    cusparseSetMatType(       descra , CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase(  descra , CUSPARSE_INDEX_BASE_ZERO );
  }

  static CudaSparseSingleton & singleton();

};

template<>
struct Multiply< CrsMatrix<double,Cuda> ,
                 View<double*,Cuda > ,
                 View<double*,Cuda > >
{
  typedef Cuda                                      execution_space ;
  typedef execution_space::size_type                    size_type ;
  typedef double                                    scalar_type ;
  typedef View< scalar_type* , execution_space >        vector_type ;
  typedef CrsMatrix< scalar_type , execution_space >    matrix_type ;

public:

  Multiply( const matrix_type & A ,
            const size_type nrow ,
            const size_type ncol ,
            const vector_type & x ,
            const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const scalar_type alpha = 1 , beta = 0 ;

    cusparseStatus_t status =
      cusparseDcsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      nrow , ncol , A.coefficients.dimension_0() ,
                      &alpha ,
                      s.descra ,
                      A.coefficients.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() ,
                      &beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};


template<>
struct Multiply< CrsMatrix<float,Cuda> ,
                 View<float*,Cuda > ,
                 View<float*,Cuda > >
{
  typedef Cuda                                      execution_space ;
  typedef execution_space::size_type                    size_type ;
  typedef float                                     scalar_type ;
  typedef View< scalar_type* , execution_space >        vector_type ;
  typedef CrsMatrix< scalar_type , execution_space >    matrix_type ;

public:

  Multiply( const matrix_type & A ,
            const size_type nrow ,
            const size_type ncol ,
            const vector_type & x ,
            const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const scalar_type alpha = 1 , beta = 0 ;

    cusparseStatus_t status =
      cusparseScsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      nrow , ncol , A.coefficients.dimension_0() ,
                      &alpha ,
                      s.descra ,
                      A.coefficients.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() ,
                      &beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef SPARSELINEARSYSTEM_HPP */

