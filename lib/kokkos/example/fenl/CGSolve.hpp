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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_CG_SOLVE
#define KOKKOS_EXAMPLE_CG_SOLVE

#include <cmath>
#include <limits>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include <WrapMPI.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< typename ValueType , class Space >
struct CrsMatrix {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  typedef Kokkos::StaticCrsGraph< unsigned , Space , void , unsigned , void >  StaticCrsGraphType ;
#else
  typedef Kokkos::StaticCrsGraph< unsigned , Space , void , void , unsigned >  StaticCrsGraphType ;
#endif
  typedef View< ValueType * , Space > coeff_type ;

  StaticCrsGraphType  graph ;
  coeff_type          coeff ;

  CrsMatrix() : graph(), coeff() {}

  CrsMatrix( const StaticCrsGraphType & arg_graph )
    : graph( arg_graph )
    , coeff( "crs_matrix_coeff" , arg_graph.entries.extent(0) )
    {}
};

template< typename MScalar 
        , typename VScalar
        , class Space >
struct Multiply {

  const Example::CrsMatrix< MScalar , Space >    m_A ;
  const Kokkos::View< const VScalar * , Space > m_x ;
  const Kokkos::View<       VScalar * , Space > m_y ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int iRow ) const
    {
      const int iEntryBegin = m_A.graph.row_map[iRow];
      const int iEntryEnd   = m_A.graph.row_map[iRow+1];

      double sum = 0 ;

      for ( int iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.coeff(iEntry) * m_x( m_A.graph.entries(iEntry) );
      }

      m_y(iRow) = sum ;
    }

  Multiply( const View<       VScalar * , Space > & y 
          , const CrsMatrix< MScalar , Space >    & A 
          , const View< const VScalar * , Space > & x 
          )
  : m_A( A ), m_x( x ), m_y( y )
  {}
};

template< typename MScalar
        , typename VScalar
        , class Space >
inline
void multiply( const int nrow
             , const Kokkos::View< VScalar * , Space >    & y
             , const Example::CrsMatrix< MScalar , Space > & A
             , const Kokkos::View< VScalar * , Space >    & x
             )
{
  Kokkos::parallel_for( Kokkos::RangePolicy<Space>(0,nrow), Multiply<MScalar,VScalar,Space>( y , A , x ) );
}

template< typename ValueType , class Space >
struct WAXPBY {
  const Kokkos::View< const ValueType * , Space >  m_x ;
  const Kokkos::View< const ValueType * , Space >  m_y ;
  const Kokkos::View<       ValueType * , Space >  m_w ;
  const double m_alpha ;
  const double m_beta ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i ) const
    { m_w(i) = m_alpha * m_x(i) + m_beta * m_y(i); }

  WAXPBY( const View< ValueType * , Space >  & arg_w
        , const double arg_alpha
        , const View< ValueType * , Space >  & arg_x
        , const double arg_beta
        , const View< ValueType * , Space >  & arg_y
        )
    : m_x( arg_x )
    , m_y( arg_y )
    , m_w( arg_w )
    , m_alpha( arg_alpha )
    , m_beta( arg_beta )
    {}
};

template< typename VScalar , class Space >
void waxpby( const int n
           , const Kokkos::View< VScalar * , Space > & arg_w
           , const double                      arg_alpha
           , const Kokkos::View< VScalar * , Space > & arg_x
           , const double                      arg_beta
           , const Kokkos::View< VScalar * , Space > & arg_y
           )
{
  Kokkos::parallel_for( Kokkos::RangePolicy<Space>(0,n), WAXPBY<VScalar,Space>(arg_w,arg_alpha,arg_x,arg_beta,arg_y) );
}

template< typename VScalar , class Space >
struct Dot {
  typedef double value_type ;

  const Kokkos::View< const VScalar * , Space >  m_x ;
  const Kokkos::View< const VScalar * , Space >  m_y ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i , value_type & update ) const
    { update += m_x(i) * m_y(i); }

  Dot( const Kokkos::View< VScalar * , Space >  & arg_x
     , const Kokkos::View< VScalar * , Space >  & arg_y
     )
    : m_x(arg_x), m_y(arg_y) {}
};

template< typename VScalar , class Space >
double dot( const int n
          , const Kokkos::View< VScalar * , Space > & arg_x
          , const Kokkos::View< VScalar * , Space > & arg_y
          )
{
  double result = 0 ;
  Kokkos::parallel_reduce( Kokkos::RangePolicy<Space>(0,n) , Dot<VScalar,Space>( arg_x , arg_y ) , result );
  return result ;
}

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

struct CGSolveResult {
  size_t  iteration ;
  double  iter_time ;
  double  matvec_time ;
  double  norm_res ;
};

template< class ImportType
        , typename MScalar
        , typename VScalar
        , class Space
        >
inline
void cgsolve( const ImportType & import
            , const CrsMatrix< MScalar , Space >      & A
            , const Kokkos::View< VScalar * , Space > & b
            , const Kokkos::View< VScalar * , Space > & x
            , const size_t  maximum_iteration = 200
            , const double  tolerance = std::numeric_limits<double>::epsilon()
            , CGSolveResult * result = 0
            )
{
  typedef View< VScalar * , Space >  VectorType ;

  const size_t count_owned = import.count_owned ;
  const size_t count_total = import.count_owned + import.count_receive;

  size_t  iteration = 0 ;
  double  iter_time = 0 ;
  double  matvec_time = 0 ;
  double  norm_res = 0 ;

  // Need input vector to matvec to be owned + received
  VectorType pAll ( "cg::p" , count_total );

  VectorType p = Kokkos::subview( pAll , std::pair<size_t,size_t>(0,count_owned) );
  VectorType r ( "cg::r" , count_owned );
  VectorType Ap( "cg::Ap", count_owned );

  /* r = b - A * x ; */

  /* p  = x       */  Kokkos::deep_copy( p , x );
  /* import p     */  import( pAll );
  /* Ap = A * p   */  multiply( count_owned , Ap , A , pAll );
  /* r = b - Ap   */  waxpby( count_owned , r , 1.0 , b , -1.0 , Ap );
  /* p  = r       */  Kokkos::deep_copy( p , r );

  double old_rdot = Kokkos::Example::all_reduce( dot( count_owned , r , r ) , import.comm );

  norm_res  = std::sqrt( old_rdot );
  iteration = 0 ;

  Kokkos::Timer wall_clock ;
  Kokkos::Timer timer;

  while ( tolerance < norm_res && iteration < maximum_iteration ) {

    /* pAp_dot = dot( p , Ap = A * p ) */

    timer.reset();
    /* import p    */  import( pAll );
    /* Ap = A * p  */  multiply( count_owned , Ap , A , pAll );
    Space::fence();
    matvec_time += timer.seconds();

    const double pAp_dot = Kokkos::Example::all_reduce( dot( count_owned , p , Ap ) , import.comm );
    const double alpha   = old_rdot / pAp_dot ;

    /* x +=  alpha * p ;  */ waxpby( count_owned , x ,  alpha, p  , 1.0 , x );
    /* r += -alpha * Ap ; */ waxpby( count_owned , r , -alpha, Ap , 1.0 , r );

    const double r_dot = Kokkos::Example::all_reduce( dot( count_owned , r , r ) , import.comm );
    const double beta  = r_dot / old_rdot ;

    /* p = r + beta * p ; */ waxpby( count_owned , p , 1.0 , r , beta , p );

    norm_res = std::sqrt( old_rdot = r_dot );

    ++iteration ;
  }

  Space::fence();
  iter_time = wall_clock.seconds();

  if ( 0 != result ) {
    result->iteration   = iteration ;
    result->iter_time   = iter_time ;
    result->matvec_time = matvec_time ;
    result->norm_res    = norm_res ;
  }
}

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_CG_SOLVE */


