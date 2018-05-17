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

#ifndef USESCASES_LINALG_BLAS_HPP
#define USESCASES_LINALG_BLAS_HPP

#include <cmath>
#include <utility>
#include <ParallelComm.hpp>
#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class Scalar , class Layout , class DeviceType > struct Dot ;

template< class Scalar , class Layout , class DeviceType > struct Dot1 ;

template< typename ScalarA ,
          typename ScalarY ,
          class Layout , class Device >
struct Scale ;

template< typename ScalarA ,
          typename ScalarY ,
          class Layout , class Device >
struct Fill ;

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarY ,
          class Layout , class Device >
struct AXPY ;

template< typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          class Layout , class Device >
struct XPBY ;

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          typename ScalarW ,
          class Layout , class Device >
struct WAXPBY ;

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_MPI )

template< typename ScalarX /* Allow mix of const and non-const */ ,
          typename ScalarY /* Allow mix of const and non-const */ ,
          class L , class D ,
          class MX /* Allow any management type */ ,
          class MY /* Allow any management type */ >
inline
double dot( const size_t n ,
            const View< ScalarX * , L , D , MX > & x ,
            const View< ScalarY * , L , D , MY > & y ,
            comm::Machine machine )
{
  double global_result = 0 ;
  double local_result = 0 ;

  Impl::Dot< ScalarX , L , D >( n , x , y , local_result );

  MPI_Allreduce( & local_result , & global_result , 1 ,
                 MPI_DOUBLE , MPI_SUM , machine.mpi_comm );

  return global_result ;
}

#else

template< typename ScalarX /* Allow mix of const and non-const */ ,
          typename ScalarY /* Allow mix of const and non-const */ ,
          class L , class D ,
          class MX /* Allow any management type */ ,
          class MY /* Allow any management type */ >
inline
double dot( const size_t n ,
            const View< ScalarX * , L , D , MX > & x ,
            const View< ScalarY * , L , D , MY > & y ,
            comm::Machine )
{
  double global_result = 0 ;

  Impl::Dot< ScalarX , L , D >( n , x , y , global_result );

  return global_result ;
}

#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_MPI )

template< typename ScalarX /* Allow mix of const and non-const */ ,
          class L , class D ,
          class MX /* Allow any management type */ >
inline
double dot( const size_t n ,
            const View< ScalarX * , L , D , MX > & x ,
            comm::Machine machine )
{
  double global_result = 0 ;
  double local_result = 0 ;

  Impl::Dot1< ScalarX , L , D >( n , x , local_result );

  MPI_Allreduce( & local_result , & global_result , 1 ,
                 MPI_DOUBLE , MPI_SUM , machine.mpi_comm );

  return global_result ;
}

#else

template< typename ScalarX /* Allow mix of const and non-const */ ,
          class L , class D ,
          class MX /* Allow any management type */ >
inline
double dot( const size_t n ,
            const View< ScalarX * , L , D , MX > & x ,
            comm::Machine )
{
  double global_result = 0 ;

  Impl::Dot1< ScalarX , L , D >( n , x , global_result );

  return global_result ;
}

#endif

//----------------------------------------------------------------------------

template< typename ScalarX /* Allow mix of const and non-const */ ,
          class L , class D ,
          class MX /* Allow any management type */ >
inline
double norm2( const size_t n ,
              const View< ScalarX * , L , D , MX > & x ,
              comm::Machine machine )
{
  return std::sqrt( dot( n , x , machine ) );
}

//----------------------------------------------------------------------------

template< typename ScalarA ,
          typename ScalarX ,
          class L ,
          class D ,
          class MX >
void scale( const size_t n ,
            const ScalarA & alpha ,
            const View< ScalarX * , L , D , MX > & x )
{
  Impl::Scale< ScalarA , ScalarX , L , D >( n , alpha , x );
}

template< typename ScalarA ,
          typename ScalarX ,
          class L ,
          class D ,
          class MX >
void fill( const size_t n ,
           const ScalarA & alpha ,
           const View< ScalarX * , L , D , MX > & x )
{
  Impl::Fill< ScalarA , ScalarX , L , D >( n , alpha , x );
}

//----------------------------------------------------------------------------

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarY ,
          class L ,
          class D ,
          class MX ,
          class MY >
void axpy( const size_t n ,
           const ScalarA & alpha ,
           const View< ScalarX *, L , D , MX > & x ,
           const View< ScalarY *, L , D , MY > & y )
{
  Impl::AXPY< ScalarA, ScalarX, ScalarY , L , D >( n, alpha, x, y );
}

//----------------------------------------------------------------------------

template< typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          class L ,
          class D ,
          class MX ,
          class MY >
void xpby( const size_t n ,
           const View< ScalarX *, L , D , MX > & x ,
           const ScalarB & beta ,
           const View< ScalarY *, L , D , MY > & y )
{
  Impl::XPBY< ScalarX, ScalarB, ScalarY , L , D >( n, x, beta, y );
}

//----------------------------------------------------------------------------
// w = alpha * x + beta * y

template< typename ScalarA ,
          typename ScalarX ,
          typename ScalarB ,
          typename ScalarY ,
          typename ScalarW ,
          class L , class D ,
          class MX , class MY , class MW >
void waxpby( const size_t n ,
             const ScalarA & alpha ,
             const View< ScalarX * , L , D , MX > & x ,
             const ScalarB & beta ,
             const View< ScalarY * , L , D , MY > & y ,
             const View< ScalarW * , L , D , MW > & w )
{
  Impl::WAXPBY<ScalarA,ScalarX,ScalarB,ScalarY,ScalarW,L,D>
    ( n , alpha , x , beta , y , w );
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Scalar , class L , class D >
struct Dot
{
private:

  typedef View< const Scalar*, L, D, MemoryUnmanaged >  vector_const_type ;

  const vector_const_type x ;
  const vector_const_type y ;

public:

  typedef typename vector_const_type::execution_space  execution_space ; // Manycore device
  typedef double      value_type ;  // Reduction value

  template< class ArgX , class ArgY >
  inline
  Dot( const size_t n , const ArgX & arg_x , const ArgY & arg_y , double & result )
    : x( arg_x ), y( arg_y )
  {
    parallel_reduce( n , *this , result );
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += x(i) * y(i); }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source;    }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
}; // Dot

//----------------------------------------------------------------------------

template< typename Scalar , class L , class D >
struct Dot1
{
private:

  typedef View< const Scalar*, L, D , MemoryUnmanaged >  vector_const_type ;

  const vector_const_type x ;

public:

  typedef typename vector_const_type::execution_space  execution_space ; // Manycore device
  typedef double      value_type ;  // Reduction value

  template< class ArgX >
  inline
  Dot1( const size_t n , const ArgX & arg_x , double & result )
    : x( arg_x )
  {
    parallel_reduce( n , *this , result );
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i , value_type & update ) const
  { update += x(i) * x(i) ; }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & source )
  { update += source ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
  { update = 0 ; }
}; // Dot

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template < typename ScalarA ,
           typename ScalarX ,
           typename ScalarB ,
           typename ScalarY ,
           typename ScalarW ,
           class L , class D >
struct WAXPBY
{
private:

  typedef View<       ScalarW *, L , D , MemoryUnmanaged > ViewW ;
  typedef View< const ScalarX *, L , D , MemoryUnmanaged > ViewX ;
  typedef View< const ScalarY *, L , D , MemoryUnmanaged > ViewY ;

  const ViewW    w ;
  const ViewX    x ;
  const ViewY    y ;
  const ScalarA  alpha ;
  const ScalarB  beta ;

public:

  typedef typename ViewW::execution_space  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType inode ) const
  {
    w(inode) = alpha * x(inode) + beta * y(inode);
  }

  template< class ArgX , class ArgY , class ArgW >
  inline
  WAXPBY( const size_t  n ,
          const ScalarA & arg_alpha ,
          const ArgX    & arg_x ,
          const ScalarB & arg_beta ,
          const ArgY    & arg_y ,
          const ArgW    & arg_w )
    : w( arg_w ), x( arg_x ), y( arg_y )
    , alpha( arg_alpha ), beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // WAXPBY

//----------------------------------------------------------------------------

template < typename ScalarB ,
           typename ScalarW ,
           class L , class D >
struct Scale
{
private:

  typedef View< ScalarW *, L , D , MemoryUnmanaged >  ViewW ;
  const ViewW    w ;
  const ScalarB  beta ;

public:

  typedef typename ViewW::execution_space  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) *= beta ; }

  template< class ArgW >
  inline
  Scale( const size_t n , const ScalarB & arg_beta , const ArgW & arg_w )
    : w( arg_w )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
};

template < typename ScalarB ,
           typename ScalarW ,
           class L , class D >
struct Fill
{
private:

  typedef View< ScalarW *, L , D , MemoryUnmanaged >  ViewW ;
  const ViewW    w ;
  const ScalarB  beta ;

public:

  typedef typename ViewW::execution_space  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) = beta ; }

  template< class ArgW >
  inline
  Fill( const size_t n , const ScalarB & arg_beta , const ArgW & arg_w )
    : w( arg_w )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
};

//----------------------------------------------------------------------------

template < typename ScalarA ,
           typename ScalarX ,
           typename ScalarW ,
           class L , class D >
struct AXPY
{
private:

  typedef View<       ScalarW *, L , D , MemoryUnmanaged >  ViewW ;
  typedef View< const ScalarX *, L , D , MemoryUnmanaged >  ViewX ;

  const ViewW    w ;
  const ViewX    x ;
  const ScalarA  alpha ;

public:

  typedef typename ViewW::execution_space  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) += alpha * x(i); }

  template< class ArgX , class ArgW >
  inline
  AXPY( const size_t  n ,
        const ScalarA & arg_alpha ,
        const ArgX    & arg_x ,
        const ArgW    & arg_w )
    : w( arg_w ), x( arg_x )
    , alpha( arg_alpha )
  {
    parallel_for( n , *this );
  }
}; // AXPY

template< typename ScalarX ,
          typename ScalarB ,
          typename ScalarW ,
          class L , class D >
struct XPBY
{
private:

  typedef View<       ScalarW *, L , D , MemoryUnmanaged >  ViewW ;
  typedef View< const ScalarX *, L , D , MemoryUnmanaged >  ViewX ;

  const ViewW    w ;
  const ViewX    x ;
  const ScalarB  beta ;

public:

  typedef typename ViewW::execution_space  execution_space ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void operator()( const iType & i ) const
  { w(i) = x(i) + beta * w(i); }

  template< class ArgX , class ArgW >
  inline
  XPBY( const size_t  n ,
        const ArgX    & arg_x ,
        const ScalarB & arg_beta ,
        const ArgW    & arg_w )
    : w( arg_w ), x( arg_x )
    , beta( arg_beta )
  {
    parallel_for( n , *this );
  }
}; // XPBY

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef USESCASES_LINALG_BLAS_HPP */


