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

#ifndef KOKKOS_IMPL_REDUCER_HPP
#define KOKKOS_IMPL_REDUCER_HPP

#include <impl/Kokkos_Traits.hpp>

//----------------------------------------------------------------------------
/*  Reducer abstraction:
 *  1) Provides 'join' operation
 *  2) Provides 'init' operation
 *  3) Provides 'copy' operation
 *  4) Optionally provides result value in a memory space
 *
 *  Created from:
 *  1) Functor::operator()( destination , source )
 *  2) Functor::{ join , init )
 */
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename value_type >
struct ReduceSum
{
  KOKKOS_INLINE_FUNCTION static
  void copy( value_type & dest
           , value_type const & src ) noexcept
    { dest = src ; }

  KOKKOS_INLINE_FUNCTION static
  void init( value_type & dest ) noexcept
    { new( &dest ) value_type(); }

  KOKKOS_INLINE_FUNCTION static
  void join( value_type volatile & dest
           , value_type const volatile & src ) noexcept
    { dest += src ; }

  KOKKOS_INLINE_FUNCTION static
  void join( value_type & dest
           , value_type const & src ) noexcept
    { dest += src ; }
};

template< typename T
        , class ReduceOp = ReduceSum< T >
        , typename MemorySpace = void >
struct Reducer
  : private ReduceOp
  , private integral_nonzero_constant
    < int , ( std::rank<T>::value == 1 ? std::extent<T>::value : 1 )>
{
private:

  // Determine if T is simple array

  enum : int { rank = std::rank<T>::value };

  static_assert( rank <= 1 , "Kokkos::Impl::Reducer type is at most rank-one" );

  using length_t =
    integral_nonzero_constant<int,( rank == 1 ? std::extent<T>::value : 1 )> ;

public:

  using reducer        = Reducer ;
  using memory_space   = MemorySpace ;
  using value_type     = typename std::remove_extent<T>::type ;
  using reference_type =
    typename std::conditional< ( rank != 0 )
                             , value_type *
                             , value_type &
                             >::type ;
private:

  //--------------------------------------------------------------------------
  // Determine what functions 'ReduceOp' provides:
  //   copy( destination , source )
  //   init( destination )
  //
  //   operator()( destination , source )
  //   join( destination , source )
  //
  // Provide defaults for missing optional operations

  template< class R , typename = void>
  struct COPY {
    KOKKOS_INLINE_FUNCTION static
    void copy( R const &
             , value_type * dst
             , value_type const * src ) { *dst = *src ; }
  };

  template< class R >
  struct COPY< R , decltype( ((R*)0)->copy( *((value_type*)0)
                                          , *((value_type const *)0) ) ) >
  {
    KOKKOS_INLINE_FUNCTION static
    void copy( R const & r
             , value_type * dst
             , value_type const * src ) { r.copy( *dst , *src ); }
  };

  template< class R , typename = void >
  struct INIT {
    KOKKOS_INLINE_FUNCTION static
    void init( R const & , value_type * dst ) { new(dst) value_type(); }
  };

  template< class R >
  struct INIT< R , decltype( ((R*)0)->init( *((value_type*)0 ) ) ) >
  {
    KOKKOS_INLINE_FUNCTION static
    void init( R const & r , value_type * dst ) { r.init( *dst ); }
  };

  template< class R , typename V , typename = void > struct JOIN
    {
      // If no join function then try operator()
      KOKKOS_INLINE_FUNCTION static
      void join( R const & r , V * dst , V const * src )
        { r.operator()(*dst,*src); }
    };

  template< class R , typename V >
  struct JOIN< R , V , decltype( ((R*)0)->join ( *((V *)0) , *((V const *)0) ) ) >
    {
      // If has join function use it
      KOKKOS_INLINE_FUNCTION static
      void join( R const & r , V * dst , V const * src )
        { r.join(*dst,*src); }
    };

  //--------------------------------------------------------------------------

  value_type * const m_result ;

  template< int Rank >
  KOKKOS_INLINE_FUNCTION
  static constexpr
  typename std::enable_if< ( 0 != Rank ) , reference_type >::type
  ref( value_type * p ) noexcept { return p ; }

  template< int Rank >
  KOKKOS_INLINE_FUNCTION
  static constexpr
  typename std::enable_if< ( 0 == Rank ) , reference_type >::type
  ref( value_type * p ) noexcept { return *p ; }

public:

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr int length() const noexcept
     { return length_t::value ; }

  KOKKOS_INLINE_FUNCTION
  value_type * data() const noexcept
    { return m_result ; }

  KOKKOS_INLINE_FUNCTION
  reference_type reference() const noexcept
    { return Reducer::template ref< rank >( m_result ); }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void copy( value_type * const dest
           , value_type const * const src ) const noexcept
    {
      for ( int i = 0 ; i < length() ; ++i ) {
        Reducer::template COPY<ReduceOp>::copy( (ReduceOp &) *this , dest + i , src + i );
      }
    }

  KOKKOS_INLINE_FUNCTION
  void init( value_type * dest ) const noexcept
    {
      for ( int i = 0 ; i < length() ; ++i ) {
        Reducer::template INIT<ReduceOp>::init( (ReduceOp &) *this , dest + i );
      }
    }

  KOKKOS_INLINE_FUNCTION
  void join( value_type * const dest
           , value_type const * const src ) const noexcept
    {
      for ( int i = 0 ; i < length() ; ++i ) {
        Reducer::template JOIN<ReduceOp,value_type>::join( (ReduceOp &) *this , dest + i , src + i );
      }
    }

  KOKKOS_INLINE_FUNCTION
  void join( value_type volatile * const dest
           , value_type volatile const * const src ) const noexcept
    {
      for ( int i = 0 ; i < length() ; ++i ) {
        Reducer::template JOIN<ReduceOp,value_type volatile>::join( (ReduceOp &) *this , dest + i , src + i );
      }
    }

  //--------------------------------------------------------------------------

  template< typename ArgT >
  KOKKOS_INLINE_FUNCTION explicit
  constexpr Reducer
    ( ArgT * arg_value
    , typename std::enable_if
        < std::is_same<ArgT,value_type>::value &&
          std::is_default_constructible< ReduceOp >::value
        , int >::type arg_length = 1
    ) noexcept
    : ReduceOp(), length_t( arg_length ), m_result( arg_value ) {}

  KOKKOS_INLINE_FUNCTION explicit
  constexpr Reducer( ReduceOp const & arg_op
                   , value_type     * arg_value = 0
                   , int arg_length = 1 ) noexcept
    : ReduceOp( arg_op ), length_t( arg_length ), m_result( arg_value ) {}

  KOKKOS_INLINE_FUNCTION explicit
  constexpr Reducer( ReduceOp      && arg_op
                   , value_type     * arg_value = 0
                   , int arg_length = 1 ) noexcept
    : ReduceOp( arg_op ), length_t( arg_length ), m_result( arg_value ) {}

  Reducer( Reducer const & ) = default ;
  Reducer( Reducer && ) = default ;
  Reducer & operator = ( Reducer const & ) = default ;
  Reducer & operator = ( Reducer && ) = default ;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< typename ValueType >
constexpr
Impl::Reducer< ValueType , Impl::ReduceSum< ValueType > >
Sum( ValueType & arg_value )
{
  static_assert( std::is_trivial<ValueType>::value
               , "Kokkos reducer requires trivial value type" );
  return Impl::Reducer< ValueType , Impl::ReduceSum< ValueType > >( & arg_value );
}

template< typename ValueType >
constexpr
Impl::Reducer< ValueType[] , Impl::ReduceSum< ValueType > >
Sum( ValueType * arg_value , int arg_length )
{
  static_assert( std::is_trivial<ValueType>::value
               , "Kokkos reducer requires trivial value type" );
  return Impl::Reducer< ValueType[] , Impl::ReduceSum< ValueType > >( arg_value , arg_length );
}

//----------------------------------------------------------------------------

template< typename ValueType , class JoinType >
Impl::Reducer< ValueType , JoinType >
reducer( ValueType & value , JoinType const & lambda )
{
  return Impl::Reducer< ValueType , JoinType >( lambda , & value );
}

} // namespace Kokkos

#endif /* #ifndef KOKKOS_IMPL_REDUCER_HPP */

