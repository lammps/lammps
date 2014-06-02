/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

#ifndef KOKKOS_VIEWSUPPORT_HPP
#define KOKKOS_VIEWSUPPORT_HPP

#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Evaluate if LHS = RHS view assignment is allowed. */
template< class ViewLHS , class ViewRHS >
struct ViewAssignable
{
  // Same memory space.
  // Same value type.
  // Compatible 'const' qualifier
  // Cannot assign managed = unmannaged
  enum { assignable_value =
    ( is_same< typename ViewLHS::value_type ,
               typename ViewRHS::value_type >::value
      ||
      is_same< typename ViewLHS::value_type ,
               typename ViewRHS::const_value_type >::value )
    &&
    is_same< typename ViewLHS::memory_space ,
             typename ViewRHS::memory_space >::value
    &&
    ( ! ( ViewLHS::is_managed && ! ViewRHS::is_managed ) )
  };

  enum { assignable_shape =
    // Compatible shape and matching layout:
    ( ShapeCompatible< typename ViewLHS::shape_type ,
                       typename ViewRHS::shape_type >::value
      &&
      is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value )
    ||
    // Matching layout, same rank, and LHS dynamic rank
    ( is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value
      &&
      int(ViewLHS::rank) == int(ViewRHS::rank)
      &&
      int(ViewLHS::rank) == int(ViewLHS::rank_dynamic) )
    ||
    // Both rank-0, any shape and layout
    ( int(ViewLHS::rank) == 0 && int(ViewRHS::rank) == 0 )
    ||
    // Both rank-1 and LHS is dynamic rank-1, any shape and layout
    ( int(ViewLHS::rank) == 1 && int(ViewRHS::rank) == 1 &&
      int(ViewLHS::rank_dynamic) == 1 )
    };

  enum { value = assignable_value && assignable_shape };
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  View tracking increment/decrement only happens when
 *          view memory is managed and executing in the host space.
 */
template< class ViewTraits , class Enable = void >
struct ViewTracking {
  KOKKOS_INLINE_FUNCTION void increment( const void * ) const {}
  KOKKOS_INLINE_FUNCTION void decrement( const void * ) const {}

  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const ViewTracking & ) { return *this ; }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const ViewTracking<T> & ) { return *this ; }

  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const bool ) { return *this ; }

  KOKKOS_INLINE_FUNCTION
  operator bool() const { return false ; }
};

template< class ViewTraits >
struct ViewTracking< ViewTraits , typename enable_if< ViewTraits::is_managed >::type >
{
private:

  enum { is_host_space = is_same< HostSpace , ExecutionSpace >::value };

  bool m_flag ;

  struct NoType {};

public:

  typedef typename ViewTraits::memory_space memory_space ;

  template< class T >
  KOKKOS_INLINE_FUNCTION
  void increment( const T * ptr
                , typename enable_if<( ! is_same<T,NoType>::value && is_host_space )>::type * = 0 ) const
    { if ( m_flag ) memory_space::increment( ptr ); }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  void increment( const T *
                , typename enable_if<( ! is_same<T,NoType>::value && ! is_host_space )>::type * = 0 ) const
    {}

  template< class T >
  KOKKOS_INLINE_FUNCTION
  void decrement( const T * ptr
                , typename enable_if<( ! is_same<T,NoType>::value && is_host_space )>::type * = 0 ) const
    { if ( m_flag ) memory_space::decrement( ptr ); }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  void decrement( const T *
                , typename enable_if<( ! is_same<T,NoType>::value && ! is_host_space )>::type * = 0 ) const
    {}

  KOKKOS_INLINE_FUNCTION
  ViewTracking() : m_flag( true ) {}

  template< class T >
  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const ViewTracking & rhs ) { m_flag = rhs.m_flag ; return *this ; }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const ViewTracking<T> & rhs ) { m_flag = rhs.operator bool(); return *this ; }

  KOKKOS_INLINE_FUNCTION
  ViewTracking & operator = ( const bool rhs ) { m_flag = rhs ; return *this ; }

  KOKKOS_INLINE_FUNCTION
  operator bool() const { return m_flag ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class OutputView , class InputView  , unsigned Rank = OutputView::Rank >
struct ViewRemap
{
  typedef typename OutputView::device_type device_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;
  const InputView  input ;
  const size_type n0 ;
  const size_type n1 ;
  const size_type n2 ;
  const size_type n3 ;
  const size_type n4 ;
  const size_type n5 ;
  const size_type n6 ;
  const size_type n7 ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    , n0( std::min( (size_t)arg_out.dimension_0() , (size_t)arg_in.dimension_0() ) )
    , n1( std::min( (size_t)arg_out.dimension_1() , (size_t)arg_in.dimension_1() ) )
    , n2( std::min( (size_t)arg_out.dimension_2() , (size_t)arg_in.dimension_2() ) )
    , n3( std::min( (size_t)arg_out.dimension_3() , (size_t)arg_in.dimension_3() ) )
    , n4( std::min( (size_t)arg_out.dimension_4() , (size_t)arg_in.dimension_4() ) )
    , n5( std::min( (size_t)arg_out.dimension_5() , (size_t)arg_in.dimension_5() ) )
    , n6( std::min( (size_t)arg_out.dimension_6() , (size_t)arg_in.dimension_6() ) )
    , n7( std::min( (size_t)arg_out.dimension_7() , (size_t)arg_in.dimension_7() ) )
    {
      parallel_for( n0 , *this );
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < n7 ; ++i7 ) {
      output.at(i0,i1,i2,i3,i4,i5,i6,i7) = input.at(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< class OutputView , class InputView  >
struct ViewRemap< OutputView ,  InputView , 0 >
{
  typedef typename OutputView::value_type   value_type ;
  typedef typename OutputView::memory_space dst_space ;
  typedef typename InputView ::memory_space src_space ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
  {
    DeepCopy< dst_space , src_space >( arg_out.ptr_on_device() ,
                                       arg_in.ptr_on_device() ,
                                       sizeof(value_type) );
  }
};

//----------------------------------------------------------------------------

template< class OutputView , unsigned Rank = OutputView::Rank >
struct ViewFill
{
  typedef typename OutputView::device_type       device_type ;
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename device_type::size_type        size_type ;

  const OutputView output ;
  const_value_type input ;

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      parallel_for( output.dimension_0() , *this );
      device_type::fence();
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      output.at(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }
};

template< class OutputView >
struct ViewFill< OutputView , 0 >
{
  typedef typename OutputView::device_type       device_type ;
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::memory_space      dst_space ;

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
  {
    DeepCopy< dst_space , dst_space >( arg_out.ptr_on_device() , & arg_in ,
                                       sizeof(const_value_type) );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWSUPPORT_HPP */


