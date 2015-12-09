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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP

#include <Kokkos_Array.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template< class DataType , class V , long N , class P , class ArrayLayout >
struct ViewDataAnalysis< DataType , Kokkos::Array<V,N,P> , ArrayLayout >
{
private:

  typedef ViewArrayAnalysis<DataType> array_analysis ;

  static_assert( std::is_same<P,void>::value , "" );
  static_assert( std::is_same<typename array_analysis::non_const_value_type , Kokkos::Array<V,N,P> >::value , "" );
  static_assert( std::is_scalar<V>::value , "View of Array type must be of a scalar type" );

public:

  typedef Kokkos::Array<>  specialize ;

  typedef typename array_analysis::dimension  dimension ;

private:

  enum { is_const = std::is_same< typename array_analysis::value_type
                                , typename array_analysis::const_value_type
                                >::value };

  typedef ViewDimension< ( dimension::rank == 0 ? N : dimension::arg_N0 )
                       , ( dimension::rank == 1 ? N : dimension::arg_N1 )
                       , ( dimension::rank == 2 ? N : dimension::arg_N2 )
                       , ( dimension::rank == 3 ? N : dimension::arg_N3 )
                       , ( dimension::rank == 4 ? N : dimension::arg_N4 )
                       , ( dimension::rank == 5 ? N : dimension::arg_N5 )
                       , ( dimension::rank == 6 ? N : dimension::arg_N6 )
                       , ( dimension::rank == 7 ? N : dimension::arg_N7 )
                       > array_scalar_dimension ;

  typedef typename std::conditional< is_const , const V , V >::type  scalar_type ;
  typedef V       non_const_scalar_type ;
  typedef const V const_scalar_type ;

public:

  typedef typename array_analysis::value_type            value_type ;
  typedef typename array_analysis::const_value_type      const_value_type ;
  typedef typename array_analysis::non_const_value_type  non_const_value_type ;

  typedef typename ViewDataType<           value_type , dimension >::type  type ;
  typedef typename ViewDataType<     const_value_type , dimension >::type  const_type ;
  typedef typename ViewDataType< non_const_value_type , dimension >::type  non_const_type ;

  typedef typename ViewDataType<           scalar_type , array_scalar_dimension >::type  array_scalar_type ;
  typedef typename ViewDataType<     const_scalar_type , array_scalar_dimension >::type  const_array_scalar_type ;
  typedef typename ViewDataType< non_const_scalar_type , array_scalar_dimension >::type  non_const_array_scalar_type ;
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

/** \brief  View mapping for non-specialized data type and standard layout */
template< class Traits >
class ViewMapping< Traits , void ,
  typename std::enable_if<( std::is_same< typename Traits::specialize , Kokkos::Array<> >::value &&
                            ( std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value ||
                              std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ||
                              std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value )
                          )>::type >
{
private:

  template< class , class , typename > friend class ViewMapping ;
  template< class , bool , bool , bool , bool , bool , bool , bool , bool , class > friend struct SubviewMapping ;
  template< class , class , class , class > friend class Kokkos::Experimental::View ;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  typedef typename Traits::value_type::pointer handle_type ;

  handle_type  m_handle ;
  offset_type  m_offset ;
  size_t       m_stride ;

  typedef typename Traits::value_type::value_type scalar_type ;

  typedef Kokkos::Array< scalar_type , ~size_t(0) , Kokkos::Array<>::contiguous >  contiguous_reference ;
  typedef Kokkos::Array< scalar_type , ~size_t(0) , Kokkos::Array<>::strided >     strided_reference ;

  enum { is_contiguous_reference =
    ( Traits::rank == 0 ) || ( std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ) };

  enum { Array_N = Traits::value_type::size() };
  enum { Array_S = is_contiguous_reference ? Array_N : 1 };

  KOKKOS_INLINE_FUNCTION
  ViewMapping( const handle_type & arg_handle , const offset_type & arg_offset )
    : m_handle( arg_handle )
    , m_offset( arg_offset )
    , m_stride( is_contiguous_reference ? 0 : arg_offset.span() )
    {}

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const { return m_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return m_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return m_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return m_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return m_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return m_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return m_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return m_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  using is_regular = typename offset_type::is_regular ;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return m_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return m_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return m_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return m_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return m_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return m_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return m_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return m_offset.stride_7(); }

  //----------------------------------------
  // Range span

  /** \brief  Span of the mapped range */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return m_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return m_offset.span_is_contiguous(); }

  typedef typename std::conditional< is_contiguous_reference , contiguous_reference , strided_reference >::type  reference_type ;

  /** \brief  If data references are lvalue_reference than can query pointer to memory */
  KOKKOS_INLINE_FUNCTION constexpr typename Traits::value_type * data() const
    { return (typename Traits::value_type *) 0 ; }

  //----------------------------------------
  // The View class performs all rank and bounds checking before
  // calling these element reference methods.

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const { return reference_type( m_handle + 0 , Array_N , 0 ); }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type
  reference( const I0 & i0 ) const
    { return reference_type( m_handle + m_offset(i0) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return reference_type( m_handle + m_offset(i0,i1) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2,i3) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2,i3,i4) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2,i3,i4,i5) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2,i3,i4,i5,i6) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 ) const
    { return reference_type( m_handle + m_offset(i0,i1,i2,i3,i4,i5,i6,i7) * Array_S , Array_N , m_stride ); }

  //----------------------------------------

private:

  enum { MemorySpanMask = 8 - 1 /* Force alignment on 8 byte boundary */ };
  enum { MemorySpanSize = sizeof(typename Traits::value_type) };

public:

  /** \brief  Span, in bytes, of the referenced memory */
  KOKKOS_INLINE_FUNCTION constexpr size_t memory_span() const
    {
      return ( m_stride * sizeof(typename Traits::value_type) + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  /** \brief  Span, in bytes, of the required memory */
  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span( const std::integral_constant<bool,AllowPadding> &
                                     , const size_t N0 , const size_t N1 , const size_t N2 , const size_t N3
                                      , const size_t N4 , const size_t N5 , const size_t N6 , const size_t N7 )
    {
      typedef std::integral_constant< unsigned , AllowPadding ? MemorySpanSize : 0 >  padding ;
      return ( offset_type( padding(), N0, N1, N2, N3, N4, N5, N6, N7 ).span() * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  /** \brief  Span, in bytes, of the required memory */
  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span( const std::integral_constant<bool,AllowPadding> &
                                       , const typename Traits::array_layout & layout )
    {
      return ( offset_type( layout ).span() * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION ~ViewMapping() {}
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_handle(), m_offset(), m_stride(0) {}
  KOKKOS_INLINE_FUNCTION ViewMapping( const ViewMapping & rhs )
    : m_handle( rhs.m_handle ), m_offset( rhs.m_offset ), m_stride( rhs.m_stride ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( const ViewMapping & rhs )
    { m_handle = rhs.m_handle ; m_offset = rhs.m_offset ; m_stride = rhs.m_stride ; ; return *this ; }

  KOKKOS_INLINE_FUNCTION ViewMapping( ViewMapping && rhs )
    : m_handle( rhs.m_handle ), m_offset( rhs.m_offset ), m_stride( rhs.m_stride ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( ViewMapping && rhs )
    { m_handle = rhs.m_handle ; m_offset = rhs.m_offset ; m_stride = rhs.m_stride ; return *this ; }

  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  ViewMapping( void * ptr
             , const std::integral_constant<bool,AllowPadding> &
             , const size_t N0 , const size_t N1 , const size_t N2 , const size_t N3
             , const size_t N4 , const size_t N5 , const size_t N6 , const size_t N7 )
    : m_handle( reinterpret_cast< handle_type >( ptr ) )
    , m_offset( std::integral_constant< unsigned , AllowPadding ? sizeof(typename Traits::value_type) : 0 >()
              , N0, N1, N2, N3, N4, N5, N6, N7 )
    , m_stride( m_offset.span() )
    {}

  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  ViewMapping( void * ptr
             , const std::integral_constant<bool,AllowPadding> &
             , const typename Traits::array_layout & layout )
    : m_handle( reinterpret_cast< handle_type >( ptr ) )
    , m_offset( layout )
    , m_stride( m_offset.span() )
    {}

  //----------------------------------------
  // If the View is to construct or destroy the elements.

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const size_t i ) const
    {
      reference_type ref( m_handle + i * Array_S , Array_N , m_stride );
      for ( size_t j = 0 ; j < Array_N ; ++j ) ref[j] = 0 ;
    }

  template< class ExecSpace >
  void construct( const ExecSpace & space ) const
    {
      typedef Kokkos::RangePolicy< ExecSpace , size_t > Policy ;

      (void) Kokkos::Impl::ParallelFor< ViewMapping , Policy >( *this , Policy( 0 , m_stride ) );
      ExecSpace::fence();
    }

  template< class ExecSpace >
  void destroy( const ExecSpace & ) const {}
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Assign compatible default mappings */

template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space , typename SrcTraits::memory_space >::value
    &&
    std::is_same< typename DstTraits::specialize , Kokkos::Array<> >::value
    &&
    (
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutStride >::value
    )
    &&
    std::is_same< typename SrcTraits::specialize , Kokkos::Array<> >::value
    &&
    (
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutStride >::value
    )
  )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Experimental::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void , void >  DstType ;
  typedef ViewMapping< SrcTraits , void , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION
  static void assign( DstType & dst , const SrcType & src , const TrackType & src_track )
    {
      static_assert( std::is_same< typename DstTraits::value_type , typename SrcTraits::value_type >::value ||
                     std::is_same< typename DstTraits::value_type , typename SrcTraits::const_value_type >::value
                   , "View assignment must have same value type or const = non-const" );

      static_assert( ViewDimensionAssignable< typename DstTraits::dimension , typename SrcTraits::dimension >::value
                   , "View assignment must have compatible dimensions" );

      static_assert( std::is_same< typename DstTraits::array_layout , typename SrcTraits::array_layout >::value ||
                     std::is_same< typename DstTraits::array_layout , Kokkos::LayoutStride >::value ||
                     ( DstTraits::dimension::rank == 0 ) ||
                     ( DstTraits::dimension::rank == 1 && DstTraits::dimension::rank_dynamic == 1 )
                   , "View assignment must have compatible layout or have rank <= 1" );

      typedef typename DstType::offset_type  dst_offset_type ;

      dst.m_offset = dst_offset_type( src.m_offset );
      dst.m_handle = src.m_handle ;
      dst.m_stride = src.m_stride ;
    }
};

/** \brief Assign Array to non-Array */

template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space , typename SrcTraits::memory_space >::value
    &&
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    (
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutStride >::value
    )
    &&
    std::is_same< typename SrcTraits::specialize , Kokkos::Array<> >::value
    &&
    (
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutStride >::value
    )
  )>::type >
{
public:

  // Can only convert to View::array_type

  enum { is_assignable = std::is_same< typename DstTraits::data_type ,    typename SrcTraits::array_scalar_type >::value &&
                         std::is_same< typename DstTraits::array_layout , typename SrcTraits::array_layout >::value };

  typedef Kokkos::Experimental::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void , void >  DstType ;
  typedef ViewMapping< SrcTraits , void , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION
  static void assign( DstType & dst , const SrcType & src , const TrackType & src_track )
    {
      static_assert( is_assignable , "Can only convert to array_type" );

      typedef typename DstType::offset_type  dst_offset_type ;

      // Array dimension becomes the last dimension.
      // Arguments beyond the destination rank are ignored.
      if ( src.span_is_contiguous() ) { // not padded
        dst.m_offset = dst_offset_type( std::integral_constant<unsigned,0>()
                                      , ( 1 < SrcType::Rank ? src.dimension_1() : SrcTraits::value_type::size() )
                                      , ( 2 < SrcType::Rank ? src.dimension_2() : SrcTraits::value_type::size() )
                                      , ( 3 < SrcType::Rank ? src.dimension_3() : SrcTraits::value_type::size() )
                                      , ( 4 < SrcType::Rank ? src.dimension_4() : SrcTraits::value_type::size() )
                                      , ( 5 < SrcType::Rank ? src.dimension_5() : SrcTraits::value_type::size() )
                                      , ( 6 < SrcType::Rank ? src.dimension_6() : SrcTraits::value_type::size() )
                                      , ( 7 < SrcType::Rank ? src.dimension_7() : SrcTraits::value_type::size() )
                                      );
      }
      else { // is padded
        typedef std::integral_constant<unsigned,sizeof(typename SrcTraits::value_type::value_type)> padded ;

        dst.m_offset = dst_offset_type( padded()
                                      , ( 0 < SrcType::Rank ? src.dimension_0() : SrcTraits::value_type::size() )
                                      , ( 1 < SrcType::Rank ? src.dimension_1() : SrcTraits::value_type::size() )
                                      , ( 2 < SrcType::Rank ? src.dimension_2() : SrcTraits::value_type::size() )
                                      , ( 3 < SrcType::Rank ? src.dimension_3() : SrcTraits::value_type::size() )
                                      , ( 4 < SrcType::Rank ? src.dimension_4() : SrcTraits::value_type::size() )
                                      , ( 5 < SrcType::Rank ? src.dimension_5() : SrcTraits::value_type::size() )
                                      , ( 6 < SrcType::Rank ? src.dimension_6() : SrcTraits::value_type::size() )
                                      , ( 7 < SrcType::Rank ? src.dimension_7() : SrcTraits::value_type::size() )
                                      );
      }

      dst.m_handle = src.m_handle ;
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief  View mapping for non-specialized data type and standard layout */
template< class Traits , bool R0 , bool R1 , bool R2 , bool R3 , bool R4 , bool R5 , bool R6 , bool R7 >
struct SubviewMapping< Traits, R0, R1, R2, R3, R4, R5, R6, R7 ,
  typename std::enable_if<(
    std::is_same< typename Traits::specialize , Kokkos::Array<> >::value
    &&
    (
      std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
    )
  )>::type >
{
private:

  // Subview's rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7) };

  // Whether right-most rank is a range.
  enum { R0_rev = 0 == Traits::rank ? false : (
                  1 == Traits::rank ? R0 : (
                  2 == Traits::rank ? R1 : (
                  3 == Traits::rank ? R2 : (
                  4 == Traits::rank ? R3 : (
                  5 == Traits::rank ? R4 : (
                  6 == Traits::rank ? R5 : (
                  7 == Traits::rank ? R6 : R7 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1 or 2, InputLayout Left, Interval 0
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0 && std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0_rev && std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value )
      ), typename Traits::array_layout , Kokkos::LayoutStride
      >::type array_layout ;

  typedef typename Traits::value_type  value_type ;

  typedef typename std::conditional< rank == 0 , value_type ,
          typename std::conditional< rank == 1 , value_type * ,
          typename std::conditional< rank == 2 , value_type ** ,
          typename std::conditional< rank == 3 , value_type *** ,
          typename std::conditional< rank == 4 , value_type **** ,
          typename std::conditional< rank == 5 , value_type ***** ,
          typename std::conditional< rank == 6 , value_type ****** ,
          typename std::conditional< rank == 7 , value_type ******* ,
                                                 value_type ********
          >::type >::type >::type >::type >::type >::type >::type >::type
     data_type ;

public:

  typedef 
    Kokkos::Experimental::ViewTraits< data_type , array_layout
                                    , typename Traits::device_type
                                    , typename Traits::memory_traits > traits_type ;

  typedef Kokkos::Experimental::View< data_type
                                    , array_layout
                                    , typename Traits::device_type
                                    , typename Traits::memory_traits > type ;

  template< class T0 , class T1 , class T2 , class T3
          , class T4 , class T5 , class T6 , class T7 >
  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< traits_type , void , void > & dst
                    , ViewMapping< Traits , void , void > const & src
                    , T0 const & arg0
                    , T1 const & arg1
                    , T2 const & arg2
                    , T3 const & arg3
                    , T4 const & arg4
                    , T5 const & arg5
                    , T6 const & arg6
                    , T7 const & arg7
                    )
    {
      typedef ViewMapping< traits_type , void , void >  DstType ;

      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T0>  V0 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T1>  V1 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T2>  V2 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T3>  V3 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T4>  V4 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T5>  V5 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T6>  V6 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T7>  V7 ;

      dst.m_offset = dst_offset_type
        ( src.m_offset
        , V0::dimension( src.m_offset.dimension_0() , arg0 )
        , V1::dimension( src.m_offset.dimension_1() , arg1 )
        , V2::dimension( src.m_offset.dimension_2() , arg2 )
        , V3::dimension( src.m_offset.dimension_3() , arg3 )
        , V4::dimension( src.m_offset.dimension_4() , arg4 )
        , V5::dimension( src.m_offset.dimension_5() , arg5 )
        , V6::dimension( src.m_offset.dimension_6() , arg6 )
        , V7::dimension( src.m_offset.dimension_7() , arg7 )
        );

      dst.m_handle = dst_handle_type( src.m_handle +
                                      src.m_offset( V0::begin( arg0 )
                                                  , V1::begin( arg1 )
                                                  , V2::begin( arg2 )
                                                  , V3::begin( arg3 )
                                                  , V4::begin( arg4 )
                                                  , V5::begin( arg5 )
                                                  , V6::begin( arg6 )
                                                  , V7::begin( arg7 )
                                                  ) );
    }
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP */

