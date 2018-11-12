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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP

#include <Kokkos_Array.hpp>

namespace Kokkos {
namespace Impl {

template< class DataType , class ArrayLayout , class V , size_t N , class P >
struct ViewDataAnalysis< DataType , ArrayLayout , Kokkos::Array<V,N,P> >
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

  typedef typename dimension::template append<N>::type array_scalar_dimension ;

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

  typedef typename ViewDataType<           scalar_type , array_scalar_dimension >::type  scalar_array_type ;
  typedef typename ViewDataType<     const_scalar_type , array_scalar_dimension >::type  const_scalar_array_type ;
  typedef typename ViewDataType< non_const_scalar_type , array_scalar_dimension >::type  non_const_scalar_array_type ;
};

}} // namespace Kokkos::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  View mapping for non-specialized data type and standard layout */
template< class Traits >
class ViewMapping< Traits ,
  typename std::enable_if<(
    std::is_same< typename Traits::specialize , Kokkos::Array<> >::value &&
    ( std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value )
  )>::type >
{
private:

  template< class , class ... > friend class ViewMapping ;
  template< class , class ... > friend class Kokkos::View ;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  typedef typename Traits::value_type::pointer handle_type ;

  handle_type  m_impl_handle ;
  offset_type  m_impl_offset ;
  size_t       m_stride ;

  typedef typename Traits::value_type::value_type scalar_type ;

  typedef Kokkos::Array< scalar_type ,KOKKOS_INVALID_INDEX , Kokkos::Array<>::contiguous >  contiguous_reference ;
  typedef Kokkos::Array< scalar_type ,KOKKOS_INVALID_INDEX , Kokkos::Array<>::strided >     strided_reference ;

  enum { is_contiguous_reference =
    ( Traits::rank == 0 ) || ( std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ) };

  enum { Array_N = Traits::value_type::size() };
  enum { Array_S = is_contiguous_reference ? Array_N : 1 };

  KOKKOS_INLINE_FUNCTION
  ViewMapping( const handle_type & arg_handle , const offset_type & arg_offset )
    : m_impl_handle( arg_handle )
    , m_impl_offset( arg_offset )
    , m_stride( is_contiguous_reference ? 0 : arg_offset.span() )
    {}

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr size_t extent( const iType & r ) const
    { return m_impl_offset.m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr
  typename Traits::array_layout layout() const
    { return m_impl_offset.layout(); }


  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const { return m_impl_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return m_impl_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return m_impl_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return m_impl_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return m_impl_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return m_impl_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return m_impl_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return m_impl_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  using is_regular = typename offset_type::is_regular ;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return m_impl_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return m_impl_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return m_impl_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return m_impl_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return m_impl_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return m_impl_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return m_impl_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return m_impl_offset.stride_7(); }

  //----------------------------------------
  // Range span

  /** \brief  Span of the mapped range */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const
    { return m_impl_offset.span() * Array_N ; }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_impl_offset.span_is_contiguous(); }

  typedef typename std::conditional< is_contiguous_reference , contiguous_reference , strided_reference >::type  reference_type ;

  typedef handle_type pointer_type ;

  /** \brief  If data references are lvalue_reference than can query pointer to memory */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const
    { return m_impl_handle ; }

  //----------------------------------------
  // The View class performs all rank and bounds checking before
  // calling these element reference methods.

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const { return reference_type( m_impl_handle + 0 , Array_N , 0 ); }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type
  reference( const I0 & i0 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2,i3) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2,i3,i4) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2,i3,i4,i5) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2,i3,i4,i5,i6) * Array_S , Array_N , m_stride ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 ) const
    { return reference_type( m_impl_handle + m_impl_offset(i0,i1,i2,i3,i4,i5,i6,i7) * Array_S , Array_N , m_stride ); }

  //----------------------------------------

private:

  enum { MemorySpanMask = 8 - 1 /* Force alignment on 8 byte boundary */ };
  enum { MemorySpanSize = sizeof(scalar_type) };

public:

  /** \brief  Span, in bytes, of the referenced memory */
  KOKKOS_INLINE_FUNCTION constexpr size_t memory_span() const
    {
      return ( m_impl_offset.span() * Array_N * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION ~ViewMapping() {}
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_impl_handle(), m_impl_offset(), m_stride(0) {}
  KOKKOS_INLINE_FUNCTION ViewMapping( const ViewMapping & rhs )
    : m_impl_handle( rhs.m_impl_handle ), m_impl_offset( rhs.m_impl_offset ), m_stride( rhs.m_stride ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( const ViewMapping & rhs )
    { m_impl_handle = rhs.m_impl_handle ; m_impl_offset = rhs.m_impl_offset ; m_stride = rhs.m_stride ; ; return *this ; }

  KOKKOS_INLINE_FUNCTION ViewMapping( ViewMapping && rhs )
    : m_impl_handle( rhs.m_impl_handle ), m_impl_offset( rhs.m_impl_offset ), m_stride( rhs.m_stride ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( ViewMapping && rhs )
    { m_impl_handle = rhs.m_impl_handle ; m_impl_offset = rhs.m_impl_offset ; m_stride = rhs.m_stride ; return *this ; }

  //----------------------------------------

  template< class ... Args >
  KOKKOS_INLINE_FUNCTION
  ViewMapping( pointer_type ptr , Args ... args )
    : m_impl_handle( ptr )
    , m_impl_offset( std::integral_constant< unsigned , 0 >() , args... )
    , m_stride( m_impl_offset.span() )
    {}

  //----------------------------------------

  template< class ... P >
  Kokkos::Impl::SharedAllocationRecord<> *
  allocate_shared( Kokkos::Impl::ViewCtorProp< P... > const & arg_prop
                 , typename Traits::array_layout const & arg_layout
                 )
  {
    typedef Kokkos::Impl::ViewCtorProp< P... > alloc_prop ;

    typedef typename alloc_prop::execution_space  execution_space ;
    typedef typename Traits::memory_space         memory_space ;
    typedef ViewValueFunctor< execution_space , scalar_type > functor_type ;
    typedef Kokkos::Impl::SharedAllocationRecord< memory_space , functor_type > record_type ;

    // Query the mapping for byte-size of allocation.
    typedef std::integral_constant< unsigned ,
      alloc_prop::allow_padding ? sizeof(scalar_type) : 0 > padding ;

    m_impl_offset = offset_type( padding(), arg_layout );

    const size_t alloc_size =
      ( m_impl_offset.span() * Array_N * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);

    // Allocate memory from the memory space and create tracking record.
    record_type * const record =
      record_type::allocate( ((Kokkos::Impl::ViewCtorProp<void,memory_space> const &) arg_prop ).value
                           , ((Kokkos::Impl::ViewCtorProp<void,std::string>  const &) arg_prop ).value
                           , alloc_size );

    if ( alloc_size ) {
      m_impl_handle =
        handle_type( reinterpret_cast< pointer_type >( record->data() ) );

      if ( alloc_prop::initialize ) {
        // The functor constructs and destroys
        record->m_destroy = functor_type( ((Kokkos::Impl::ViewCtorProp<void,execution_space> const & )arg_prop).value
                                        , (pointer_type) m_impl_handle
                                        , m_impl_offset.span() * Array_N
                                        );

        record->m_destroy.construct_shared_allocation();
      }
    }

    return record ;
  }
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

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcType ;

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

      dst.m_impl_offset = dst_offset_type( src.m_impl_offset );
      dst.m_impl_handle = src.m_impl_handle ;
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

  enum { is_assignable = std::is_same< typename DstTraits::data_type ,    typename SrcTraits::scalar_array_type >::value &&
                         std::is_same< typename DstTraits::array_layout , typename SrcTraits::array_layout >::value };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void >  DstType ;
  typedef ViewMapping< SrcTraits , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION
  static void assign( DstType & dst , const SrcType & src , const TrackType & src_track )
    {
      static_assert( is_assignable , "Can only convert to array_type" );

      typedef typename DstType::offset_type  dst_offset_type ;

      // Array dimension becomes the last dimension.
      // Arguments beyond the destination rank are ignored.
      if ( src.span_is_contiguous() ) { // not padded
        dst.m_impl_offset = dst_offset_type( std::integral_constant<unsigned,0>() ,
          typename DstTraits::array_layout
            ( ( 0 < SrcType::Rank ? src.dimension_0() : SrcTraits::value_type::size() )
            , ( 1 < SrcType::Rank ? src.dimension_1() : SrcTraits::value_type::size() )
            , ( 2 < SrcType::Rank ? src.dimension_2() : SrcTraits::value_type::size() )
            , ( 3 < SrcType::Rank ? src.dimension_3() : SrcTraits::value_type::size() )
            , ( 4 < SrcType::Rank ? src.dimension_4() : SrcTraits::value_type::size() )
            , ( 5 < SrcType::Rank ? src.dimension_5() : SrcTraits::value_type::size() )
            , ( 6 < SrcType::Rank ? src.dimension_6() : SrcTraits::value_type::size() )
            , ( 7 < SrcType::Rank ? src.dimension_7() : SrcTraits::value_type::size() )
            ) );
      }
      else { // is padded
        typedef std::integral_constant<unsigned,sizeof(typename SrcTraits::value_type::value_type)> padded ;

        dst.m_impl_offset = dst_offset_type( padded() ,
          typename DstTraits::array_layout
            ( ( 0 < SrcType::Rank ? src.dimension_0() : SrcTraits::value_type::size() )
            , ( 1 < SrcType::Rank ? src.dimension_1() : SrcTraits::value_type::size() )
            , ( 2 < SrcType::Rank ? src.dimension_2() : SrcTraits::value_type::size() )
            , ( 3 < SrcType::Rank ? src.dimension_3() : SrcTraits::value_type::size() )
            , ( 4 < SrcType::Rank ? src.dimension_4() : SrcTraits::value_type::size() )
            , ( 5 < SrcType::Rank ? src.dimension_5() : SrcTraits::value_type::size() )
            , ( 6 < SrcType::Rank ? src.dimension_6() : SrcTraits::value_type::size() )
            , ( 7 < SrcType::Rank ? src.dimension_7() : SrcTraits::value_type::size() )
            ) );
      }

      dst.m_impl_handle = src.m_impl_handle ;
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class SrcTraits , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      std::is_same< typename SrcTraits::specialize , Kokkos::Array<> >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutStride >::value
      )
    )>::type
  , SrcTraits
  , Args ... >
{
private:

  static_assert( SrcTraits::rank == sizeof...(Args) , "" );

  enum : bool
    { R0 = is_integral_extent<0,Args...>::value
    , R1 = is_integral_extent<1,Args...>::value
    , R2 = is_integral_extent<2,Args...>::value
    , R3 = is_integral_extent<3,Args...>::value
    , R4 = is_integral_extent<4,Args...>::value
    , R5 = is_integral_extent<5,Args...>::value
    , R6 = is_integral_extent<6,Args...>::value
    , R7 = is_integral_extent<7,Args...>::value
    };

  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7) };

  // Whether right-most rank is a range.
  enum { R0_rev = 0 == SrcTraits::rank ? false : (
                  1 == SrcTraits::rank ? R0 : (
                  2 == SrcTraits::rank ? R1 : (
                  3 == SrcTraits::rank ? R2 : (
                  4 == SrcTraits::rank ? R3 : (
                  5 == SrcTraits::rank ? R4 : (
                  6 == SrcTraits::rank ? R5 : (
                  7 == SrcTraits::rank ? R6 : R7 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1 or 2, InputLayout Left, Interval 0
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0 && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0_rev && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value )
      ), typename SrcTraits::array_layout , Kokkos::LayoutStride
      >::type array_layout ;

  typedef typename SrcTraits::value_type  value_type ;

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

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;

  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< traits_type , void > & dst
                    , ViewMapping< SrcTraits , void > const & src
                    , Args ... args )
    {
      typedef ViewMapping< traits_type , void >  DstType ;

      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      const SubviewExtents< SrcTraits::rank , rank >
        extents( src.m_impl_offset.m_dim , args... );

      dst.m_impl_offset = dst_offset_type( src.m_impl_offset , extents );
      dst.m_impl_handle = dst_handle_type( src.m_impl_handle +
                                      src.m_impl_offset( extents.domain_offset(0)
                                                  , extents.domain_offset(1)
                                                  , extents.domain_offset(2)
                                                  , extents.domain_offset(3)
                                                  , extents.domain_offset(4)
                                                  , extents.domain_offset(5)
                                                  , extents.domain_offset(6)
                                                  , extents.domain_offset(7)
                                                  ) );
    }
};

}} // namespace Kokkos::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_ARRAY_MAPPING_HPP */

