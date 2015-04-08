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

#ifndef KOKKOS_VIEWDEFAULT_HPP
#define KOKKOS_VIEWDEFAULT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct ViewAssignment< ViewDefault , ViewDefault , void >
{
  typedef ViewDefault Specialize ;

  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    ||
                    ( ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                      ViewTraits<ST,SL,SD,SM> >::assignable_value
                      &&
                      ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                       typename ViewTraits<ST,SL,SD,SM>::shape_type >::value
                      &&
                      is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout,LayoutStride>::value )
                  )>::type * = 0 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = typename ViewDataManagement< ViewTraits<DT,DL,DD,DM> >::handle_type( src.m_ptr_on_device );

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-1 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                  ), unsigned >::type i0 )
  {
    assert_shape_bounds( src.m_offset_map , 1 , i0 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + i0 ;

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-2 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 )
  {
    assert_shape_bounds( src.m_offset_map , 2 , i0 , i1 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-3 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 )
  {
    assert_shape_bounds( src.m_offset_map, 3, i0, i1, i2 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-4 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 )
  {
    assert_shape_bounds( src.m_offset_map, 4, i0, i1, i2, i3 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2,i3);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-5 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 )
  {
    assert_shape_bounds( src.m_offset_map, 5, i0, i1, i2, i3, i4);

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2,i3,i4);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-6 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 )
  {
    assert_shape_bounds( src.m_offset_map, 6, i0, i1, i2, i3, i4, i5);

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2,i3,i4,i5);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-7 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 )
  {
    assert_shape_bounds( src.m_offset_map, 7, i0, i1, i2, i3, i4, i5, i6 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2,i3,i4,i5,i6);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-0 from Rank-8 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::assignable_value &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 0 ) &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                  ), unsigned >::type i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 ,
                  const unsigned i7 )
  {
    assert_shape_bounds( src.m_offset_map, 8, i0, i1, i2, i3, i4, i5, i6, i7 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,i1,i2,i3,i4,i5,i6,i7);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from range of Rank-1 array, either layout */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ) >::type * = 0 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_offset_map.N0 = 0 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_offset_map , 1 , range.first );
      assert_shape_bounds( src.m_offset_map , 1 , range.second - 1 );

      dst.m_management      = src.m_management ;
      dst.m_offset_map.N0 = range.second - range.first ;
      dst.m_ptr_on_device = src.ptr_on_device() + range.first ;

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array, using ALL as first argument. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ), unsigned >::type i1 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N0 ;
    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(0,i1);

    dst.m_management.increment( dst.m_ptr_on_device );
  }


  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutLeft Rank-2 array, using a row range as first argument. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename IndexType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<IndexType, IndexType>& rowRange,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ), IndexType >::type columnIndex )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    if (rowRange.first < rowRange.second) { // valid row range
      dst.m_management = src.m_management;
      dst.m_offset_map.N0 = rowRange.second - rowRange.first;
      dst.m_ptr_on_device = src.ptr_on_device () +
        src.m_offset_map (rowRange.first, columnIndex);

      dst.m_management.increment( dst.m_ptr_on_device );
    }
    else { // not a valid row range
      dst.m_offset_map.N0 = 0;
      dst.m_ptr_on_device = 0;
    }
  }


  //------------------------------------
  /** \brief  Extract Rank-1 array from LayoutRight Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ), ALL >::type & )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N1 ;
    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), unsigned >::type i1 )
  {
    assert_shape_bounds( src.m_offset_map , 2 , range.first , i1 );
    assert_shape_bounds( src.m_offset_map , 2 , range.second - 1 , i1 );

    dst.m_management.decrement( dst.m_ptr_on_device );

    if ( range.first < range.second ) {
      dst.m_management      = src.m_management ;
      dst.m_offset_map.N0 = range.second - range.first ;
      dst.m_offset_map.N1 = 1 ;
      dst.m_offset_map.S0 = range.second - range.first ;
      dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(range.first,i1);

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutLeft Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), unsigned >::type i1 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N0 ;
    dst.m_offset_map.N1 = 1 ;

    dst.m_offset_map.S0 = src.m_offset_map.N0 ;
    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(0,i1);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-2 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = 1 ;
    dst.m_offset_map.N1 = src.m_offset_map.N1 ;
    dst.m_offset_map.SR = src.m_offset_map.SR ;
    dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(i0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }
  //------------------------------------
  /** \brief  Extract LayoutRight Rank-N array from range of LayoutRight Rank-N array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank > 1 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic > 0 )
                  )>::type * = 0 )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused
    //typedef typename traits_type::shape_type shape_type ; // unused
    //typedef typename View<DT,DL,DD,DM,Specialize>::stride_type stride_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_offset_map.assign( 0, 0, 0, 0, 0, 0, 0, 0 );

    dst.m_ptr_on_device = 0 ;

    if ( ( range.first == range.second ) ||
         ( (src.capacity()==0u) && (range.second<src.m_offset_map.N0) )) {
      dst.m_offset_map.assign( 0 , src.m_offset_map.N1 , src.m_offset_map.N2 , src.m_offset_map.N3 ,
                                   src.m_offset_map.N4 , src.m_offset_map.N5 , src.m_offset_map.N6 , src.m_offset_map.N7 );
      dst.m_offset_map.SR = src.m_offset_map.SR ;
    }
    else if ( (range.first < range.second) ) {
      assert_shape_bounds( src.m_offset_map , 8 , range.first ,      0,0,0,0,0,0,0);
      assert_shape_bounds( src.m_offset_map , 8 , range.second - 1 , 0,0,0,0,0,0,0);

      dst.m_offset_map.assign( range.second - range.first
                             , src.m_offset_map.N1 , src.m_offset_map.N2 , src.m_offset_map.N3
                             , src.m_offset_map.N4 , src.m_offset_map.N5 , src.m_offset_map.N6 , src.m_offset_map.N7 );

      dst.m_offset_map.SR = src.m_offset_map.SR ;

      dst.m_management      = src.m_management ;

      dst.m_ptr_on_device = src.ptr_on_device() + range.first * src.m_offset_map.SR ;

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType0,iType0> & range0 ,
                  const std::pair<iType1,iType1> & range1 ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_offset_map.assign(0,0,0,0, 0,0,0,0);
    dst.m_ptr_on_device = 0 ;

    if ( (range0.first == range0.second) ||
         (range1.first == range1.second) ||
         ( ( src.capacity() == 0u ) &&
           ( long(range0.second) < long(src.m_offset_map.N0) ) &&
           ( long(range1.second) < long(src.m_offset_map.N1) ) ) ) {

      dst.m_offset_map.assign( src.m_offset_map );
      dst.m_offset_map.N0 = range0.second - range0.first ;
      dst.m_offset_map.N1 = range1.second - range1.first ;
    }
    else if ( (range0.first < range0.second && range1.first < range1.second) ) {

      assert_shape_bounds( src.m_offset_map , 2 , range0.first , range1.first );
      assert_shape_bounds( src.m_offset_map , 2 , range0.second - 1 , range1.second - 1 );

      dst.m_offset_map.assign( src.m_offset_map );
      dst.m_offset_map.N0 = range0.second - range0.first ;
      dst.m_offset_map.N1 = range1.second - range1.first ;

      dst.m_management = src.m_management ;

      dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(range0.first,range1.first);

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  ALL ,
                  const std::pair<iType,iType> & range1 ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_offset_map.assign(0,0,0,0, 0,0,0,0);
    dst.m_ptr_on_device = 0 ;

    if ( (range1.first == range1.second) || ( (src.capacity()==0) && (range1.second<src.m_offset_map.N1) )) {
      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N1 = range1.second - range1.first ;
    }
    else if ( (range1.first < range1.second) ) {
      assert_shape_bounds( src.m_offset_map , 2 , 0 , range1.first );
      assert_shape_bounds( src.m_offset_map , 2 , src.m_offset_map.N0 - 1 , range1.second - 1 );

      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N1 = range1.second - range1.first ;
      dst.m_management      = src.m_management ;

      dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(0,range1.first);

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range0 ,
                  ALL ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2
                  ) >::type * = 0 )
  {
    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_offset_map.assign(0,0,0,0, 0,0,0,0);
    dst.m_ptr_on_device = 0 ;

    if ( (range0.first == range0.second) || ( (src.capacity()==0) && (range0.second<src.m_offset_map.N0) )) {
      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N0 = range0.second - range0.first ;
    }
    else if ( (range0.first < range0.second) ) {
      assert_shape_bounds( src.m_offset_map , 2 , range0.first , 0 );
      assert_shape_bounds( src.m_offset_map , 2 , range0.second - 1 , src.m_offset_map.N1 - 1 );

      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N0 = range0.second - range0.first ;
      dst.m_management = src.m_management ;

      dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(range0.first,0);

      dst.m_management.increment( dst.m_ptr_on_device );
    }
  }

  //------------------------------------
  /** \brief  Extract rank-2 from rank-2 array */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const std::pair<iType,iType> & range0 ,
                  ALL ,
                  typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank == 2
                    &&
                    ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1
                  ) >::type * = 0 )
  {
    dst.m_tracking.decrement( dst.ptr_on_device() );
    dst.m_offset_map.assign(0,0,0,0, 0,0,0,0);
    dst.m_ptr_on_device = 0 ;

    if ( (range0.first == range0.second) || ( (src.capacity()==0) && (range0.second<src.m_offset_map.N0) )) {
      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N0 = range0.second - range0.first ;
    }
    else if ( (range0.first < range0.second) ) {
      assert_shape_bounds( src.m_offset_map , 2 , range0.first , 0 );
      assert_shape_bounds( src.m_offset_map , 2 , range0.second - 1 , src.m_offset_map.N1 - 1 );

      dst.m_offset_map.assign(src.m_offset_map);
      dst.m_offset_map.N0 = range0.second - range0.first ;
      dst.m_tracking = src.m_tracking ;

      dst.m_ptr_on_device = src.ptr_on_device() + src.m_offset_map(range0.first,0);

      dst.m_tracking.increment( dst.ptr_on_device() );
    }
  }
  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-3 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N1 ;
    dst.m_offset_map.N1 = src.m_offset_map.N2 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-4 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N2 ;
    dst.m_offset_map.N1 = src.m_offset_map.N3 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,i1,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-5 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N3 ;
    dst.m_offset_map.N1 = src.m_offset_map.N4 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }


  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-6 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N4 ;
    dst.m_offset_map.N1 = src.m_offset_map.N5 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-7 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N5 ;
    dst.m_offset_map.N1 = src.m_offset_map.N6 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,i4,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }


  //------------------------------------
  /** \brief  Extract Rank-2 array from LayoutRight Rank-8 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 2 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 2 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N6 ;
    dst.m_offset_map.N1 = src.m_offset_map.N7 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,i4,i5,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-3 array from LayoutRight Rank-4 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 3 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N1 ;
    dst.m_offset_map.N1 = src.m_offset_map.N2 ;
    dst.m_offset_map.N2 = src.m_offset_map.N3 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 ;
    dst.m_ptr_on_device = &src(i0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-3 array from LayoutRight Rank-5 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 3 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N2 ;
    dst.m_offset_map.N1 = src.m_offset_map.N3 ;
    dst.m_offset_map.N2 = src.m_offset_map.N4 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 ;
    dst.m_ptr_on_device = &src(i0,i1,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-3 array from LayoutRight Rank-6 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 3 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N3 ;
    dst.m_offset_map.N1 = src.m_offset_map.N4 ;
    dst.m_offset_map.N2 = src.m_offset_map.N5 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-3 array from LayoutRight Rank-7 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 3 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N4 ;
    dst.m_offset_map.N1 = src.m_offset_map.N5 ;
    dst.m_offset_map.N2 = src.m_offset_map.N6 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-3 array from LayoutRight Rank-8 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 3 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 3 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N5 ;
    dst.m_offset_map.N1 = src.m_offset_map.N6 ;
    dst.m_offset_map.N2 = src.m_offset_map.N7 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,i4,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-4 array from LayoutRight Rank-5 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 4 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N1 ;
    dst.m_offset_map.N1 = src.m_offset_map.N2 ;
    dst.m_offset_map.N2 = src.m_offset_map.N3 ;
    dst.m_offset_map.N3 = src.m_offset_map.N4 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 ;
    dst.m_ptr_on_device = &src(i0,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-4 array from LayoutRight Rank-6 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 4 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N2 ;
    dst.m_offset_map.N1 = src.m_offset_map.N3 ;
    dst.m_offset_map.N2 = src.m_offset_map.N4 ;
    dst.m_offset_map.N3 = src.m_offset_map.N5 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 ;
    dst.m_ptr_on_device = &src(i0,i1,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-4 array from LayoutRight Rank-7 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 4 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N3 ;
    dst.m_offset_map.N1 = src.m_offset_map.N4 ;
    dst.m_offset_map.N2 = src.m_offset_map.N5 ;
    dst.m_offset_map.N3 = src.m_offset_map.N6 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-4 array from LayoutRight Rank-8 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 4 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 4 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N4 ;
    dst.m_offset_map.N1 = src.m_offset_map.N5 ;
    dst.m_offset_map.N2 = src.m_offset_map.N6 ;
    dst.m_offset_map.N3 = src.m_offset_map.N7 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,i3,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-5 array from LayoutRight Rank-6 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 5 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N1 ;
    dst.m_offset_map.N1 = src.m_offset_map.N2 ;
    dst.m_offset_map.N2 = src.m_offset_map.N3 ;
    dst.m_offset_map.N3 = src.m_offset_map.N4 ;
    dst.m_offset_map.N4 = src.m_offset_map.N5 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 * dst.m_offset_map.N4 ;
    dst.m_ptr_on_device = &src(i0,0,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-5 array from LayoutRight Rank-7 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 5 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N2 ;
    dst.m_offset_map.N1 = src.m_offset_map.N3 ;
    dst.m_offset_map.N2 = src.m_offset_map.N4 ;
    dst.m_offset_map.N3 = src.m_offset_map.N5 ;
    dst.m_offset_map.N4 = src.m_offset_map.N6 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 * dst.m_offset_map.N4 ;
    dst.m_ptr_on_device = &src(i0,i1,0,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Extract Rank-5 array from LayoutRight Rank-8 array. */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const ALL & ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutRight >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank == 5 )
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 5 )
                  ), ALL >::type & )
  {
    //typedef ViewTraits<DT,DL,DD,DM> traits_type ; // unused

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management      = src.m_management ;
    dst.m_offset_map.N0 = src.m_offset_map.N3 ;
    dst.m_offset_map.N1 = src.m_offset_map.N4 ;
    dst.m_offset_map.N2 = src.m_offset_map.N5 ;
    dst.m_offset_map.N3 = src.m_offset_map.N6 ;
    dst.m_offset_map.N4 = src.m_offset_map.N7 ;
    dst.m_offset_map.SR = dst.m_offset_map.N1 * dst.m_offset_map.N2 *
                          dst.m_offset_map.N3 * dst.m_offset_map.N4 ;
    dst.m_ptr_on_device = &src(i0,i1,i2,0,0,0,0,0);

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 1 };

    size_t str[2] = {0,0};

    src.m_offset_map.stride( str );

    const size_t offset = ViewOffsetRange< Type0 >::begin( arg0 ) * str[0] ;

    LayoutStride spec ;

    // Collapse dimension for non-ranges
    if ( ViewOffsetRange< Type0 >::is_range ) {
      spec.dimension[0] = ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 );
      spec.stride[0]    = str[0] ;
    }
    else {
      spec.dimension[0] = 1 ;
      spec.stride[0]    = 1 ;
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 2 };

    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { static_cast<unsigned>(ViewOffsetRange< Type0 >::begin( arg0 ))
      , static_cast<unsigned>(ViewOffsetRange< Type1 >::begin( arg1 ))
      };

    size_t stride[9] ;

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    spec.dimension[0] = ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 );
    spec.dimension[1] = ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 );
    spec.stride[0]    = stride[0] ;
    spec.stride[1]    = stride[1] ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = spec.dimension[i] ;
      spec.stride[j]    = spec.stride[i] ;
      offset += begin[i] * spec.stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 3 };

    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      };

    // FIXME (mfh 26 Oct 2014) Should use size_type typedef here
    // instead of unsigned.  If we did that, the static_casts would be
    // unnecessary.
    const unsigned begin[ src_rank ] = {
      static_cast<unsigned> (ViewOffsetRange< Type0 >::begin (arg0))
      , static_cast<unsigned> (ViewOffsetRange< Type1 >::begin (arg1))
      , static_cast<unsigned> (ViewOffsetRange< Type2 >::begin (arg2))
    };

    // FIXME (mfh 26 Oct 2014) Should use size_type typedef here
    // instead of unsigned.  If we did that, the static_casts would be
    // unnecessary.
    unsigned dim[ src_rank ] = {
      static_cast<unsigned> (ViewOffsetRange< Type0 >::dimension (src.m_offset_map.N0, arg0))
      , static_cast<unsigned> (ViewOffsetRange< Type1 >::dimension (src.m_offset_map.N1, arg1))
      , static_cast<unsigned> (ViewOffsetRange< Type2 >::dimension (src.m_offset_map.N2, arg2))
    };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          , class Type3
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const Type3 & arg3 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type3 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 4 };
    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      , ViewOffsetRange< Type3 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { static_cast<unsigned>(ViewOffsetRange< Type0 >::begin( arg0 ))
      , static_cast<unsigned>(ViewOffsetRange< Type1 >::begin( arg1 ))
      , static_cast<unsigned>(ViewOffsetRange< Type2 >::begin( arg2 ))
      , static_cast<unsigned>(ViewOffsetRange< Type3 >::begin( arg3 ))
      };

    unsigned dim[ src_rank ] =
      { static_cast<unsigned>(ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 ))
      , static_cast<unsigned>(ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 ))
      , static_cast<unsigned>(ViewOffsetRange< Type2 >::dimension( src.m_offset_map.N2 , arg2 ))
      , static_cast<unsigned>(ViewOffsetRange< Type3 >::dimension( src.m_offset_map.N3 , arg3 ))
      };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          , class Type3
          , class Type4
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const Type3 & arg3 ,
                  const Type4 & arg4 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type3 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type4 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 5 };
    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      , ViewOffsetRange< Type3 >::is_range
      , ViewOffsetRange< Type4 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { ViewOffsetRange< Type0 >::begin( arg0 )
      , ViewOffsetRange< Type1 >::begin( arg1 )
      , ViewOffsetRange< Type2 >::begin( arg2 )
      , ViewOffsetRange< Type3 >::begin( arg3 )
      , ViewOffsetRange< Type4 >::begin( arg4 )
      };

    unsigned dim[ src_rank ] =
      { ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 )
      , ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 )
      , ViewOffsetRange< Type2 >::dimension( src.m_offset_map.N2 , arg2 )
      , ViewOffsetRange< Type3 >::dimension( src.m_offset_map.N3 , arg3 )
      , ViewOffsetRange< Type4 >::dimension( src.m_offset_map.N4 , arg4 )
      };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          , class Type3
          , class Type4
          , class Type5
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const Type3 & arg3 ,
                  const Type4 & arg4 ,
                  const Type5 & arg5 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type3 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type4 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type5 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 6 };
    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      , ViewOffsetRange< Type3 >::is_range
      , ViewOffsetRange< Type4 >::is_range
      , ViewOffsetRange< Type5 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { ViewOffsetRange< Type0 >::begin( arg0 )
      , ViewOffsetRange< Type1 >::begin( arg1 )
      , ViewOffsetRange< Type2 >::begin( arg2 )
      , ViewOffsetRange< Type3 >::begin( arg3 )
      , ViewOffsetRange< Type4 >::begin( arg4 )
      , ViewOffsetRange< Type5 >::begin( arg5 )
      };

    unsigned dim[ src_rank ] =
      { ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 )
      , ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 )
      , ViewOffsetRange< Type2 >::dimension( src.m_offset_map.N2 , arg2 )
      , ViewOffsetRange< Type3 >::dimension( src.m_offset_map.N3 , arg3 )
      , ViewOffsetRange< Type4 >::dimension( src.m_offset_map.N4 , arg4 )
      , ViewOffsetRange< Type5 >::dimension( src.m_offset_map.N5 , arg5 )
      };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          , class Type3
          , class Type4
          , class Type5
          , class Type6
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const Type3 & arg3 ,
                  const Type4 & arg4 ,
                  const Type5 & arg5 ,
                  const Type6 & arg6 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type3 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type4 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type5 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type6 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 7 };
    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      , ViewOffsetRange< Type3 >::is_range
      , ViewOffsetRange< Type4 >::is_range
      , ViewOffsetRange< Type5 >::is_range
      , ViewOffsetRange< Type6 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { ViewOffsetRange< Type0 >::begin( arg0 )
      , ViewOffsetRange< Type1 >::begin( arg1 )
      , ViewOffsetRange< Type2 >::begin( arg2 )
      , ViewOffsetRange< Type3 >::begin( arg3 )
      , ViewOffsetRange< Type4 >::begin( arg4 )
      , ViewOffsetRange< Type5 >::begin( arg5 )
      , ViewOffsetRange< Type6 >::begin( arg6 )
      };

    unsigned dim[ src_rank ] =
      { ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 )
      , ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 )
      , ViewOffsetRange< Type2 >::dimension( src.m_offset_map.N2 , arg2 )
      , ViewOffsetRange< Type3 >::dimension( src.m_offset_map.N3 , arg3 )
      , ViewOffsetRange< Type4 >::dimension( src.m_offset_map.N4 , arg4 )
      , ViewOffsetRange< Type5 >::dimension( src.m_offset_map.N5 , arg5 )
      , ViewOffsetRange< Type6 >::dimension( src.m_offset_map.N6 , arg6 )
      };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );
    dst.m_management = src.m_management ;
    dst.m_offset_map.assign( spec );
    dst.m_ptr_on_device = src.ptr_on_device() + offset ;
    dst.m_management.increment( dst.m_ptr_on_device );
  }

  template< class DT , class DL , class DD , class DM
          , class ST , class SL , class SD , class SM
          , class Type0
          , class Type1
          , class Type2
          , class Type3
          , class Type4
          , class Type5
          , class Type6
          , class Type7
          >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const Type0 & arg0 ,
                  const Type1 & arg1 ,
                  const Type2 & arg2 ,
                  const Type3 & arg3 ,
                  const Type4 & arg4 ,
                  const Type5 & arg5 ,
                  const Type6 & arg6 ,
                  const Type7 & arg7 ,
                  const typename enable_if< (
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::assignable_value
                    &&
                    is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout , LayoutStride >::value
                    &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) ==
                      ( ViewOffsetRange< Type0 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type1 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type2 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type3 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type4 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type5 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type6 >::is_range ? 1u : 0 ) +
                      ( ViewOffsetRange< Type7 >::is_range ? 1u : 0 ) )
                  )>::type * = 0 )
  {
    enum { src_rank = 8 };

    const bool is_range[ src_rank ] =
      { ViewOffsetRange< Type0 >::is_range
      , ViewOffsetRange< Type1 >::is_range
      , ViewOffsetRange< Type2 >::is_range
      , ViewOffsetRange< Type3 >::is_range
      , ViewOffsetRange< Type4 >::is_range
      , ViewOffsetRange< Type5 >::is_range
      , ViewOffsetRange< Type6 >::is_range
      , ViewOffsetRange< Type7 >::is_range
      };

    const unsigned begin[ src_rank ] =
      { ViewOffsetRange< Type0 >::begin( arg0 )
      , ViewOffsetRange< Type1 >::begin( arg1 )
      , ViewOffsetRange< Type2 >::begin( arg2 )
      , ViewOffsetRange< Type3 >::begin( arg3 )
      , ViewOffsetRange< Type4 >::begin( arg4 )
      , ViewOffsetRange< Type5 >::begin( arg5 )
      , ViewOffsetRange< Type6 >::begin( arg6 )
      , ViewOffsetRange< Type7 >::begin( arg7 )
      };

    unsigned dim[ src_rank ] =
      { ViewOffsetRange< Type0 >::dimension( src.m_offset_map.N0 , arg0 )
      , ViewOffsetRange< Type1 >::dimension( src.m_offset_map.N1 , arg1 )
      , ViewOffsetRange< Type2 >::dimension( src.m_offset_map.N2 , arg2 )
      , ViewOffsetRange< Type3 >::dimension( src.m_offset_map.N3 , arg3 )
      , ViewOffsetRange< Type4 >::dimension( src.m_offset_map.N4 , arg4 )
      , ViewOffsetRange< Type5 >::dimension( src.m_offset_map.N5 , arg5 )
      , ViewOffsetRange< Type6 >::dimension( src.m_offset_map.N6 , arg6 )
      , ViewOffsetRange< Type7 >::dimension( src.m_offset_map.N7 , arg7 )
      };

    size_t stride[9] = {0,0,0,0,0,0,0,0,0};

    src.m_offset_map.stride( stride );

    LayoutStride spec ;

    size_t offset = 0 ;

    // Collapse dimension for non-ranges
    for ( int i = 0 , j = 0 ; i < int(src_rank) ; ++i ) {
      spec.dimension[j] = dim[i] ;
      spec.stride[j]    = stride[i] ;
      offset += begin[i] * stride[i] ;
      if ( is_range[i] ) { ++j ; }
    }

    dst.m_management.decrement( dst.m_ptr_on_device );

    dst.m_management = src.m_management ;

    dst.m_offset_map.assign( spec );

    dst.m_ptr_on_device = src.ptr_on_device() + offset ;

    dst.m_management.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Deep copy data from compatible value type, layout, rank, and specialization.
   *          Check the dimensions and allocation lengths at runtime.
   */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline static
  void deep_copy( const View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename Impl::enable_if<(
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::value_type ,
                                   typename ViewTraits<ST,SL,SD,SM>::non_const_value_type >::value
                    &&
                    Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                                   typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                    &&
                    ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                  )>::type * = 0 )
  {
    typedef typename ViewTraits<DT,DL,DD,DM>::memory_space dst_memory_space ;
    typedef typename ViewTraits<ST,SL,SD,SM>::memory_space src_memory_space ;

    if ( dst.ptr_on_device() != src.ptr_on_device() ) {

      Impl::assert_shapes_are_equal( dst.m_offset_map , src.m_offset_map );

      const size_t nbytes = dst.m_offset_map.scalar_size * dst.m_offset_map.capacity();

      DeepCopy< dst_memory_space , src_memory_space >( dst.ptr_on_device() , src.ptr_on_device() , nbytes );
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type , class SubArg7_type
        >
struct ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                  , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                  , SubArg4_type , SubArg5_type , SubArg6_type , SubArg7_type >
{
private:

  typedef View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >  SrcViewType ;

  enum { V0 = Impl::is_same< SubArg0_type , void >::value ? 1 : 0 };
  enum { V1 = Impl::is_same< SubArg1_type , void >::value ? 1 : 0 };
  enum { V2 = Impl::is_same< SubArg2_type , void >::value ? 1 : 0 };
  enum { V3 = Impl::is_same< SubArg3_type , void >::value ? 1 : 0 };
  enum { V4 = Impl::is_same< SubArg4_type , void >::value ? 1 : 0 };
  enum { V5 = Impl::is_same< SubArg5_type , void >::value ? 1 : 0 };
  enum { V6 = Impl::is_same< SubArg6_type , void >::value ? 1 : 0 };
  enum { V7 = Impl::is_same< SubArg7_type , void >::value ? 1 : 0 };

  // The source view rank must be equal to the input argument rank
  // Once a void argument is encountered all subsequent arguments must be void.
  enum { InputRank =
    Impl::StaticAssert<( SrcViewType::rank ==
                         ( V0 ? 0 : (
                           V1 ? 1 : (
                           V2 ? 2 : (
                           V3 ? 3 : (
                           V4 ? 4 : (
                           V5 ? 5 : (
                           V6 ? 6 : (
                           V7 ? 7 : 8 ))))))) ))
                       &&
                       ( SrcViewType::rank ==
                         ( 8 - ( V0 + V1 + V2 + V3 + V4 + V5 + V6 + V7 ) ) )
    >::value ? SrcViewType::rank : 0 };

  enum { R0 = Impl::ViewOffsetRange< SubArg0_type >::is_range ? 1 : 0 };
  enum { R1 = Impl::ViewOffsetRange< SubArg1_type >::is_range ? 1 : 0 };
  enum { R2 = Impl::ViewOffsetRange< SubArg2_type >::is_range ? 1 : 0 };
  enum { R3 = Impl::ViewOffsetRange< SubArg3_type >::is_range ? 1 : 0 };
  enum { R4 = Impl::ViewOffsetRange< SubArg4_type >::is_range ? 1 : 0 };
  enum { R5 = Impl::ViewOffsetRange< SubArg5_type >::is_range ? 1 : 0 };
  enum { R6 = Impl::ViewOffsetRange< SubArg6_type >::is_range ? 1 : 0 };
  enum { R7 = Impl::ViewOffsetRange< SubArg7_type >::is_range ? 1 : 0 };

  enum { OutputRank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
                    + unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7) };

  // Reverse
  enum { R0_rev = 0 == InputRank ? 0u : (
                  1 == InputRank ? unsigned(R0) : (
                  2 == InputRank ? unsigned(R1) : (
                  3 == InputRank ? unsigned(R2) : (
                  4 == InputRank ? unsigned(R3) : (
                  5 == InputRank ? unsigned(R4) : (
                  6 == InputRank ? unsigned(R5) : (
                  7 == InputRank ? unsigned(R6) : unsigned(R7) ))))))) };

  typedef typename SrcViewType::array_layout  SrcViewLayout ;

  // Choose array layout, attempting to preserve original layout if at all possible.
  typedef typename Impl::if_c<
     ( // Same Layout IF
       // OutputRank 0
       ( OutputRank == 0 )
       ||
       // OutputRank 1 or 2, InputLayout Left, Interval 0
       // because single stride one or second index has a stride.
       ( OutputRank <= 2 && R0 && Impl::is_same<SrcViewLayout,LayoutLeft>::value )
       ||
       // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
       // because single stride one or second index has a stride.
       ( OutputRank <= 2 && R0_rev && Impl::is_same<SrcViewLayout,LayoutRight>::value )
     ), SrcViewLayout , Kokkos::LayoutStride >::type OutputViewLayout ;

  // Choose data type as a purely dynamic rank array to accomodate a runtime range.
  typedef typename Impl::if_c< OutputRank == 0 , typename SrcViewType::value_type ,
          typename Impl::if_c< OutputRank == 1 , typename SrcViewType::value_type *,
          typename Impl::if_c< OutputRank == 2 , typename SrcViewType::value_type **,
          typename Impl::if_c< OutputRank == 3 , typename SrcViewType::value_type ***,
          typename Impl::if_c< OutputRank == 4 , typename SrcViewType::value_type ****,
          typename Impl::if_c< OutputRank == 5 , typename SrcViewType::value_type *****,
          typename Impl::if_c< OutputRank == 6 , typename SrcViewType::value_type ******,
          typename Impl::if_c< OutputRank == 7 , typename SrcViewType::value_type *******,
                                                 typename SrcViewType::value_type ********
  >::type >::type >::type >::type >::type >::type >::type >::type  OutputData ;

  // Choose space.
  // If the source view's template arg1 or arg2 is a space then use it,
  // otherwise use the source view's execution space.

  typedef typename Impl::if_c< Impl::is_space< SrcArg1Type >::value , SrcArg1Type ,
          typename Impl::if_c< Impl::is_space< SrcArg2Type >::value , SrcArg2Type , typename SrcViewType::execution_space 
  >::type >::type OutputSpace ;

public:

  // If keeping the layout then match non-data type arguments
  // else keep execution space and memory traits.
  typedef typename
    Impl::if_c< Impl::is_same< SrcViewLayout , OutputViewLayout >::value
              , Kokkos::View< OutputData , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
              , Kokkos::View< OutputData , OutputViewLayout , OutputSpace
                            , typename SrcViewType::memory_traits
                            , Impl::ViewDefault >
              >::type  type ;
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

// Construct subview of a Rank 8 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type , class SubArg7_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    , const SubArg6_type & arg6
    , const SubArg7_type & arg7
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , SubArg6_type , SubArg7_type >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;
    typedef Impl::ViewOffsetRange< SubArg6_type > R6 ;
    typedef Impl::ViewOffsetRange< SubArg7_type > R7 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , R6::dimension( src.m_offset_map.N6 , arg6 )
                                 , R7::dimension( src.m_offset_map.N7 , arg7 )
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , R6::begin( arg6 )
                                        , R7::begin( arg7 ) );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 7 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type , class SubArg6_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    , const SubArg6_type & arg6
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , SubArg6_type , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;
    typedef Impl::ViewOffsetRange< SubArg6_type > R6 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , R6::dimension( src.m_offset_map.N6 , arg6 )
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , R6::begin( arg6 )
                                        , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 6 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type , class SubArg5_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    , const SubArg5_type & arg5
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , SubArg5_type , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;
    typedef Impl::ViewOffsetRange< SubArg5_type > R5 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , R5::dimension( src.m_offset_map.N5 , arg5 )
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , R5::begin( arg5 )
                                        , 0
                                        , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 5 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        , class SubArg4_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    , const SubArg4_type & arg4
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , SubArg4_type , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;
    typedef Impl::ViewOffsetRange< SubArg4_type > R4 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , R4::dimension( src.m_offset_map.N4 , arg4 )
                                 , 0
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , R4::begin( arg4 )
                                        , 0
                                        , 0
                                        , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 4 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type , class SubArg3_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    , const SubArg3_type & arg3
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , SubArg3_type
                           , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;
    typedef Impl::ViewOffsetRange< SubArg3_type > R3 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , R3::dimension( src.m_offset_map.N3 , arg3 )
                                 , 0
                                 , 0
                                 , 0
                                 , 0
                                 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , R3::begin( arg3 )
                                        , 0
                                        , 0
                                        , 0
                                        , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 3 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type , class SubArg2_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    , const SubArg2_type & arg2
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , SubArg2_type , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;
    typedef Impl::ViewOffsetRange< SubArg2_type > R2 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , R2::dimension( src.m_offset_map.N2 , arg2 )
                                 , 0 , 0 , 0 , 0 , 0);

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , R2::begin( arg2 )
                                        , 0 , 0 , 0 , 0 , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 2 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type , class SubArg1_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    , const SubArg1_type & arg1
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , SubArg1_type , void , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;
    typedef Impl::ViewOffsetRange< SubArg1_type > R1 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , R1::dimension( src.m_offset_map.N1 , arg1 )
                                 , 0 , 0 , 0 , 0 , 0 , 0 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , R1::begin( arg1 )
                                        , 0 , 0 , 0 , 0 , 0 , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

// Construct subview of a Rank 1 view
template< class DstDataType , class DstArg1Type , class DstArg2Type , class DstArg3Type >
template< class SrcDataType , class SrcArg1Type , class SrcArg2Type , class SrcArg3Type
        , class SubArg0_type
        >
KOKKOS_INLINE_FUNCTION
View< DstDataType , DstArg1Type , DstArg2Type , DstArg3Type , Impl::ViewDefault >::
View( const View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault > & src
    , const SubArg0_type & arg0
    )
  : m_ptr_on_device( (typename traits::value_type*) NULL)
  , m_offset_map()
  , m_management()
{
  // This constructor can only be used to construct a subview
  // from the source view.  This type must match the subview type
  // deduced from the source view and subview arguments.

  typedef Impl::ViewSubview< View< SrcDataType , SrcArg1Type , SrcArg2Type , SrcArg3Type , Impl::ViewDefault >
                           , SubArg0_type , void , void , void , void , void , void , void >
    ViewSubviewDeduction ;

  enum { is_a_valid_subview_constructor =
    Impl::StaticAssert<
      Impl::is_same< View , typename ViewSubviewDeduction::type >::value
    >::value
  };

  if ( is_a_valid_subview_constructor ) {

    typedef Impl::ViewOffsetRange< SubArg0_type > R0 ;

    // 'assign_subview' returns whether the subview offset_map
    // introduces noncontiguity in the view.
    const bool introduce_noncontiguity =
      m_offset_map.assign_subview( src.m_offset_map
                                 , R0::dimension( src.m_offset_map.N0 , arg0 )
                                 , 0 , 0 , 0 , 0 , 0 , 0 , 0 );

    if ( m_offset_map.capacity() ) {

      m_management = src.m_management ;

      if ( introduce_noncontiguity ) m_management.set_noncontiguous();

      m_ptr_on_device = src.m_ptr_on_device +
                        src.m_offset_map( R0::begin( arg0 )
                                        , 0 , 0 , 0 , 0 , 0 , 0 , 0 );
      m_management.increment( m_ptr_on_device );
    }
  }
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWDEFAULT_HPP */

