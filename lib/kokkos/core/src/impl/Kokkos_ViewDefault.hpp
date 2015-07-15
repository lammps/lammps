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
  /** \brief  Compatible value and shape and LayoutLeft/Right to LayoutStride*/

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
                      is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout,LayoutStride>::value
                      && (is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout,LayoutLeft>::value ||
                          is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout,LayoutRight>::value))
                  )>::type * = 0 )
  {
    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = ViewDataManagement< ViewTraits<DT,DL,DD,DM> >::create_handle( src.m_ptr_on_device, src.m_tracker );

    dst.m_tracker = src.m_tracker ;

  }


  /** \brief  Assign 1D Strided View to LayoutLeft or LayoutRight if stride[0]==1 */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,LayoutStride,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    (
                      ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,LayoutStride,SD,SM> >::value
                      ||
                      ( ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                      ViewTraits<ST,LayoutStride,SD,SM> >::assignable_value
                        &&
                        ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                       typename ViewTraits<ST,LayoutStride,SD,SM>::shape_type >::value
                      )
                     )
                     &&
                      (View<DT,DL,DD,DM,Specialize>::rank==1)
                     && (is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout,LayoutLeft>::value ||
                          is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout,LayoutRight>::value)
                  )>::type * = 0 )
  {
    size_t strides[8];
    src.stride(strides);
    if(strides[0]!=1) {
      abort("Trying to assign strided 1D View to LayoutRight or LayoutLeft which is not stride-1");
    }
    dst.m_offset_map.assign( src.dimension_0(), 0, 0, 0, 0, 0, 0, 0, 0 );

    dst.m_management = src.m_management ;

    dst.m_ptr_on_device = ViewDataManagement< ViewTraits<DT,DL,DD,DM> >::create_handle( src.m_ptr_on_device, src.m_tracker );

    dst.m_tracker = src.m_tracker ;

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

template< class ExecSpace , class DT , class DL, class DD, class DM, class DS >
struct ViewDefaultConstruct< ExecSpace , Kokkos::View<DT,DL,DD,DM,DS> , true >
{
  Kokkos::View<DT,DL,DD,DM,DS> * const m_ptr ;

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const typename ExecSpace::size_type& i ) const
    { new(m_ptr+i) Kokkos::View<DT,DL,DD,DM,DS>(); }

  ViewDefaultConstruct( Kokkos::View<DT,DL,DD,DM,DS> * pointer , size_t capacity )
    : m_ptr( pointer )
    {
      Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
      parallel_for( range , *this );
      ExecSpace::fence();
    }
};

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
          typename Impl::if_c< Impl::is_space< SrcArg2Type >::value , SrcArg2Type , typename SrcViewType::device_type
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
  , m_tracker()
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
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
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
  , m_tracker()
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
                                        );
      m_tracker = src.m_tracker ;
    }
  }
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWDEFAULT_HPP */

