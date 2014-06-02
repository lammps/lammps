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

#ifndef KOKKOS_VIEWTILELEFT_HPP
#define KOKKOS_VIEWTILELEFT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct ViewTileLeftFast ;
struct ViewTileLeftSlow ;

template< class ValueType , unsigned N0 , unsigned N1 , bool B , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ValueType , void ,
                       LayoutTileLeft<N0,N1,B> ,
                       MemorySpace , MemoryTraits >
{ typedef typename if_c< B , ViewTileLeftFast , ViewTileLeftSlow >::type type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewTileLeftFast , void , void >
{
private:

  template< class DT , class DL , class DD , class DM >
  inline
  void allocate( View<DT,DL,DD,DM,ViewTileLeftFast> & dst , const std::string label )
  {
    typedef View<DT,DL,DD,DM,ViewTileLeftFast>  DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;

    dst.m_tracking.decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = (typename DstViewType::value_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::value_type) ,
                              sizeof(typename DstViewType::value_type) ,
                              dst.capacity() );

    ViewFill< DstViewType > init( dst , typename DstViewType::value_type() );
  }

public:

  template< class DT , class DL , class DD , class DM >
  inline
  ViewAssignment( View<DT,DL,DD,DM,ViewTileLeftFast> & dst ,
                  const typename enable_if< ViewTraits<DT,DL,DD,DM>::is_managed , std::string >::type & label ,
                  const size_t n0 ,
                  const size_t n1 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 ,
                  const size_t = 0 )
  {
    typedef View<DT,DL,DD,DM,ViewTileLeftFast>  DstViewType ;

    dst.m_shape.N0 = n0 ;
    dst.m_shape.N1 = n1 ;
    dst.m_tile_N0  = ( n0 + DstViewType::MASK_0 ) >> DstViewType::SHIFT_0 ;

    allocate( dst , label );
  }


  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,ViewTileLeftFast> & dst ,
                  const View<ST,SL,SD,SM,ViewTileLeftFast> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,ViewTileLeftFast> ,
                             typename View<ST,SL,SD,SM,ViewTileLeftFast>::HostMirror >::value
                  >::type * = 0 )
  {
    dst.m_shape   = src.m_shape ;
    dst.m_tile_N0 = src.m_tile_N0 ;
    allocate( dst , "mirror" );
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewTileLeftFast , ViewTileLeftFast, void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewTileLeftFast> & dst ,
                  const View<ST,SL,SD,SM,ViewTileLeftFast> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,ViewTileLeftFast> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    //typedef typename DstViewType::memory_space  memory_space ; // unused
    //typedef typename DstViewType::memory_traits memory_traits ; // unused

    dst.m_tracking.decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape, src.m_shape.N0 , src.m_shape.N1 );

    dst.m_tracking       = src.m_tracking ;
    dst.m_tile_N0       = src.m_tile_N0 ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    dst.m_tracking.increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  /** \brief  Deep copy data from compatible value type, layout, rank, and specialization.
   *          Check the dimensions and allocation lengths at runtime.
   */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline static
  void deep_copy( const View<DT,DL,DD,DM,Impl::ViewTileLeftFast> & dst ,
                  const View<ST,SL,SD,SM,Impl::ViewTileLeftFast> & src ,
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
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef ViewTraits<ST,SL,SD,SM> src_traits ;

    if ( dst.m_ptr_on_device != src.m_ptr_on_device ) {

      Impl::assert_shapes_are_equal( dst.m_shape , src.m_shape );

      const size_t n_dst = sizeof(typename dst_traits::value_type) * dst.capacity();
      const size_t n_src = sizeof(typename src_traits::value_type) * src.capacity();

      Impl::assert_counts_are_equal( n_dst , n_src );

      DeepCopy< typename dst_traits::memory_space ,
                typename src_traits::memory_space >( dst.m_ptr_on_device , src.m_ptr_on_device , n_dst );
    }
  }
};

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewDefault , ViewTileLeftFast, void >
{
  /** \brief Extracting a single tile from a tiled view */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewDefault> & dst ,
                  const View<ST,SL,SD,SM,ViewTileLeftFast> & src ,
                  const unsigned i0 ,
                  const typename enable_if<(
                    is_same< View<DT,DL,DD,DM,ViewDefault> ,
                             typename View<ST,SL,SD,SM,ViewTileLeftFast>::tile_type >::value
                  ), unsigned >::type i1 )
  {
    //typedef View<DT,DL,DD,DM,ViewDefault> DstViewType ; // unused
    //typedef typename DstViewType::shape_type    shape_type ; // unused
    //typedef typename DstViewType::memory_space  memory_space ; // unused
    //typedef typename DstViewType::memory_traits memory_traits ; // unused

    dst.m_tracking.decrement( dst.m_ptr_on_device );

    enum { N0 = SL::N0 };
    enum { N1 = SL::N1 };
    enum { SHIFT_0 = power_of_two<N0>::value };
    enum { MASK_0 = N0 - 1 };
    enum { SHIFT_1 = power_of_two<N1>::value };

    const unsigned NT0 = ( src.dimension_0() + MASK_0 ) >> SHIFT_0 ;

    dst.m_tracking      = src.m_tracking ;
    dst.m_ptr_on_device = src.m_ptr_on_device + (( i0 + i1 * NT0 ) << ( SHIFT_0 + SHIFT_1 ));

    dst.m_tracking.increment( dst.m_ptr_on_device );
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< class DataType , class Arg1Type , class Arg2Type , class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::ViewTileLeftFast >
  : public ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type >
{
private:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type > traits ;

  typedef Impl::ViewAssignment<Impl::ViewTileLeftFast> alloc ;

  typedef Impl::ViewAssignment<Impl::ViewTileLeftFast,
                               Impl::ViewTileLeftFast> assign ;

  typename traits::value_type * m_ptr_on_device ;
  typename traits::shape_type   m_shape ;
  unsigned                      m_tile_N0 ;
  Impl::ViewTracking< traits >  m_tracking ;

  typedef typename traits::array_layout layout ;

  enum { SHIFT_0 = Impl::power_of_two<layout::N0>::value };
  enum { SHIFT_1 = Impl::power_of_two<layout::N1>::value };
  enum { MASK_0  = layout::N0 - 1 };
  enum { MASK_1  = layout::N1 - 1 };

public:

  typedef Impl::ViewTileLeftFast specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  enum { Rank = 2 };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0) {}

  KOKKOS_INLINE_FUNCTION
  ~View() { m_tracking.decrement( m_ptr_on_device ); }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0) { (void)assign( *this , rhs ); }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { (void)assign( *this , rhs ); return *this ; }

  //------------------------------------
  // Array allocator and member access operator:

  View( const std::string & label , const size_t n0 , const size_t n1 )
    : m_ptr_on_device(0) { (void)alloc( *this , label , n0 , n1 ); }

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename traits::value_type & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );

      // Use care to insert necessary parentheses as the
      // shift operators have lower precedence than the arithmatic operators.

      return m_ptr_on_device[
        // ( ( Tile offset                               ) *  ( Tile size       ) )
         + ( ( (i0>>SHIFT_0) + m_tile_N0 * (i1>>SHIFT_1) ) << (SHIFT_0 + SHIFT_1) )
        // ( Offset within tile                       )
         + ( (i0 & MASK_0) + ((i1 & MASK_1)<<SHIFT_0) ) ] ;
    }

  //------------------------------------
  // Accept but ignore extra indices, they should be zero.

  template< typename iType0 , typename iType1 >
  KOKKOS_INLINE_FUNCTION
  typename traits::value_type &
    at( const iType0 & i0 , const iType1 & i1 , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );

      // Use care to insert necessary parentheses as the
      // shift operators have lower precedence than the arithmatic operators.

      return m_ptr_on_device[
        // ( ( Tile offset                               ) *  ( Tile size       ) )
         + ( ( (i0>>SHIFT_0) + m_tile_N0 * (i1>>SHIFT_1) ) << (SHIFT_0 + SHIFT_1) )
        // ( Offset within tile                       )
         + ( (i0 & MASK_0) + ((i1 & MASK_1)<<SHIFT_0) ) ] ;
    }

  //------------------------------------
  // Tile specialization specific declarations and functions:

  typedef View< typename traits::value_type [ layout::N0 ][ layout::N1 ] ,
                LayoutLeft ,
                typename traits::device_type ,
                MemoryUnmanaged >
    tile_type ;

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  KOKKOS_INLINE_FUNCTION
  size_t tiles_in_dimension_0() const { return m_tile_N0 ; }

  KOKKOS_INLINE_FUNCTION
  size_t tiles_in_dimension_1() const { return ( m_shape.N1 + MASK_1 ) >> SHIFT_1 ; }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t global_to_tile_index_0( const iType & global_i0 ) const
    { return global_i0 >> SHIFT_0 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t global_to_tile_index_1( const iType & global_i1 ) const
    { return global_i1 >> SHIFT_1 ; }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t global_to_local_tile_index_0( const iType & global_i0 ) const
    { return global_i0 & MASK_0 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t global_to_local_tile_index_1( const iType & global_i1 ) const
    { return global_i1 & MASK_1 ; }


  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const
  {
    return ( m_tile_N0 * ( ( m_shape.N1 + MASK_1 ) >> SHIFT_1 ) ) << ( SHIFT_0 + SHIFT_1 );
  }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWTILELEFT_HPP */

