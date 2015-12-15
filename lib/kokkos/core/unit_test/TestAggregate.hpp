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

#ifndef TEST_AGGREGATE_HPP
#define TEST_AGGREGATE_HPP

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

/*--------------------------------------------------------------------------*/

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace Test {

struct EmbedArray {};

struct ArrayProxyContiguous {};
struct ArrayProxyStrided {};

template< typename T , unsigned N = 0 , class Proxy = void >
struct Array ;

template< typename T >
struct Array<T,0,ArrayProxyContiguous>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * const value ;
  const unsigned count ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned n ) : value(v), count(n) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,0,Proxy> & rhs ) { return *this ; }
};

template< typename T , unsigned N >
struct Array<T,N,ArrayProxyContiguous>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T * const value ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned ) : value(v) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & rhs ) { return *this ; }
};

template< typename T , unsigned N >
struct Array<T,N,ArrayProxyStrided>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T * const value ;
  const unsigned stride ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned , unsigned s ) : value(v), stride(s) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & rhs ) { return *this ; }
};

template< typename T >
struct Array<T,0,ArrayProxyStrided>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * const value ;
  const unsigned count ;
  const unsigned stride ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned n , unsigned s ) : value(v), count(n), stride(s) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,0,Proxy> & rhs ) { return *this ; }
};

template< typename T >
struct Array<T,0,void>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * value ;
  const unsigned count ;

  KOKKOS_INLINE_FUNCTION
  Array() : value(0) , count(0) {}

  template< unsigned N , class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array( const Array<T,N,Proxy> & rhs ) : value(rhs.value), count(N) {}
};

template< typename T , unsigned N >
struct Array<T,N,void>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T value[N] ;

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & ) { return *this ; }
};

} // namespace Test

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template< typename T , unsigned N >
struct AnalyzeShape< Test::Array< T , N > >
  : public ShapeInsert< typename AnalyzeShape< T >::shape , N >::type
{
private:
  typedef AnalyzeShape< T > nested ;
public:

  typedef Test::EmbedArray specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_intrinsic_type   array_intrinsic_type[ N ];
  typedef Test::Array< T , N >          value_type ;
  typedef Test::Array< T , N >          type ;

  typedef const array_intrinsic_type  const_array_intrinsic_type ;
  typedef const value_type  const_value_type ;
  typedef const type        const_type ;

  typedef typename nested::non_const_array_intrinsic_type          non_const_array_intrinsic_type[ N ];
  typedef Test::Array< typename nested::non_const_value_type , N > non_const_value_type ;
  typedef Test::Array< typename nested::non_const_value_type , N > non_const_type ;
};

template< typename T >
struct AnalyzeShape< Test::Array< T , 0 > >
  : public ShapeInsert< typename AnalyzeShape< T >::shape , 0 >::type
{
private:
  typedef AnalyzeShape< T > nested ;
public:

  typedef Test::EmbedArray specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type * array_intrinsic_type ;
  typedef Test::Array< T , 0 >          value_type ;
  typedef Test::Array< T , 0 >          type ;

  typedef const array_intrinsic_type  const_array_intrinsic_type ;
  typedef const value_type  const_value_type ;
  typedef const type        const_type ;

  typedef typename nested::non_const_array_intrinsic_type  * non_const_array_intrinsic_type ;
  typedef Test::Array< typename nested::non_const_value_type , 0 > non_const_value_type ;
  typedef Test::Array< typename nested::non_const_value_type , 0 > non_const_type ;
};

/*--------------------------------------------------------------------------*/

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ValueType
                     , Test::EmbedArray
                     , LayoutLeft
                     , MemorySpace
                     , MemoryTraits >
{ typedef Test::EmbedArray type ; };

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ValueType
                     , Test::EmbedArray
                     , LayoutRight
                     , MemorySpace
                     , MemoryTraits >
{ typedef Test::EmbedArray type ; };

/*--------------------------------------------------------------------------*/

template<>
struct ViewAssignment< Test::EmbedArray , Test::EmbedArray , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Test::EmbedArray> & dst
                , const View<ST,SL,SD,SM,Test::EmbedArray> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    dst.m_tracker = src.m_tracker;
  }
};

template<>
struct ViewAssignment< ViewDefault , Test::EmbedArray , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,Test::EmbedArray>::array_type & dst
                , const View<ST,SL,SD,SM,Test::EmbedArray> & src
                )
  {
    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    dst.m_tracker = src.m_tracker;
  }
};


} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Test::EmbedArray >
  : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:

  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef Impl::ViewOffset< typename traits::shape_type ,
                            typename traits::array_layout > offset_map_type ;

  typedef Impl::ViewDataManagement< traits > view_data_management ;

  // traits::value_type = Test::Array< T , N >

  typename traits::value_type::value_type * m_ptr_on_device ;
  offset_map_type                           m_offset_map ;
  view_data_management                      m_management ;
  Impl::AllocationTracker                   m_tracker ;

public:

  typedef View< typename traits::array_intrinsic_type ,
                typename traits::array_layout ,
                typename traits::execution_space ,
                typename traits::memory_traits > array_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::execution_space ,
                typename traits::memory_traits > non_const_type ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::execution_space ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::host_mirror_space ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  enum { Rank = traits::rank - 1 };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_offset_map ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_offset_map.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_offset_map.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_offset_map.N2 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_offset_map.N3 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_offset_map.N4 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_offset_map.N5 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_offset_map.N6 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_offset_map.N7 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_offset_map.N0
           * m_offset_map.N1
           * m_offset_map.N2
           * m_offset_map.N3
           * m_offset_map.N4
           * m_offset_map.N5
           * m_offset_map.N6
           * m_offset_map.N7
           ;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_offset_map , i ); }

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() {}

  KOKKOS_INLINE_FUNCTION
  View()
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  { m_offset_map.assing(0,0,0,0,0,0,0,0); }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs )
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  {
    (void) Impl::ViewAssignment<
      typename traits::specialize ,
      typename traits::specialize >( *this , rhs );
  }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Construct or assign compatible view:

  template< class RT , class RL , class RD , class RM , class RS >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,RS> & rhs )
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  {
    (void) Impl::ViewAssignment<
      typename traits::specialize , RS >( *this , rhs );
  }

  template< class RT , class RL , class RD , class RM , class RS >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,RS> & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize , RS >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Allocation of a managed view with possible alignment padding.

  template< class AllocationProperties >
  explicit inline
  View( const AllocationProperties & prop ,
        const typename Impl::ViewAllocProp< traits , AllocationProperties >::size_type n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  {
    typedef Impl::ViewAllocProp< traits , AllocationProperties > Alloc ;

    typedef typename traits::memory_space  memory_space ;
    typedef typename traits::value_type::value_type   scalar_type ;

    m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
    m_offset_map.set_padding();

    m_tracker = memory_space::allocate_and_track( Alloc::label( prop ), sizeof(scalar_type) * m_offset_map.capacity() );

    m_ptr_on_device = reinterpret_cast<scalar_type *>(m_tracker.alloc_ptr());

    (void) Impl::ViewDefaultConstruct< typename traits::execution_space , scalar_type , Alloc::Initialize >( m_ptr_on_device , m_offset_map.capacity() );
  }

  //------------------------------------
  // Assign an unmanaged View from pointer, can be called in functors.
  // No alignment padding is performed.

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::value_type::value_type * ,
                      Impl::ViewError::user_pointer_constructor_requires_unmanaged >
    if_user_pointer_constructor ;

  View( typename if_user_pointer_constructor::type ptr ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  {
    m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
    m_ptr_on_device = if_user_pointer_constructor::select( ptr );
    m_management.set_unmanaged();
  }

  //------------------------------------
  // Assign unmanaged View to portion of Device shared memory

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::execution_space ,
                      Impl::ViewError::device_shmem_constructor_requires_unmanaged >
      if_device_shmem_constructor ;

  explicit KOKKOS_INLINE_FUNCTION
  View( typename if_device_shmem_constructor::type & dev ,
        const unsigned n0 = 0 ,
        const unsigned n1 = 0 ,
        const unsigned n2 = 0 ,
        const unsigned n3 = 0 ,
        const unsigned n4 = 0 ,
        const unsigned n5 = 0 ,
        const unsigned n6 = 0 ,
        const unsigned n7 = 0 )
    : m_ptr_on_device(0)
    , m_offset_map()
    , m_management()
    , m_tracker()
  {
    typedef typename traits::value_type::value_type   scalar_type ;

    enum { align = 8 };
    enum { mask  = align - 1 };

    typedef Impl::if_c< ! traits::is_managed ,
                        scalar_type * ,
                        Impl::ViewError::device_shmem_constructor_requires_unmanaged >
      if_device_shmem_pointer ;

    m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );

    // Select the first argument:
    m_ptr_on_device = if_device_shmem_pointer::select(
     (scalar_type *) dev.get_shmem( unsigned( sizeof(scalar_type) * m_offset_map.capacity() + unsigned(mask) ) & ~unsigned(mask) ) );
  }

  static inline
  unsigned shmem_size( const unsigned n0 = 0 ,
                       const unsigned n1 = 0 ,
                       const unsigned n2 = 0 ,
                       const unsigned n3 = 0 ,
                       const unsigned n4 = 0 ,
                       const unsigned n5 = 0 ,
                       const unsigned n6 = 0 ,
                       const unsigned n7 = 0 )
  {
    enum { align = 8 };
    enum { mask  = align - 1 };

    typedef typename traits::value_type::value_type   scalar_type ;

    offset_map_type offset_map ;

    offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );

    return unsigned( sizeof(scalar_type) * offset_map.capacity() + unsigned(mask) ) & ~unsigned(mask) ;
  }

  //------------------------------------
  // Is not allocated

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  // LayoutLeft, rank 2:

  typedef Test::Array< typename traits::value_type::value_type ,
                       traits::value_type::StaticLength ,
                       Test::ArrayProxyStrided > LeftValue ;

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_offset_map.N1 , m_offset_map.S0 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_offset_map.N1 , m_offset_map.S0 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_offset_map.N1 , m_offset_map.S0 );
    }

  //------------------------------------
  // LayoutRight, rank 2:

  typedef Test::Array< typename traits::value_type::value_type ,
                       traits::value_type::StaticLength ,
                       Test::ArrayProxyContiguous > RightValue ;

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_offset_map.N1 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_offset_map.N1 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_offset_map.N1 );
    }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    { m_offset_map.stride( s ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const
    { return m_offset_map.capacity(); }
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Test {

template< class DeviceType >
int TestViewAggregate()
{
  typedef Kokkos::View< Test::Array<double,32> * , DeviceType > a32_type ;
  typedef typename a32_type::array_type a32_base_type ;

  typedef Kokkos::View< Test::Array<double> * , DeviceType > a0_type ;
  typedef typename a0_type::array_type a0_base_type ;

  a32_type      a32("a32",100);
  a32_base_type a32_base ;

  a0_type       a0("a0",100,32);
  a0_base_type  a0_base ;

  a32_base = a32 ;
  a0_base = a0 ;


  return 0 ;
}

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#else /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

#include <impl/KokkosExp_ViewArray.hpp>

namespace Test {

template< class DeviceType >
void TestViewAggregate()
{
  typedef Kokkos::Array<double,32>  value_type ;

  typedef Kokkos::Experimental::Impl::
    ViewDataAnalysis< value_type * , Kokkos::LayoutLeft , value_type >
      analysis_1d ;

  static_assert( std::is_same< typename analysis_1d::specialize , Kokkos::Array<> >::value , "" );


  typedef Kokkos::ViewTraits< value_type ** , DeviceType > a32_traits ;
  typedef Kokkos::ViewTraits< typename a32_traits::array_scalar_type , DeviceType > flat_traits ;

  static_assert( std::is_same< typename a32_traits::specialize , Kokkos::Array<> >::value , "" );
  static_assert( std::is_same< typename a32_traits::value_type , value_type >::value , "" );
  static_assert( a32_traits::rank == 2 , "" );
  static_assert( a32_traits::rank_dynamic == 2 , "" );

  static_assert( std::is_same< typename flat_traits::specialize , void >::value , "" );
  static_assert( flat_traits::rank == 3 , "" );
  static_assert( flat_traits::rank_dynamic == 2 , "" );
  static_assert( flat_traits::dimension::N2 == 32 , "" );


  typedef Kokkos::View< Kokkos::Array<double,32> ** , DeviceType > a32_type ;

  typedef typename a32_type::array_type  a32_flat_type ;

  static_assert( std::is_same< typename a32_type::value_type , value_type >::value , "" );
  static_assert( std::is_same< typename a32_type::pointer_type , double * >::value , "" );
  static_assert( a32_type::Rank == 2 , "" );
  static_assert( a32_flat_type::Rank == 3 , "" );

  a32_type x("test",4,5);
  a32_flat_type y( x );

  ASSERT_EQ( x.extent(0) , 4 );
  ASSERT_EQ( x.extent(1) , 5 );
  ASSERT_EQ( y.extent(0) , 4 );
  ASSERT_EQ( y.extent(1) , 5 );
  ASSERT_EQ( y.extent(2) , 32 );
}

}

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* #ifndef TEST_AGGREGATE_HPP */
