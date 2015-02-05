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

#include <Kokkos_ExecPolicy.hpp>
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

template< class ExecSpace , class Type , bool Initialize >
struct ViewDefaultConstruct
{ ViewDefaultConstruct( Type * , size_t ) {} };


/** \brief  ViewDataHandle provides the type of the 'data handle' which the view
 *          uses to access data with the [] operator. It also provides
 *          an allocate function and a function to extract a raw ptr from the
 *          data handle. ViewDataHandle also defines an enum ReferenceAble which
 *          specifies whether references/pointers to elements can be taken and a
 *          'return_type' which is what the view operators will give back.
 *          Specialisation of this object allows three things depending
 *          on ViewTraits and compiler options:
 *          (i)   Use special allocator (e.g. huge pages/small pages and pinned memory)
 *          (ii)  Use special data handle type (e.g. add Cuda Texture Object)
 *          (iii) Use special access intrinsics (e.g. texture fetch and non-caching loads)
 */
template< class StaticViewTraits , class Enable = void >
struct ViewDataHandle {

  enum { ReturnTypeIsReference = true };

  typedef typename StaticViewTraits::value_type * handle_type;
  typedef typename StaticViewTraits::value_type & return_type;
};

template< class StaticViewTraits , class Enable = void >
class ViewDataManagement : public ViewDataHandle< StaticViewTraits > {
private:

  template< class , class > friend class ViewDataManagement ;

  struct PotentiallyManaged  {};
  struct StaticallyUnmanaged {};

  /* Statically unmanaged if traits or not executing in host-accessible memory space */
  typedef typename
    Impl::if_c< StaticViewTraits::is_managed &&
                Impl::is_same< Kokkos::HostSpace
                             , Kokkos::Impl::ActiveExecutionMemorySpace >::value
              , PotentiallyManaged
              , StaticallyUnmanaged
              >::type StaticManagementTag ;

  enum { Unmanaged     = 0x01
       , Noncontiguous = 0x02
       };

  enum { DefaultTraits = Impl::is_same< StaticManagementTag , StaticallyUnmanaged >::value ? Unmanaged : 0 };

  unsigned m_traits ; ///< Runtime traits


  template< class T >
  inline static
  unsigned assign( const ViewDataManagement<T> & rhs , const PotentiallyManaged & )
    { return rhs.m_traits | ( rhs.is_managed() && Kokkos::HostSpace::in_parallel() ? unsigned(Unmanaged) : 0u ); }

  template< class T >
  KOKKOS_INLINE_FUNCTION static
  unsigned assign( const ViewDataManagement<T> & rhs , const StaticallyUnmanaged & )
    { return rhs.m_traits | Unmanaged ; }

  inline
  void increment( const void * ptr , const PotentiallyManaged & ) const
    { if ( is_managed() ) StaticViewTraits::memory_space::increment( ptr ); }
  
  inline
  void decrement( const void * ptr , const PotentiallyManaged & ) const
    { if ( is_managed() ) StaticViewTraits::memory_space::decrement( ptr ); }
  
  KOKKOS_INLINE_FUNCTION
  void increment( const void * , const StaticallyUnmanaged & ) const {}
  
  KOKKOS_INLINE_FUNCTION
  void decrement( const void * , const StaticallyUnmanaged & ) const {}

public:

  typedef typename ViewDataHandle< StaticViewTraits >::handle_type handle_type;

  KOKKOS_INLINE_FUNCTION
  ViewDataManagement() : m_traits( DefaultTraits ) {}

  KOKKOS_INLINE_FUNCTION
  ViewDataManagement( const ViewDataManagement & rhs )
    : m_traits( assign( rhs , StaticManagementTag() ) ) {}

  KOKKOS_INLINE_FUNCTION
  ViewDataManagement & operator = ( const ViewDataManagement & rhs )
    { m_traits = assign( rhs , StaticManagementTag() ); return *this ; }

  template< class SVT >
  KOKKOS_INLINE_FUNCTION
  ViewDataManagement( const ViewDataManagement<SVT> & rhs )
    : m_traits( assign( rhs , StaticManagementTag() ) ) {}

  template< class SVT >
  KOKKOS_INLINE_FUNCTION
  ViewDataManagement & operator = ( const ViewDataManagement<SVT> & rhs )
    { m_traits = assign( rhs , StaticManagementTag() ); return *this ; }

  KOKKOS_INLINE_FUNCTION
  bool is_managed() const { return ! ( m_traits & Unmanaged ); }

  KOKKOS_INLINE_FUNCTION
  bool is_contiguous() const { return ! ( m_traits & Noncontiguous ); }

  KOKKOS_INLINE_FUNCTION
  void set_unmanaged() { m_traits |= Unmanaged ; }

  KOKKOS_INLINE_FUNCTION
  void set_noncontiguous() { m_traits |= Noncontiguous ; }


  KOKKOS_INLINE_FUNCTION
  void increment( handle_type handle ) const
    { increment( ( typename StaticViewTraits::value_type *) handle , StaticManagementTag() ); }

  KOKKOS_INLINE_FUNCTION
  void decrement( handle_type handle ) const
    { decrement( ( typename StaticViewTraits::value_type *) handle , StaticManagementTag() ); }


  KOKKOS_INLINE_FUNCTION
  void increment( const void * ptr ) const
    { increment( ptr , StaticManagementTag() ); }

  KOKKOS_INLINE_FUNCTION
  void decrement( const void * ptr ) const
    { decrement( ptr , StaticManagementTag() ); }


  template< bool Initialize >
  static
  handle_type allocate( const std::string & label
                      , const Impl::ViewOffset< typename StaticViewTraits::shape_type
                                              , typename StaticViewTraits::array_layout > & offset_map )
    {
      typedef typename StaticViewTraits::execution_space  execution_space ;
      typedef typename StaticViewTraits::memory_space     memory_space ;
      typedef typename StaticViewTraits::value_type       value_type ;

      const size_t count = offset_map.capacity();

      value_type * ptr = (value_type*) memory_space::allocate( label , sizeof(value_type) * count );

        // Default construct within the view's execution space.
      (void) ViewDefaultConstruct< execution_space , value_type , Initialize >( ptr , count );

      return typename ViewDataHandle< StaticViewTraits >::handle_type(ptr);
    }
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
  typedef typename OutputView::size_type   size_type ;

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
      typedef typename OutputView::execution_space execution_space ;
      Kokkos::RangePolicy< execution_space > range( 0 , n0 );
      parallel_for( range , *this );
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

template< class ExecSpace , class Type >
struct ViewDefaultConstruct< ExecSpace , Type , true >
{
  Type * const m_ptr ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename ExecSpace::size_type i ) const
    { new( m_ptr + i ) Type(); }

  ViewDefaultConstruct( Type * pointer , size_t capacity )
    : m_ptr( pointer )
    {
      Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
      parallel_for( range , *this );
      ExecSpace::fence();
    }
};

template< class OutputView , unsigned Rank = OutputView::Rank >
struct ViewFill
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::size_type         size_type ;

  const OutputView output ;
  const_value_type input ;

  ViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      typedef typename OutputView::execution_space execution_space ;
      Kokkos::RangePolicy< execution_space > range( 0 , output.dimension_0() );
      parallel_for( range , *this );
      execution_space::fence();
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

namespace Kokkos {

struct ViewAllocateWithoutInitializing {

  const std::string label ;

  ViewAllocateWithoutInitializing() : label() {}
  ViewAllocateWithoutInitializing( const std::string & arg_label ) : label( arg_label ) {}
  ViewAllocateWithoutInitializing( const char * const  arg_label ) : label( arg_label ) {}
};

struct ViewAllocate {

  const std::string  label ;

  ViewAllocate() : label() {}
  ViewAllocate( const std::string & arg_label ) : label( arg_label ) {}
  ViewAllocate( const char * const  arg_label ) : label( arg_label ) {}
};

}

namespace Kokkos {
namespace Impl {

template< class Traits , class AllocationProperties , class Enable = void >
struct ViewAllocProp : public Kokkos::Impl::false_type {};

template< class Traits >
struct ViewAllocProp< Traits , Kokkos::ViewAllocate
  , typename Kokkos::Impl::enable_if<(
      Traits::is_managed && ! Kokkos::Impl::is_const< typename Traits::value_type >::value
    )>::type >
  : public Kokkos::Impl::true_type
{
  typedef size_t               size_type ;
  typedef const ViewAllocate & property_type ;

  enum { Initialize = true };
  enum { AllowPadding = false };

  inline
  static const std::string & label( property_type p ) { return p.label ; }
};

template< class Traits >
struct ViewAllocProp< Traits , std::string
  , typename Kokkos::Impl::enable_if<(
      Traits::is_managed && ! Kokkos::Impl::is_const< typename Traits::value_type >::value
    )>::type >
  : public Kokkos::Impl::true_type
{
  typedef size_t              size_type ;
  typedef const std::string & property_type ;

  enum { Initialize = true };
  enum { AllowPadding = false };

  inline
  static const std::string & label( property_type s ) { return s ; }
};

template< class Traits , unsigned N >
struct ViewAllocProp< Traits , char[N]
  , typename Kokkos::Impl::enable_if<(
      Traits::is_managed && ! Kokkos::Impl::is_const< typename Traits::value_type >::value
    )>::type >
  : public Kokkos::Impl::true_type
{
private:
  typedef char label_type[N] ;
public:

  typedef size_t             size_type ;
  typedef const label_type & property_type ;

  enum { Initialize = true };
  enum { AllowPadding = false };

  inline
  static std::string label( property_type s ) { return std::string(s) ; }
};

template< class Traits >
struct ViewAllocProp< Traits , Kokkos::ViewAllocateWithoutInitializing
  , typename Kokkos::Impl::enable_if<(
      Traits::is_managed && ! Kokkos::Impl::is_const< typename Traits::value_type >::value
    )>::type >
  : public Kokkos::Impl::true_type
{
  typedef size_t size_type ;
  typedef const Kokkos::ViewAllocateWithoutInitializing & property_type ;

  enum { Initialize = false };
  enum { AllowPadding = false };

  inline
  static std::string label( property_type s ) { return s.label ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class Traits , class PointerProperties , class Enable = void >
struct ViewRawPointerProp : public Kokkos::Impl::false_type {};

template< class Traits , typename T >
struct ViewRawPointerProp< Traits , T ,
  typename Kokkos::Impl::enable_if<(
    Impl::is_same< T , typename Traits::value_type >::value ||
    Impl::is_same< T , typename Traits::non_const_value_type >::value
  )>::type >
  : public Kokkos::Impl::true_type
{
  typedef size_t size_type ; 
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWSUPPORT_HPP */


