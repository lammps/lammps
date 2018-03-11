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

#ifndef KOKKOS_DYNAMIC_VIEW_HPP
#define KOKKOS_DYNAMIC_VIEW_HPP

#include <cstdio>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Experimental {

/** \brief Dynamic views are restricted to rank-one and no layout.
 *         Subviews are not allowed.
 */
template< typename DataType , typename ... P >
class DynamicView : public Kokkos::ViewTraits< DataType , P ... >
{
public:

  typedef Kokkos::ViewTraits< DataType , P ... >  traits ;

private:

  template< class , class ... > friend class DynamicView ;

  typedef Kokkos::Experimental::Impl::SharedAllocationTracker   track_type ;

  static_assert( traits::rank == 1 && traits::rank_dynamic == 1
               , "DynamicView must be rank-one" );

  static_assert( std::is_trivial< typename traits::value_type >::value &&
                 std::is_same< typename traits::specialize , void >::value &&
                 Kokkos::Impl::is_power_of_two
                   <sizeof(typename traits::value_type)>::value
               , "DynamicView must have trivial value_type and sizeof(value_type) is a power-of-two");


  template< class Space , bool = Kokkos::Impl::MemorySpaceAccess< Space , typename traits::memory_space >::accessible > struct verify_space
    { KOKKOS_FORCEINLINE_FUNCTION static void check() {} };

  template< class Space > struct verify_space<Space,false>
    { KOKKOS_FORCEINLINE_FUNCTION static void check()
        { Kokkos::abort("Kokkos::DynamicView ERROR: attempt to access inaccessible memory space"); };
    };

public:

  typedef Kokkos::MemoryPool< typename traits::device_type > memory_pool ;

private:

  memory_pool                    m_pool ;
  track_type                     m_track ;
  typename traits::value_type ** m_chunks ;
  unsigned                       m_chunk_shift ;
  unsigned                       m_chunk_mask ;
  unsigned                       m_chunk_max ;

public:

  //----------------------------------------------------------------------

  /** \brief  Compatible view of array of scalar types */
  typedef DynamicView< typename traits::data_type ,
                       typename traits::device_type >
    array_type ;

  /** \brief  Compatible view of const data type */
  typedef DynamicView< typename traits::const_data_type ,
                       typename traits::device_type >
    const_type ;

  /** \brief  Compatible view of non-const data type */
  typedef DynamicView< typename traits::non_const_data_type ,
                       typename traits::device_type >
    non_const_type ;

  /** \brief  Must be accessible everywhere */
  typedef DynamicView  HostMirror ;

  //----------------------------------------------------------------------

  enum { Rank = 1 };

  KOKKOS_INLINE_FUNCTION
  size_t size() const noexcept
    {
      uintptr_t n = 0 ;

      if ( Kokkos::Impl::MemorySpaceAccess
            < Kokkos::Impl::ActiveExecutionMemorySpace
            , typename traits::memory_space
            >::accessible ) {
        n = *reinterpret_cast<const uintptr_t*>( m_chunks + m_chunk_max );
      }
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      else {
        Kokkos::Impl::DeepCopy< Kokkos::HostSpace
                              , typename traits::memory_space
                              , Kokkos::HostSpace::execution_space >
          ( & n
          , reinterpret_cast<const uintptr_t*>( m_chunks + m_chunk_max )
          , sizeof(uintptr_t) );
      }
#endif
      return n << m_chunk_shift ;
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t extent( const iType & r ) const
    { return r == 0 ? size() : 1 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  size_t extent_int( const iType & r ) const
    { return r == 0 ? size() : 1 ; }

  KOKKOS_INLINE_FUNCTION size_t dimension_0() const { return size(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return 0 ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const { *s = 0 ; }

  //----------------------------------------------------------------------
  // Range span is the span which contains all members.

  typedef typename traits::value_type &  reference_type ;
  typedef typename traits::value_type *  pointer_type ;

  enum { reference_type_is_lvalue_reference = std::is_lvalue_reference< reference_type >::value };

  KOKKOS_INLINE_FUNCTION constexpr bool   span_is_contiguous() const { return false ; }
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const { return 0 ; }

  //----------------------------------------

  template< typename I0 , class ... Args >
  KOKKOS_INLINE_FUNCTION
  reference_type operator()( const I0 & i0 , const Args & ... args ) const
    {
      static_assert( Kokkos::Impl::are_integral<I0,Args...>::value
                   , "Indices must be integral type" );

      DynamicView::template verify_space< Kokkos::Impl::ActiveExecutionMemorySpace >::check();

      // Which chunk is being indexed.
      const uintptr_t ic = uintptr_t( i0 >> m_chunk_shift );

      typename traits::value_type * volatile * const ch = m_chunks + ic ;

      // Do bounds checking if enabled or if the chunk pointer is zero.
      // If not bounds checking then we assume a non-zero pointer is valid.

#if ! defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )
      if ( 0 == *ch )
#endif
      {
        // Verify that allocation of the requested chunk in in progress.

        // The allocated chunk counter is m_chunks[ m_chunk_max ]
        const uintptr_t n =
          *reinterpret_cast<uintptr_t volatile *>( m_chunks + m_chunk_max );

        if ( n <= ic ) {
          Kokkos::abort("Kokkos::DynamicView array bounds error");
        }

        // Allocation of this chunk is in progress
        // so wait for allocation to complete.
        while ( 0 == *ch );
      }

      return (*ch)[ i0 & m_chunk_mask ];
    }

  //----------------------------------------
  /** \brief  Resizing in parallel only increases the array size,
   *          never decrease.
   */
  KOKKOS_INLINE_FUNCTION
  void resize_parallel( size_t n ) const
    {
      typedef typename traits::value_type value_type ;

      DynamicView::template verify_space< Kokkos::Impl::ActiveExecutionMemorySpace >::check();

      const uintptr_t NC = ( n + m_chunk_mask ) >> m_chunk_shift ;

      if ( m_chunk_max < NC ) {
#if defined( KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK )
        printf("DynamicView::resize_parallel(%lu) m_chunk_max(%u) NC(%lu)\n"
              , n , m_chunk_max , NC );
#endif
        Kokkos::abort("DynamicView::resize_parallel exceeded maximum size");
      }

      typename traits::value_type * volatile * const ch = m_chunks ;

      // The allocated chunk counter is m_chunks[ m_chunk_max ]
      uintptr_t volatile * const pc =
        reinterpret_cast<uintptr_t volatile*>( m_chunks + m_chunk_max );

      // Potentially concurrent iteration of allocation to the required size.

      for ( uintptr_t jc = *pc ; jc < NC ; ) {

        // Claim the 'jc' chunk to-be-allocated index

        const uintptr_t jc_try = jc ;

        // Jump iteration to the chunk counter.

        jc = atomic_compare_exchange( pc , jc_try , jc_try + 1 );

        if ( jc_try == jc ) {

          ch[jc_try] = reinterpret_cast<value_type*>(
            m_pool.allocate( sizeof(value_type) << m_chunk_shift ));

          if ( 0 == ch[jc_try] ) {
            Kokkos::abort("DynamicView::resize_parallel exhausted memory pool");
          }

          Kokkos::memory_fence();
        }
      }
    }

  /** \brief  Resizing in serial can grow or shrink the array size, */
  template< typename IntType >
  inline
  typename std::enable_if
    < std::is_integral<IntType>::value &&
      Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace
                                     , typename traits::memory_space
                                     >::accessible
    >::type
  resize_serial( IntType const & n )
    {
      typedef typename traits::value_type value_type ;
      typedef value_type * pointer_type ;

      const uintptr_t NC = ( n + m_chunk_mask ) >> m_chunk_shift ;

      if ( m_chunk_max < NC ) {
        Kokkos::abort("DynamicView::resize_serial exceeded maximum size");
      }

      uintptr_t * const pc =
        reinterpret_cast<uintptr_t*>( m_chunks + m_chunk_max );

      if ( *pc < NC ) {
        while ( *pc < NC ) {
          m_chunks[*pc] = reinterpret_cast<pointer_type>
            ( m_pool.allocate( sizeof(value_type) << m_chunk_shift ) );
          ++*pc ;
        }
      }
      else {
        while ( NC + 1 <= *pc ) {
          --*pc ;
          m_pool.deallocate( m_chunks[*pc]
                           , sizeof(value_type) << m_chunk_shift );
          m_chunks[*pc] = 0 ;
        }
      }
    }

  //----------------------------------------

  struct ResizeSerial {
    memory_pool                    m_pool ;
    typename traits::value_type ** m_chunks ;
    uintptr_t                    * m_pc ;
    uintptr_t                      m_nc ;
    unsigned                       m_chunk_shift ;

    KOKKOS_INLINE_FUNCTION
    void operator()( int ) const
      {
        typedef typename traits::value_type value_type ;
        typedef value_type * pointer_type ;

        if ( *m_pc < m_nc ) {
          while ( *m_pc < m_nc ) {
            m_chunks[*m_pc] = reinterpret_cast<pointer_type>
              ( m_pool.allocate( sizeof(value_type) << m_chunk_shift ) );
            ++*m_pc ;
          }
        }
        else {
          while ( m_nc + 1 <= *m_pc ) {
            --*m_pc ;
            m_pool.deallocate( m_chunks[*m_pc]
                             , sizeof(value_type) << m_chunk_shift );
            m_chunks[*m_pc] = 0 ;
          }
        }
      }

    ResizeSerial( memory_pool            const & arg_pool
                , typename traits::value_type ** arg_chunks
                , uintptr_t                    * arg_pc
                , uintptr_t                      arg_nc
                , unsigned                       arg_chunk_shift
                )
      : m_pool( arg_pool )
      , m_chunks( arg_chunks )
      , m_pc( arg_pc )
      , m_nc( arg_nc )
      , m_chunk_shift( arg_chunk_shift )
      {}
  };

  template< typename IntType >
  inline
  typename std::enable_if
    < std::is_integral<IntType>::value &&
      ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace
                                       , typename traits::memory_space
                                       >::accessible
    >::type
  resize_serial( IntType const & n )
    {
      const uintptr_t NC = ( n + m_chunk_mask ) >> m_chunk_shift ;

      if ( m_chunk_max < NC ) {
        Kokkos::abort("DynamicView::resize_serial exceeded maximum size");
      }

      // Must dispatch kernel

      typedef Kokkos::RangePolicy< typename traits::execution_space > Range ;

      uintptr_t * const pc =
        reinterpret_cast<uintptr_t*>( m_chunks + m_chunk_max );

      Kokkos::Impl::ParallelFor<ResizeSerial,Range>
        closure( ResizeSerial( m_pool, m_chunks, pc, NC, m_chunk_shift )
               , Range(0,1) );

      closure.execute();

      traits::execution_space::fence();
    }

  //----------------------------------------------------------------------

  ~DynamicView() = default ;
  DynamicView() = default ;
  DynamicView( DynamicView && ) = default ;
  DynamicView( const DynamicView & ) = default ;
  DynamicView & operator = ( DynamicView && ) = default ;
  DynamicView & operator = ( const DynamicView & ) = default ;

  template< class RT , class ... RP >
  DynamicView( const DynamicView<RT,RP...> & rhs )
    : m_pool( rhs.m_pool )
    , m_track( rhs.m_track )
    , m_chunks( (typename traits::value_type **) rhs.m_chunks )
    , m_chunk_shift( rhs.m_chunk_shift )
    , m_chunk_mask( rhs.m_chunk_mask )
    , m_chunk_max( rhs.m_chunk_max )
    {
      typedef typename DynamicView<RT,RP...>::traits  SrcTraits ;
      typedef Kokkos::Impl::ViewMapping< traits , SrcTraits , void >  Mapping ;
      static_assert( Mapping::is_assignable , "Incompatible DynamicView copy construction" );
    }

  //----------------------------------------------------------------------

  struct Destroy {
    memory_pool                    m_pool ;
    typename traits::value_type ** m_chunks ;
    unsigned                       m_chunk_max ;
    bool                           m_destroy ;

    // Initialize or destroy array of chunk pointers.
    // Two entries beyond the max chunks are allocation counters.

    KOKKOS_INLINE_FUNCTION
    void operator()( unsigned i ) const
      {
        if ( m_destroy && i < m_chunk_max && 0 != m_chunks[i] ) {
          m_pool.deallocate( m_chunks[i] , m_pool.min_block_size() );
        }
        m_chunks[i] = 0 ;
      }

    void execute( bool arg_destroy )
      {
        typedef Kokkos::RangePolicy< typename traits::execution_space > Range ;

        m_destroy = arg_destroy ;

        Kokkos::Impl::ParallelFor<Destroy,Range>
          closure( *this , Range(0, m_chunk_max + 1) );

        closure.execute();

        traits::execution_space::fence();
      }

    void construct_shared_allocation()
      { execute( false ); }

    void destroy_shared_allocation()
      { execute( true ); }

    Destroy() = default ;
    Destroy( Destroy && ) = default ;
    Destroy( const Destroy & ) = default ;
    Destroy & operator = ( Destroy && ) = default ;
    Destroy & operator = ( const Destroy & ) = default ;

    Destroy( const memory_pool & arg_pool
           , typename traits::value_type ** arg_chunk
           , const unsigned arg_chunk_max )
     : m_pool( arg_pool )
     , m_chunks( arg_chunk )
     , m_chunk_max( arg_chunk_max )
     , m_destroy( false )
     {}
  };


  /**\brief  Allocation constructor
   *
   *  Memory is allocated in chunks from the memory pool.
   *  The chunk size conforms to the memory pool's chunk size.
   *  A maximum size is required in order to allocate a
   *  chunk-pointer array.
   */
  explicit inline
  DynamicView( const std::string & arg_label
             , const memory_pool & arg_pool
             , const size_t        arg_size_max )
    : m_pool( arg_pool )
    , m_track()
    , m_chunks(0)
    // The memory pool chunk is guaranteed to be a power of two
    , m_chunk_shift(
        Kokkos::Impl::integral_power_of_two(
          m_pool.min_block_size()/sizeof(typename traits::value_type)) )
    , m_chunk_mask( ( 1 << m_chunk_shift ) - 1 )
    , m_chunk_max( ( arg_size_max + m_chunk_mask ) >> m_chunk_shift )
    {
      // A functor to deallocate all of the chunks upon final destruction

      typedef typename traits::memory_space  memory_space ;
      typedef Kokkos::Experimental::Impl::SharedAllocationRecord< memory_space , Destroy > record_type ;

      // Allocate chunk pointers and allocation counter
      record_type * const record =
        record_type::allocate( memory_space()
                             , arg_label
                             , ( sizeof(pointer_type) * ( m_chunk_max + 1 ) ) );

      m_chunks = reinterpret_cast<pointer_type*>( record->data() );

      record->m_destroy = Destroy( m_pool , m_chunks , m_chunk_max );

      // Initialize to zero

      record->m_destroy.construct_shared_allocation();

      m_track.assign_allocated_record_to_uninitialized( record );
    }
};

} // namespace Experimental
} // namespace Kokkos

namespace Kokkos {
namespace Experimental {

template< class T , class ... P >
inline
typename Kokkos::Experimental::DynamicView<T,P...>::HostMirror
create_mirror_view( const Kokkos::Experimental::DynamicView<T,P...> & src )
{
  return src ;
}

template< class T , class ... DP , class ... SP >
inline
void deep_copy( const View<T,DP...> & dst
              , const DynamicView<T,SP...> & src
              )
{
  typedef View<T,DP...>        dst_type ;
  typedef DynamicView<T,SP...> src_type ;

  typedef typename ViewTraits<T,DP...>::execution_space  dst_execution_space ;
  typedef typename ViewTraits<T,SP...>::memory_space     src_memory_space ;

  enum { DstExecCanAccessSrc =
   Kokkos::Impl::SpaceAccessibility< dst_execution_space , src_memory_space >::accessible };

  if ( DstExecCanAccessSrc ) {
    // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap< dst_type , src_type >( dst , src );
  }
  else {
    Kokkos::Impl::throw_runtime_exception("deep_copy given views that would require a temporary allocation");
  }
}

template< class T , class ... DP , class ... SP >
inline
void deep_copy( const DynamicView<T,DP...> & dst
              , const View<T,SP...> & src
              )
{
  typedef DynamicView<T,SP...> dst_type ;
  typedef View<T,DP...>        src_type ;

  typedef typename ViewTraits<T,DP...>::execution_space  dst_execution_space ;
  typedef typename ViewTraits<T,SP...>::memory_space     src_memory_space ;

  enum { DstExecCanAccessSrc =
   Kokkos::Impl::SpaceAccessibility< dst_execution_space , src_memory_space >::accessible };

  if ( DstExecCanAccessSrc ) {
    // Copying data between views in accessible memory spaces and either non-contiguous or incompatible shape.
    Kokkos::Impl::ViewRemap< dst_type , src_type >( dst , src );
  }
  else {
    Kokkos::Impl::throw_runtime_exception("deep_copy given views that would require a temporary allocation");
  }
}

} // namespace Experimental
} // namespace Kokkos

#endif /* #ifndef KOKKOS_DYNAMIC_VIEW_HPP */

