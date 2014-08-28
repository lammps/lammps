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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Serial.hpp
/// \brief Declaration and definition of Kokkos::Serial device.

#ifndef KOKKOS_SERIAL_HPP
#define KOKKOS_SERIAL_HPP

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class Serial
/// \brief Kokkos device for non-parallel execution
///
/// A "device" represents a parallel execution model.  It tells Kokkos
/// how to parallelize the execution of kernels in a parallel_for or
/// parallel_reduce.  For example, the Threads device uses Pthreads or
/// C++11 threads on a CPU, the OpenMP device uses the OpenMP language
/// extensions, and the Cuda device uses NVIDIA's CUDA programming
/// model.  The Serial device executes "parallel" kernels
/// sequentially.  This is useful if you really do not want to use
/// threads, or if you want to explore different combinations of MPI
/// and shared-memory parallel programming models.
class Serial {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! The tag (what type of kokkos_object is this).
  typedef Impl::ExecutionSpaceTag  kokkos_tag ;

  //! The device type (same as this class).
  typedef Serial                device_type ;
  typedef Serial                execution_space ;
  //! The size_type typedef best suited for this device.
  typedef HostSpace::size_type  size_type ;
  //! This device's preferred memory space.
  typedef HostSpace             memory_space ;
  //! This device's preferred array layout.
  typedef LayoutRight           array_layout ;
  /// \brief This device's host mirror type.
  ///
  /// Serial is a host device, so the host mirror type is the same as
  /// the device type itself.
  typedef Serial                host_mirror_device_type ;

  /// \brief  Scratch memory space
  typedef ScratchMemorySpace< Kokkos::Serial >  scratch_memory_space ;

  //@}

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  ///
  /// For the Serial device, this method <i>always</i> returns false,
  /// because parallel_for or parallel_reduce with the Serial device
  /// always execute sequentially.
  inline static int in_parallel() { return false ; }

  /** \brief  Set the device in a "sleep" state.
   *
   * This function sets the device in a "sleep" state in which it is
   * not ready for work.  This may consume less resources than if the
   * device were in an "awake" state, but it may also take time to
   * bring the device from a sleep state to be ready for work.
   *
   * \return True if the device is in the "sleep" state, else false if
   *   the device is actively working and could not enter the "sleep"
   *   state.
   */
  static bool sleep();

  /// \brief Wake the device from the 'sleep' state so it is ready for work.
  ///
  /// \return True if the device is in the "ready" state, else "false"
  ///  if the device is actively working (which also means that it's
  ///  awake).
  static bool wake();

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void fence() {}

  static void initialize( unsigned threads_count = 1 ,
                          unsigned use_numa_count = 0 ,
                          unsigned use_cores_per_numa = 0 ,
                          bool allow_asynchronous_threadpool = false) {
    (void) threads_count;
    (void) use_numa_count;
    (void) use_cores_per_numa;
    (void) allow_asynchronous_threadpool;
  }

  static int is_initialized() { return 1 ; }

  //! Free any resources being consumed by the device.
  static void finalize() {}

  //! Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool detail = false );

  //--------------------------------------------------------------------------

  inline static int thread_pool_size( int = 0 ) { return 1 ; }
  KOKKOS_INLINE_FUNCTION static int thread_pool_rank() { return 0 ; }

  //--------------------------------------------------------------------------

  inline static unsigned hardware_thread_id() { return thread_pool_rank(); }
  inline static unsigned max_hardware_threads() { return thread_pool_size(0); }

  static inline int team_max()         { return thread_pool_size(1) ; }
  static inline int team_recommended() { return thread_pool_size(2); }

  //--------------------------------------------------------------------------

  static void * scratch_memory_resize( unsigned reduce_size , unsigned shared_size );

  //--------------------------------------------------------------------------
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<>
struct VerifyExecutionCanAccessMemorySpace
  < Kokkos::Serial::memory_space
  , Kokkos::Serial::scratch_memory_space
  >
{
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

class SerialTeamMember {
private:
  typedef Kokkos::ScratchMemorySpace< Kokkos::Serial > scratch_memory_space ;
  const scratch_memory_space  m_space ;
  const int                   m_league_rank ;
  const int                   m_league_size ;

  SerialTeamMember & operator = ( const SerialTeamMember & );

public:

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space & team_shmem() const { return m_space ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {}

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value , Type * const global_accum ) const
    {
      const Type tmp = global_accum ? *global_accum : Type(0) ;
      if ( global_accum ) { *global_accum += value ; }
      return tmp ;
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & ) const
    { return Type(0); }

  //----------------------------------------
  // Execution space specific:

  SerialTeamMember( int arg_league_rank
                  , int arg_league_size 
                  , int arg_shared_size
                  );
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

template < class WorkArgTag >
class TeamPolicy< Kokkos::Serial , WorkArgTag > {
private:

  const int m_league_size ;

public:

  typedef Impl::ExecutionPolicyTag   kokkos_tag ;      ///< Concept tag
  typedef Kokkos::Serial             execution_space ; ///< Execution space

  inline int team_size() const { return 1 ; }
  inline int league_size() const { return m_league_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space & , int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  TeamPolicy( int league_size_request , int /* team_size_request */ )
    : m_league_size( league_size_request )
    { }

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & ) { return 1 ; }

  typedef Impl::SerialTeamMember  member_type ;
};

} /* namespace Kokkos */

namespace Kokkos {
namespace Impl {

template< class FunctorType , typename IntType , unsigned P >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P >
                 , Kokkos::Serial
                 >
{
public:
  typedef Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P > Policy ;

  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
    {
      const typename Policy::member_type e = policy.end();
      for ( typename Policy::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i );
      }
    }
};

template< class FunctorType , typename IntType , unsigned P >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P >
                    , Kokkos::Serial
                    >
{
public:
  typedef Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P > Policy ;

  typedef ReduceAdapter< FunctorType >  Reduce ;
  typedef typename Reduce::pointer_type pointer_type ;

  template< class ViewType >
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const ViewType     & result
                , const typename enable_if<
                   ( is_view< ViewType >::value &&
                     is_same< typename ViewType::memory_space , HostSpace >::value
                   )>::type * = 0
                )
    {
      pointer_type result_ptr = result.ptr_on_device();

      if ( ! result_ptr ) {
        result_ptr = (pointer_type)
          Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );
      }

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename Policy::member_type e = policy.end();
      for ( typename Policy::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i , update );
      }

      Reduce::final( functor , result_ptr );
    }
};

template< class FunctorType , typename IntType , unsigned P >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P >
                  , Kokkos::Serial
                  >
{
public:
  typedef Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P > Policy ;

  typedef ReduceAdapter< FunctorType >  Reduce ;
  typedef typename Reduce::pointer_type pointer_type ;

  ParallelScan( const FunctorType  & functor
               , const Policy      & policy
               )
    {
      pointer_type result_ptr = (pointer_type)
        Kokkos::Serial::scratch_memory_resize( Reduce::value_size( functor ) , 0 );

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      const typename Policy::member_type e = policy.end();
      for ( typename Policy::member_type i = policy.begin() ; i < e ; ++i ) {
        functor( i , update , true );
      }

      Reduce::final( functor , result_ptr );
    }
};

//----------------------------------------------------------------------------

template< class FunctorType >
class ParallelFor< FunctorType , Kokkos::TeamPolicy< Kokkos::Serial , void > , Kokkos::Serial >
{
public:
  typedef Kokkos::TeamPolicy< Kokkos::Serial , void > Policy ;

  ParallelFor( const FunctorType & functor
             , const Policy      & policy )
    {
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );

      Kokkos::Serial::scratch_memory_resize( 0 , shared_size );

      for ( int ileague = 0 ; ileague < policy.league_size() ; ++ileague ) {
        functor( typename Policy::member_type(ileague,policy.league_size(),shared_size) );
      }
    }
};

template< class FunctorType >
class ParallelReduce< FunctorType , Kokkos::TeamPolicy< Kokkos::Serial , void > , Kokkos::Serial > {
public:

  typedef Kokkos::TeamPolicy< Kokkos::Serial , void > Policy ;
  typedef ReduceAdapter< FunctorType >  Reduce ;
  typedef typename Reduce::pointer_type pointer_type ;

  template< class ViewType >
  ParallelReduce( const FunctorType  & functor
                , const Policy       & policy
                , const ViewType     & result
                )
    {
      const int reduce_size = Reduce::value_size( functor );
      const int shared_size = FunctorTeamShmemSize< FunctorType >::value( functor , policy.team_size() );
      void * const scratch_reduce = Kokkos::Serial::scratch_memory_resize( reduce_size , shared_size );

      const pointer_type result_ptr =
        result.ptr_on_device() ? result.ptr_on_device() 
                               : (pointer_type) scratch_reduce ;

      typename Reduce::reference_type update = Reduce::init( functor , result_ptr );
      
      for ( int ileague = 0 ; ileague < policy.league_size() ; ++ileague ) {
        functor( typename Policy::member_type(ileague,policy.league_size(),shared_size) , update );
      }

      Reduce::final( functor , result_ptr );
    }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_SERIAL_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

