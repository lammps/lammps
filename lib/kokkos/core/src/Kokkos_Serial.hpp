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
#include <impl/Kokkos_FunctorAdapter.hpp>

#if defined( KOKKOS_HAVE_SERIAL )

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

  //! Tag this class as an execution space:
  typedef Serial                execution_space ;
  //! The size_type typedef best suited for this device.
  typedef HostSpace::size_type  size_type ;
  //! This device's preferred memory space.
  typedef HostSpace             memory_space ;
  //! This execution space preferred device_type
  typedef Kokkos::Device<execution_space,memory_space> device_type;

  //! This device's preferred array layout.
  typedef LayoutRight           array_layout ;

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

    // Init the array of locks used for arbitrarily sized atomics
    Impl::init_lock_array_host_space();

  }

  static int is_initialized() { return 1 ; }

  //! Free any resources being consumed by the device.
  static void finalize() {}

  //! Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool /* detail */ = false ) {}

  //--------------------------------------------------------------------------

  inline static int thread_pool_size( int = 0 ) { return 1 ; }
  KOKKOS_INLINE_FUNCTION static int thread_pool_rank() { return 0 ; }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION static unsigned hardware_thread_id() { return thread_pool_rank(); }
  inline static unsigned max_hardware_threads() { return thread_pool_size(0); }

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
  enum { value = true };
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

namespace SerialImpl {

struct Sentinel {

  void *   m_scratch ;
  unsigned m_reduce_end ;
  unsigned m_shared_end ;

  Sentinel();
  ~Sentinel();
  static Sentinel & singleton();
};

inline
unsigned align( unsigned n );
}
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

  template<class ValueType>
  KOKKOS_INLINE_FUNCTION
  void team_broadcast(const ValueType& , const int& ) const {}

  template< class ValueType, class JoinOp >
  KOKKOS_INLINE_FUNCTION
  ValueType team_reduce( const ValueType & value , const JoinOp & ) const
    {
      return value ;
    }

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


/*
 * < Kokkos::Serial , WorkArgTag >
 * < WorkArgTag , Impl::enable_if< Impl::is_same< Kokkos::Serial , Kokkos::DefaultExecutionSpace >::value >::type >
 *
 */
template< class Arg0 , class Arg1 >
class TeamPolicy< Arg0 , Arg1 , Kokkos::Serial >
{
private:

  const int m_league_size ;
  const int m_scratch_size ;

public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicy      execution_policy ;

  //! Execution space of this execution policy:
  typedef Kokkos::Serial  execution_space ;

  typedef typename
    Impl::if_c< ! Impl::is_same< Kokkos::Serial , Arg0 >::value , Arg0 , Arg1 >::type
      work_tag ;

  //----------------------------------------

  template< class FunctorType >
  static
  int team_size_max( const FunctorType & ) { return 1 ; }

  template< class FunctorType >
  static
  int team_size_recommended( const FunctorType & ) { return 1 ; }

  template< class FunctorType >
  static
  int team_size_recommended( const FunctorType & , const int& ) { return 1 ; }

  //----------------------------------------

  inline int team_size() const { return 1 ; }
  inline int league_size() const { return m_league_size ; }
  inline size_t scratch_size() const { return m_scratch_size ; }

  /** \brief  Specify league size, request team size */
  TeamPolicy( execution_space &
            , int league_size_request
            , int /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size( league_size_request )
    , m_scratch_size ( 0 )
    {}

  TeamPolicy( execution_space &
            , int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size( league_size_request )
    , m_scratch_size ( 0 )
    {}

  TeamPolicy( int league_size_request
            , int /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size( league_size_request )
    , m_scratch_size ( 0 )
    {}

  TeamPolicy( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_league_size( league_size_request )
    , m_scratch_size ( 0 )
    {}

  template<class MemorySpace>
  TeamPolicy( int league_size_request
            , int /* team_size_request */
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size(league_size_request)
    , m_scratch_size(scratch_request.total(1))
    {}


  template<class MemorySpace>
  TeamPolicy( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , const Experimental::TeamScratchRequest<MemorySpace> & scratch_request )
    : m_league_size(league_size_request)
    , m_scratch_size(scratch_request.total(1))
    {}

  typedef Impl::SerialTeamMember  member_type ;
};

} /* namespace Kokkos */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Parallel patterns for Kokkos::Serial with RangePolicy */

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelFor< FunctorType
                 , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial >
                 >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial > Policy ;

  const FunctorType m_functor ;
  const Policy      m_policy ;

  template< class TagType >
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec() const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i );
      }
    }

  template< class TagType >
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec() const
    {
      const TagType t{} ;
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i );
      }
    }

public:

  inline
  void execute() const
    { this-> template exec< typename Policy::work_tag >(); }

  inline
  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    {}
};

/*--------------------------------------------------------------------------*/

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial >
                    >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial > Policy ;
  typedef typename Policy::work_tag                                  WorkTag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , WorkTag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const pointer_type  m_result_ptr ;


  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      reference_type update = ValueInit::init( m_functor , ptr );

      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i , update );
      }

      Kokkos::Impl::FunctorFinal< FunctorType , TagType >::
        final( m_functor , ptr );
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      const TagType t{} ;
      reference_type update = ValueInit::init( m_functor , ptr );

      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i , update );
      }

      Kokkos::Impl::FunctorFinal< FunctorType , TagType >::
        final( m_functor , ptr );
    }

public:

  inline
  void execute() const
    {
      pointer_type ptr = (pointer_type) Kokkos::Serial::scratch_memory_resize
           ( ValueTraits::value_size( m_functor ) , 0 );

      this-> template exec< WorkTag >( m_result_ptr ? m_result_ptr : ptr );
    }

  template< class ViewType >
  ParallelReduce( const FunctorType & arg_functor
                , const Policy      & arg_policy
                , const ViewType    & arg_result )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_result_ptr( arg_result.ptr_on_device() )
    {
      static_assert( Kokkos::is_view< ViewType >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View" );

      static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View in HostSpace" );
    }
};

/*--------------------------------------------------------------------------*/

template< class FunctorType , class Arg0 , class Arg1 , class Arg2 >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial >
                  >
{
private:

  typedef Kokkos::RangePolicy< Arg0 , Arg1 , Arg2 , Kokkos::Serial > Policy ;
  typedef typename Policy::work_tag                                  WorkTag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , WorkTag >  ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag >  ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      reference_type update = ValueInit::init( m_functor , ptr );

      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i , update , true );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      const TagType t{} ;
      reference_type update = ValueInit::init( m_functor , ptr );

      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i , update , true );
      }
    }

public:

  inline
  void execute() const
    {
      pointer_type ptr = (pointer_type)
        Kokkos::Serial::scratch_memory_resize( ValueTraits::value_size( m_functor ) , 0 );
      this-> template exec< WorkTag >( ptr );
    }

  inline
  ParallelScan( const FunctorType & arg_functor
              , const Policy      & arg_policy
              )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    {}
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Parallel patterns for Kokkos::Serial with TeamPolicy */

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Serial >
                 >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Serial > Policy ;
  typedef typename Policy::member_type                       Member ;

  const FunctorType  m_functor ;
  const int          m_league ;
  const int          m_shared ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec() const
    {
      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( Member(ileague,m_league,m_shared) );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec() const
    {
      const TagType t{} ;
      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( t , Member(ileague,m_league,m_shared) );
      }
    }

public:

  inline
  void execute() const
    {
      Kokkos::Serial::scratch_memory_resize( 0 , m_shared );
      this-> template exec< typename Policy::work_tag >();
    }

  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_league(  arg_policy.league_size() )
    , m_shared( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( arg_functor , 1 ) )
    { }
};

/*--------------------------------------------------------------------------*/

template< class FunctorType , class Arg0 , class Arg1 >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Serial >
                    >
{
private:

  typedef Kokkos::TeamPolicy< Arg0 , Arg1 , Kokkos::Serial > Policy ;
  typedef typename Policy::member_type                       Member ;
  typedef typename Policy::work_tag                          WorkTag ;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , WorkTag > ValueTraits ;
  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag > ValueInit ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const int          m_league ;
  const int          m_shared ;
        pointer_type m_result_ptr ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      reference_type update = ValueInit::init( m_functor , ptr );

      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( Member(ileague,m_league,m_shared) , update );
      }

      Kokkos::Impl::FunctorFinal< FunctorType , TagType >::
        final( m_functor , ptr );
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( pointer_type ptr ) const
    {
      const TagType t{} ;

      reference_type update = ValueInit::init( m_functor , ptr );

      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( t , Member(ileague,m_league,m_shared) , update );
      }

      Kokkos::Impl::FunctorFinal< FunctorType , TagType >::
        final( m_functor , ptr );
    }

public:

  inline
  void execute() const
    {
      pointer_type ptr = (pointer_type) Kokkos::Serial::scratch_memory_resize
           ( ValueTraits::value_size( m_functor ) , m_shared );

      this-> template exec< WorkTag >( m_result_ptr ? m_result_ptr : ptr );
    }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType     & arg_result
                )
    : m_functor( arg_functor )
    , m_league( arg_policy.league_size() )
    , m_shared( arg_policy.scratch_size() + FunctorTeamShmemSize< FunctorType >::value( m_functor , 1 ) )
    , m_result_ptr( arg_result.ptr_on_device() )
    {
      static_assert( Kokkos::is_view< ViewType >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View" );

      static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View in HostSpace" );
    }

};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Nested parallel patterns for Kokkos::Serial with TeamPolicy */

namespace Kokkos {
namespace Impl {

template<typename iType>
struct TeamThreadRangeBoundariesStruct<iType,SerialTeamMember> {
  typedef iType index_type;
  const iType begin ;
  const iType end ;
  enum {increment = 1};
  const SerialTeamMember& thread;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct (const SerialTeamMember& arg_thread, const iType& arg_count)
    : begin(0)
    , end(arg_count)
    , thread(arg_thread)
    {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct (const SerialTeamMember& arg_thread, const iType& arg_begin, const iType & arg_end )
    : begin( arg_begin )
    , end(   arg_end)
    , thread( arg_thread )
    {}
};

  template<typename iType>
  struct ThreadVectorRangeBoundariesStruct<iType,SerialTeamMember> {
    typedef iType index_type;
    enum {start = 0};
    const iType end;
    enum {increment = 1};

    KOKKOS_INLINE_FUNCTION
    ThreadVectorRangeBoundariesStruct (const SerialTeamMember& thread, const iType& count):
      end( count )
    {}
  };

} // namespace Impl

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>
TeamThreadRange( const Impl::SerialTeamMember& thread, const iType & count )
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>
TeamThreadRange( const Impl::SerialTeamMember& thread, const iType & begin , const iType & end )
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>(thread,begin,end);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >
  ThreadVectorRange(const Impl::SerialTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >(thread,count);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::SerialTeamMember> PerTeam(const Impl::SerialTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::SerialTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::SerialTeamMember> PerThread(const Impl::SerialTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::SerialTeamMember>(thread);
}

} // namespace Kokkos

namespace Kokkos {

  /** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
   *
   * The range i=0..N-1 is mapped to all threads of the the calling thread team.
   * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>& loop_boundaries, const Lambda& lambda) {
  for( iType i = loop_boundaries.begin; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>& loop_boundaries,
                     const Lambda & lambda, ValueType& result) {

  result = ValueType();

  for( iType i = loop_boundaries.begin; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }

  result = loop_boundaries.thread.team_reduce(result,Impl::JoinAdd<ValueType>());
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::SerialTeamMember>& loop_boundaries,
                     const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;

  for( iType i = loop_boundaries.begin; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }

  init_result = loop_boundaries.thread.team_reduce(result,Impl::JoinLambdaAdapter<ValueType,JoinType>(join));
}

} //namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >&
    loop_boundaries, const Lambda& lambda) {
  #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
  #pragma ivdep
  #endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a summation of
 * val is performed and put into result. This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >&
      loop_boundaries, const Lambda & lambda, ValueType& result) {
  result = ValueType();
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    result+=tmp;
  }
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i, ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread and a reduction of
 * val is performed using JoinType(ValueType& val, const ValueType& update) and put into init_result.
 * The input value of init_result is used as initializer for temporary variables of ValueType. Therefore
 * the input value should be the neutral element with respect to the join operation (e.g. '0 for +-' or
 * '1 for *'). This functionality requires C++11 support.*/
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >&
      loop_boundaries, const Lambda & lambda, const JoinType& join, ValueType& init_result) {

  ValueType result = init_result;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i,tmp);
    join(result,tmp);
  }
  init_result = result;
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes lambda(iType i, ValueType & val, bool final)
 *          for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan operation is performed.
 * Depending on the target execution space the operator might be called twice: once with final=false
 * and once with final=true. When final==true val contains the prefix sum value. The contribution of this
 * "i" needs to be added to val no matter whether final==true or not. In a serial execution
 * (i.e. team_size==1) the operator is only called once with final==true. Scan_val will be set
 * to the final sum value over all vector lanes.
 * This functionality requires C++11 support.*/
template< typename iType, class FunctorType >
KOKKOS_INLINE_FUNCTION
void parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::SerialTeamMember >&
      loop_boundaries, const FunctorType & lambda) {

  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void > ValueTraits ;
  typedef typename ValueTraits::value_type value_type ;

  value_type scan_val = value_type();

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,scan_val,true);
  }
}

} // namespace Kokkos

namespace Kokkos {

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::SerialTeamMember>& , const FunctorType& lambda) {
  lambda();
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::SerialTeamMember>& , const FunctorType& lambda) {
  lambda();
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::SerialTeamMember>& , const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::SerialTeamMember>& , const FunctorType& lambda, ValueType& val) {
  lambda(val);
}
}

#endif // defined( KOKKOS_HAVE_SERIAL )
#endif /* #define KOKKOS_SERIAL_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

