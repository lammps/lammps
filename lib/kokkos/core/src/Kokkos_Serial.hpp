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

/// \file Kokkos_Serial.hpp
/// \brief Declaration and definition of Kokkos::Serial device.

#ifndef KOKKOS_SERIAL_HPP
#define KOKKOS_SERIAL_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_SERIAL )

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

#include <Kokkos_UniqueToken.hpp>

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

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void fence() {}

  /** \brief  Return the maximum amount of concurrency.  */
  static int concurrency() {return 1;};

  //! Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool /* detail */ = false ) {}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  static bool sleep();
  static bool wake();

  static void initialize( unsigned threads_count = 1 ,
                          unsigned use_numa_count = 0 ,
                          unsigned use_cores_per_numa = 0 ,
                          bool allow_asynchronous_threadpool = false);

  static bool is_initialized();

  //! Free any resources being consumed by the device.
  static void finalize();

  //--------------------------------------------------------------------------

  inline static int thread_pool_size( int = 0 ) { return 1 ; }
  KOKKOS_INLINE_FUNCTION static int thread_pool_rank() { return 0 ; }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION static unsigned hardware_thread_id() { return thread_pool_rank(); }
  inline static unsigned max_hardware_threads() { return thread_pool_size(0); }
#else
  static void impl_initialize();

  static bool impl_is_initialized();

  //! Free any resources being consumed by the device.
  static void impl_finalize();

  //--------------------------------------------------------------------------

  inline static int impl_thread_pool_size( int = 0 ) { return 1 ; }
  KOKKOS_INLINE_FUNCTION static int impl_thread_pool_rank() { return 0 ; }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION static unsigned impl_hardware_thread_id() { return impl_thread_pool_rank(); }
  inline static unsigned impl_max_hardware_threads() { return impl_thread_pool_size(0); }
#endif

  static const char* name();
  //--------------------------------------------------------------------------
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<>
struct MemorySpaceAccess
  < Kokkos::Serial::memory_space
  , Kokkos::Serial::scratch_memory_space
  >
{
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy   = false };
};

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

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

// Resize thread team data scratch memory
void serial_resize_thread_team_data( size_t pool_reduce_bytes
                                   , size_t team_reduce_bytes
                                   , size_t team_shared_bytes
                                   , size_t thread_local_bytes );

HostThreadTeamData * serial_get_thread_team_data();

} /* namespace Impl */
} /* namespace Kokkos */


namespace Kokkos {
namespace Impl {

/*
 * < Kokkos::Serial , WorkArgTag >
 * < WorkArgTag , Impl::enable_if< std::is_same< Kokkos::Serial , Kokkos::DefaultExecutionSpace >::value >::type >
 *
 */
template< class ... Properties >
class TeamPolicyInternal< Kokkos::Serial , Properties ... >:public PolicyTraits<Properties...>
{
private:

  size_t m_team_scratch_size[2] ;
  size_t m_thread_scratch_size[2] ;
  int    m_league_size ;
  int    m_chunk_size;

public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicyInternal      execution_policy ;

  typedef PolicyTraits<Properties ... > traits;

  //! Execution space of this execution policy:
  typedef Kokkos::Serial  execution_space ;

  TeamPolicyInternal& operator = (const TeamPolicyInternal& p) {
    m_league_size = p.m_league_size;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    return *this;
  }

  //----------------------------------------
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  template< class FunctorType >
  static
  int team_size_max( const FunctorType & ) { return 1 ; }

  template< class FunctorType >
  static
  int team_size_recommended( const FunctorType & ) { return 1 ; }

  template< class FunctorType >
  static
  int team_size_recommended( const FunctorType & , const int& ) { return 1 ; }
#endif

  template<class FunctorType>
  int team_size_max( const FunctorType&, const ParallelForTag& ) const { return 1 ; }
  template<class FunctorType>
  int team_size_max( const FunctorType&, const ParallelReduceTag& ) const { return 1 ; }
  template<class FunctorType>
  int team_size_recommended( const FunctorType&, const ParallelForTag& ) const { return 1 ; }
  template<class FunctorType>
  int team_size_recommended( const FunctorType&, const ParallelReduceTag& ) const { return 1 ; }

  //----------------------------------------

  inline int team_size() const { return 1 ; }
  inline int league_size() const { return m_league_size ; }
  inline size_t scratch_size(const int& level, int = 0) const { return m_team_scratch_size[level] + m_thread_scratch_size[level]; }

  inline static
  int vector_length_max()
    { return 1024; } // Use arbitrary large number, is meant as a vectorizable length

  inline static
  int scratch_size_max(int level)
  { return (level==0?
        1024*32:
        20*1024*1024);
  }
  /** \brief  Specify league size, request team size */
  TeamPolicyInternal( execution_space &
            , int league_size_request
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE
            , int team_size_request
#else
            , int /* team_size_request */
#endif
            , int /* vector_length_request */ = 1 )
    : m_team_scratch_size { 0 , 0 }
    , m_thread_scratch_size { 0 , 0 }
    , m_league_size( league_size_request )
    , m_chunk_size ( 32 )
    {
      #ifndef KOKKOS_ENABLE_DEPRECATED_CODE
      if(team_size_request > 1) Kokkos::abort("Kokkos::abort: Requested Team Size is too large!");
      #endif
    }

  TeamPolicyInternal( execution_space &
            , int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_team_scratch_size { 0 , 0 }
    , m_thread_scratch_size { 0 , 0 }
    , m_league_size( league_size_request )
    , m_chunk_size ( 32 )
    {}

  TeamPolicyInternal( int league_size_request
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE
            , int team_size_request
#else
            , int /* team_size_request */
#endif
            , int /* vector_length_request */ = 1 )
    : m_team_scratch_size { 0 , 0 }
    , m_thread_scratch_size { 0 , 0 }
    , m_league_size( league_size_request )
    , m_chunk_size ( 32 )
    {
      #ifndef KOKKOS_ENABLE_DEPRECATED_CODE
      if(team_size_request > 1) Kokkos::abort("Kokkos::abort: Requested Team Size is too large!");
      #endif
    }

  TeamPolicyInternal( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
    : m_team_scratch_size { 0 , 0 }
    , m_thread_scratch_size { 0 , 0 }
    , m_league_size( league_size_request )
    , m_chunk_size ( 32 )
    {}

  inline int chunk_size() const { return m_chunk_size ; }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
    p.m_chunk_size = chunk_size_;
    return p;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    return p;
  }

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  }

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  }
#else
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal& set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) {
    m_team_scratch_size[level] = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
#endif

  typedef Impl::HostThreadTeamMember< Kokkos::Serial >  member_type ;

protected:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal internal_set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) {
    m_team_scratch_size[level] = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
#endif
};
} /* namespace Impl */
} /* namespace Kokkos */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Parallel patterns for Kokkos::Serial with RangePolicy */

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType ,
                   Kokkos::RangePolicy< Traits ... > ,
                   Kokkos::Serial
                 >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;

  const FunctorType m_functor ;
  const Policy      m_policy ;

  template< class TagType >
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec() const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i );
      }
    }

  template< class TagType >
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

template< class FunctorType , class ReducerType , class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::RangePolicy< Traits ... >
                    , ReducerType
                    , Kokkos::Serial
                    >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;
  typedef typename Policy::work_tag                                  WorkTag ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;

  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef FunctorAnalysis< FunctorPatternInterface::REDUCE , Policy , FunctorType > Analysis ;

  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd >  ValueInit ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i , update );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const TagType t{} ;

      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i , update );
      }
    }

public:

  inline
  void execute() const
    {
      const size_t pool_reduce_size =
        Analysis::value_size( ReducerConditional::select(m_functor , m_reducer) );
      const size_t team_reduce_size  = 0 ; // Never shrinks
      const size_t team_shared_size  = 0 ; // Never shrinks
      const size_t thread_local_size = 0 ; // Never shrinks

      serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

      HostThreadTeamData & data = *serial_get_thread_team_data();

      pointer_type ptr =
        m_result_ptr ? m_result_ptr : pointer_type(data.pool_reduce_local());

      reference_type update =
        ValueInit::init(  ReducerConditional::select(m_functor , m_reducer) , ptr );

      this-> template exec< WorkTag >( update );

      Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::
        final(  ReducerConditional::select(m_functor , m_reducer) , ptr );
    }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & arg_functor ,
                  const Policy       & arg_policy ,
                  const HostViewType & arg_result_view ,
                  typename std::enable_if<
                               Kokkos::is_view< HostViewType >::value &&
                              !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    : m_functor( arg_functor )
    , m_policy( arg_policy )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result_view.data() )
    {
      static_assert( Kokkos::is_view< HostViewType >::value
        , "Kokkos::Serial reduce result must be a View" );

      static_assert( std::is_same< typename HostViewType::memory_space , HostSpace >::value
        , "Kokkos::Serial reduce result must be a View in HostSpace" );
    }

  inline
  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ReducerType& reducer )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_reducer( reducer )
    , m_result_ptr(  reducer.view().data() )
    {
      /*static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
    }
};


/*--------------------------------------------------------------------------*/

template< class FunctorType , class ... Traits >
class ParallelScan< FunctorType
                  , Kokkos::RangePolicy< Traits ... >
                  , Kokkos::Serial
                  >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;
  typedef typename Policy::work_tag                                  WorkTag ;

  typedef FunctorAnalysis< FunctorPatternInterface::SCAN , Policy , FunctorType > Analysis ;

  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag >  ValueInit ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i , update , true );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const TagType t{} ;
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i , update , true );
      }
    }

public:

  inline
  void execute() const
    {
      const size_t pool_reduce_size = Analysis::value_size( m_functor );
      const size_t team_reduce_size  = 0 ; // Never shrinks
      const size_t team_shared_size  = 0 ; // Never shrinks
      const size_t thread_local_size = 0 ; // Never shrinks

      serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

      HostThreadTeamData & data = *serial_get_thread_team_data();

      reference_type update =
        ValueInit::init( m_functor , pointer_type(data.pool_reduce_local()) );

      this-> template exec< WorkTag >( update );
    }

  inline
  ParallelScan( const FunctorType & arg_functor
              , const Policy      & arg_policy
              )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    {}
};

/*--------------------------------------------------------------------------*/
template< class FunctorType , class ReturnType, class ... Traits >
class ParallelScanWithTotal< FunctorType
                           , Kokkos::RangePolicy< Traits ... >
                           , ReturnType
                           , Kokkos::Serial
                           >
{
private:

  typedef Kokkos::RangePolicy< Traits ... > Policy ;
  typedef typename Policy::work_tag                                  WorkTag ;

  typedef FunctorAnalysis< FunctorPatternInterface::SCAN , Policy , FunctorType > Analysis ;

  typedef Kokkos::Impl::FunctorValueInit<   FunctorType , WorkTag >  ValueInit ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType   m_functor ;
  const Policy        m_policy ;
  ReturnType & m_returnvalue;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( i , update , true );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( reference_type update ) const
    {
      const TagType t{} ;
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        m_functor( t , i , update , true );
      }
    }

public:

  inline
  void execute()
    {
      const size_t pool_reduce_size = Analysis::value_size( m_functor );
      const size_t team_reduce_size  = 0 ; // Never shrinks
      const size_t team_shared_size  = 0 ; // Never shrinks
      const size_t thread_local_size = 0 ; // Never shrinks

     serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

      HostThreadTeamData & data = *serial_get_thread_team_data();

      reference_type update =
        ValueInit::init( m_functor , pointer_type(data.pool_reduce_local()) );

      this-> template exec< WorkTag >( update );

      m_returnvalue = update;
    }

  inline
  ParallelScanWithTotal( const FunctorType & arg_functor
                       , const Policy      & arg_policy
                       , ReturnType        & arg_returnvalue
                       )
    : m_functor( arg_functor )
    , m_policy(  arg_policy )
    , m_returnvalue(  arg_returnvalue )
    {}
};

} // namespace Impl
} // namespace Kokkos


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Parallel patterns for Kokkos::Serial with MDRangePolicy */

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType ,
                   Kokkos::MDRangePolicy< Traits ... > ,
                   Kokkos::Serial
                 >
{
private:

  typedef Kokkos::MDRangePolicy< Traits ... > MDRangePolicy ;
  typedef typename MDRangePolicy::impl_range_policy Policy ;

  typedef typename Kokkos::Impl::HostIterateTile< MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void > iterate_type;

  const FunctorType   m_functor ;
  const MDRangePolicy m_mdr_policy ;
  const Policy        m_policy ;

  void
  exec() const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        iterate_type( m_mdr_policy, m_functor )( i );
      }
    }

public:

  inline
  void execute() const
    { this->exec(); }

  inline
  ParallelFor( const FunctorType   & arg_functor
             , const MDRangePolicy & arg_policy )
    : m_functor( arg_functor )
    , m_mdr_policy(  arg_policy )
    , m_policy( Policy(0, m_mdr_policy.m_num_tiles).set_chunk_size(1) )
    {}
};


template< class FunctorType , class ReducerType , class ... Traits >
class ParallelReduce< FunctorType
                    , Kokkos::MDRangePolicy< Traits ... >
                    , ReducerType
                    , Kokkos::Serial
                    >
{
private:

  typedef Kokkos::MDRangePolicy< Traits ... > MDRangePolicy ;
  typedef typename MDRangePolicy::impl_range_policy Policy ;

  typedef typename MDRangePolicy::work_tag                                  WorkTag ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef FunctorAnalysis< FunctorPatternInterface::REDUCE , MDRangePolicy , FunctorType > Analysis ;

  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd >  ValueInit ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::value_type      value_type ;
  typedef typename Analysis::reference_type  reference_type ;


  using iterate_type = typename Kokkos::Impl::HostIterateTile< MDRangePolicy
                                                             , FunctorType
                                                             , WorkTag
                                                             , reference_type
                                                             >;


  const FunctorType   m_functor ;
  const MDRangePolicy m_mdr_policy ;
  const Policy        m_policy ;
  const ReducerType   m_reducer ;
  const pointer_type  m_result_ptr ;

  inline
  void
  exec( reference_type update ) const
    {
      const typename Policy::member_type e = m_policy.end();
      for ( typename Policy::member_type i = m_policy.begin() ; i < e ; ++i ) {
        iterate_type( m_mdr_policy, m_functor, update )( i );
      }
    }

public:

  inline
  void execute() const
    {
      const size_t pool_reduce_size =
        Analysis::value_size( ReducerConditional::select(m_functor , m_reducer) );
      const size_t team_reduce_size  = 0 ; // Never shrinks
      const size_t team_shared_size  = 0 ; // Never shrinks
      const size_t thread_local_size = 0 ; // Never shrinks

      serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

      HostThreadTeamData & data = *serial_get_thread_team_data();

      pointer_type ptr =
        m_result_ptr ? m_result_ptr : pointer_type(data.pool_reduce_local());

      reference_type update =
        ValueInit::init(  ReducerConditional::select(m_functor , m_reducer) , ptr );

      this-> exec( update );

      Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::
        final(  ReducerConditional::select(m_functor , m_reducer) , ptr );
    }

  template< class HostViewType >
  ParallelReduce( const FunctorType  & arg_functor ,
                  const MDRangePolicy       & arg_policy ,
                  const HostViewType & arg_result_view ,
                  typename std::enable_if<
                               Kokkos::is_view< HostViewType >::value &&
                              !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    : m_functor( arg_functor )
    , m_mdr_policy( arg_policy )
    , m_policy( Policy(0, m_mdr_policy.m_num_tiles).set_chunk_size(1) )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result_view.data() )
    {
      static_assert( Kokkos::is_view< HostViewType >::value
        , "Kokkos::Serial reduce result must be a View" );

      static_assert( std::is_same< typename HostViewType::memory_space , HostSpace >::value
        , "Kokkos::Serial reduce result must be a View in HostSpace" );
    }

  inline
  ParallelReduce( const FunctorType & arg_functor
                , MDRangePolicy       arg_policy
                , const ReducerType& reducer )
    : m_functor( arg_functor )
    , m_mdr_policy(  arg_policy )
    , m_policy( Policy(0, m_mdr_policy.m_num_tiles).set_chunk_size(1) )
    , m_reducer( reducer )
    , m_result_ptr(  reducer.view().data() )
    {
      /*static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
    }
};



} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Parallel patterns for Kokkos::Serial with TeamPolicy */

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Properties >
class ParallelFor< FunctorType
                 , Kokkos::TeamPolicy< Properties ... >
                 , Kokkos::Serial
                 >
{
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef TeamPolicyInternal< Kokkos::Serial , Properties ...> Policy ;
  typedef typename Policy::member_type                       Member ;

  const FunctorType  m_functor ;
  const int          m_league ;
  const int          m_shared ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( HostThreadTeamData & data ) const
    {
      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( Member(data,ileague,m_league) );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( HostThreadTeamData & data ) const
    {
      const TagType t{} ;
      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( t , Member(data,ileague,m_league) );
      }
    }

public:

  inline
  void execute() const
    {
      const size_t pool_reduce_size  = 0 ; // Never shrinks
      const size_t team_reduce_size  = TEAM_REDUCE_SIZE ;
      const size_t team_shared_size  = m_shared ;
      const size_t thread_local_size = 0 ; // Never shrinks

      serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );

      HostThreadTeamData & data = *serial_get_thread_team_data();

      this->template exec< typename Policy::work_tag >( data );
    }

  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_functor( arg_functor )
    , m_league(  arg_policy.league_size() )
    , m_shared( arg_policy.scratch_size(0) +
                arg_policy.scratch_size(1) +
                FunctorTeamShmemSize< FunctorType >::value( arg_functor , 1 ) )
    { }
};

/*--------------------------------------------------------------------------*/

template< class FunctorType , class ReducerType , class ... Properties >
class ParallelReduce< FunctorType
                    , Kokkos::TeamPolicy< Properties ... >
                    , ReducerType
                    , Kokkos::Serial
                    >
{
private:

  enum { TEAM_REDUCE_SIZE = 512 };

  typedef TeamPolicyInternal< Kokkos::Serial, Properties ... > Policy ;

  typedef FunctorAnalysis< FunctorPatternInterface::REDUCE , Policy , FunctorType > Analysis ;

  typedef typename Policy::member_type                       Member ;
  typedef typename Policy::work_tag                          WorkTag ;

  typedef Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, FunctorType, ReducerType> ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef typename Kokkos::Impl::if_c< std::is_same<InvalidType,ReducerType>::value, WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueInit<   ReducerTypeFwd , WorkTagFwd >  ValueInit ;

  typedef typename Analysis::pointer_type    pointer_type ;
  typedef typename Analysis::reference_type  reference_type ;

  const FunctorType  m_functor ;
  const int          m_league ;
  const ReducerType  m_reducer ;
        pointer_type m_result_ptr ;
  const int          m_shared ;

  template< class TagType >
  inline
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec( HostThreadTeamData & data , reference_type update ) const
    {
      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( Member(data,ileague,m_league) , update );
      }
    }

  template< class TagType >
  inline
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec( HostThreadTeamData & data , reference_type update ) const
    {
      const TagType t{} ;

      for ( int ileague = 0 ; ileague < m_league ; ++ileague ) {
        m_functor( t , Member(data,ileague,m_league) , update );
      }
    }

public:

  inline
  void execute() const
    {
      const size_t pool_reduce_size  =
        Analysis::value_size( ReducerConditional::select(m_functor, m_reducer));

      const size_t team_reduce_size  = TEAM_REDUCE_SIZE ;
      const size_t team_shared_size  = m_shared ;
      const size_t thread_local_size = 0 ; // Never shrinks

      serial_resize_thread_team_data( pool_reduce_size
                                    , team_reduce_size
                                    , team_shared_size
                                    , thread_local_size );


      HostThreadTeamData & data = *serial_get_thread_team_data();

      pointer_type ptr =
        m_result_ptr ? m_result_ptr : pointer_type(data.pool_reduce_local());

      reference_type update =
        ValueInit::init( ReducerConditional::select(m_functor , m_reducer) , ptr );

      this-> template exec< WorkTag >( data , update );

      Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::
        final(  ReducerConditional::select(m_functor , m_reducer) , ptr );
    }

  template< class ViewType >
  ParallelReduce( const FunctorType  & arg_functor
                , const Policy       & arg_policy
                , const ViewType     & arg_result ,
                typename std::enable_if<
                  Kokkos::is_view< ViewType >::value &&
                  !Kokkos::is_reducer_type<ReducerType>::value
                  ,void*>::type = NULL)
    : m_functor( arg_functor )
    , m_league( arg_policy.league_size() )
    , m_reducer( InvalidType() )
    , m_result_ptr( arg_result.data() )
    , m_shared( arg_policy.scratch_size(0) +
                arg_policy.scratch_size(1) +
                FunctorTeamShmemSize< FunctorType >::value( m_functor , 1 ) )
    {
      static_assert( Kokkos::is_view< ViewType >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View" );

      static_assert( std::is_same< typename ViewType::memory_space
                                      , Kokkos::HostSpace >::value
        , "Reduction result on Kokkos::Serial must be a Kokkos::View in HostSpace" );
    }

  inline
  ParallelReduce( const FunctorType & arg_functor
                , Policy       arg_policy
                , const ReducerType& reducer )
    : m_functor( arg_functor )
    , m_league(  arg_policy.league_size() )
    , m_reducer( reducer )
    , m_result_ptr(  reducer.view().data() )
    , m_shared( arg_policy.scratch_size(0) +
                arg_policy.scratch_size(1) +
                FunctorTeamShmemSize< FunctorType >::value( arg_functor , 1 ) )
  {
  /*static_assert( std::is_same< typename ViewType::memory_space
                          , Kokkos::HostSpace >::value
  , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace" );*/
  }

};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos { namespace Experimental {

template<>
class UniqueToken< Serial, UniqueTokenScope::Instance>
{
public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken( execution_space const& = execution_space() ) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const  noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release( int ) const noexcept {}
};

template<>
class UniqueToken< Serial, UniqueTokenScope::Global>
{
public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken( execution_space const& = execution_space() ) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const  noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release( int ) const noexcept {}
};

}} // namespace Kokkos::Experimental

#include <impl/Kokkos_Serial_Task.hpp>

#endif // defined( KOKKOS_ENABLE_SERIAL )
#endif /* #define KOKKOS_SERIAL_HPP */

