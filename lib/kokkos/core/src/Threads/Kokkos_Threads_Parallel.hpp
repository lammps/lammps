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

#ifndef KOKKOS_THREADS_PARALLEL_HPP
#define KOKKOS_THREADS_PARALLEL_HPP

#include <vector>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_ParallelReduce.hpp>

#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelFor< FunctorType , WorkSpec , Kokkos::Threads >
{
public:

  const FunctorType  m_func ;
  const size_t       m_work ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    const std::pair<size_t,size_t> work = exec.work_range( self.m_work );

    for ( size_t iwork = work.first, work_end = work.second ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork );
    }

    exec.fan_in();
  }

  ParallelFor( const FunctorType & functor , const size_t work )
    : m_func( functor ), m_work( work )
    {
      ThreadsExec::start( & ParallelFor::execute , this );
      ThreadsExec::fence();
    }

  inline void wait() {}

  inline ~ParallelFor() { wait(); }
};

template< class FunctorType >
class ParallelFor< FunctorType , ParallelWorkRequest , Kokkos::Threads >
{
public:

  const FunctorType  m_func ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelFor & self = * ((const ParallelFor *) arg );

    for ( ; exec.team_work_avail() ; exec.team_work_next() ) {
      self.m_func( Threads( exec ) );
    }

    exec.fan_in();
  }

  ParallelFor( const FunctorType & functor , const ParallelWorkRequest & work )
    : m_func( functor )
    {
      ThreadsExec::resize_shared_scratch( FunctorShmemSize< FunctorType >::value( functor ) );
      ThreadsExec::start( & ParallelFor::execute , this , work.league_size , work.team_size );
      ThreadsExec::fence();
    }

  inline void wait() {}

  inline ~ParallelFor() { wait(); }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class WorkSpec >
class ParallelReduce< FunctorType , WorkSpec , Kokkos::Threads >
{
public:

  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  const FunctorType  m_func ;
  const size_t       m_work ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    typename Reduce::reference_type update = Reduce::reference( exec.reduce_base() );

    self.m_func.init( update ); // Initialize thread-local value

    const std::pair<size_t,size_t> work = exec.work_range( self.m_work );

    for ( size_t iwork = work.first, work_end = work.second ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update );
    }

    exec.fan_in_reduce( self.m_func );
  }

  ParallelReduce( const FunctorType & functor ,
                  const size_t        work ,
                  const pointer_type  result_ptr = 0 )
    : m_func( functor ), m_work( work )
    {
      ThreadsExec::resize_reduce_scratch( Reduce::value_size( m_func ) );

      ThreadsExec::start( & ParallelReduce::execute , this );

      const pointer_type data = (pointer_type) ThreadsExec::root_reduce_scratch();

      ThreadsExec::fence();

      Reduce::final( m_func , data );

      if ( result_ptr ) {
        const unsigned n = Reduce::value_count( m_func );
        for ( unsigned i = 0 ; i < n ; ++i ) { result_ptr[i] = data[i]; }
      }
    }

  inline void wait() {}

  inline ~ParallelReduce() { wait(); }
};

template< class FunctorType >
class ParallelReduce< FunctorType , ParallelWorkRequest , Kokkos::Threads >
{
public:

  typedef ReduceAdapter< FunctorType >   Reduce ;
  typedef typename Reduce::pointer_type  pointer_type ;

  const FunctorType  m_func ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelReduce & self = * ((const ParallelReduce *) arg );

    typename Reduce::reference_type update = Reduce::reference( exec.reduce_base() );

    self.m_func.init( update ); // Initialize thread-local value

    for ( ; exec.team_work_avail() ; exec.team_work_next() ) {
      self.m_func( Threads( exec ) , update );
    }

    exec.fan_in_reduce( self.m_func );
  }

  ParallelReduce( const FunctorType & functor ,
                  const ParallelWorkRequest & work ,
                  const pointer_type  result_ptr = 0 )
    : m_func( functor )
    {
      ThreadsExec::resize_shared_scratch( FunctorShmemSize< FunctorType >::value( functor ) );
      ThreadsExec::resize_reduce_scratch( Reduce::value_size( m_func ) );

      ThreadsExec::start( & ParallelReduce::execute , this , work.league_size , work.team_size );

      const pointer_type data = (pointer_type) ThreadsExec::root_reduce_scratch();

      ThreadsExec::fence();

      Reduce::final( m_func , data );

      if ( result_ptr ) {
        const unsigned n = Reduce::value_count( m_func );
        for ( unsigned i = 0 ; i < n ; ++i ) { result_ptr[i] = data[i]; }
      }
    }

  inline void wait() {}

  inline ~ParallelReduce() { wait(); }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

struct ThreadsExecUseScanSmall {
  size_t nwork ;
  operator size_t () const { return nwork ; }
  ThreadsExecUseScanSmall( size_t n ) : nwork( n ) {}
};

template< class FunctorType , class WorkSpec >
class ParallelScan< FunctorType , WorkSpec , Kokkos::Threads >
{
public:

  typedef ReduceAdapter< FunctorType > Reduce ;
  typedef typename Reduce::pointer_type pointer_type ;

  const FunctorType  m_func ;
  const size_t       m_work ;

  static void execute( ThreadsExec & exec , const void * arg )
  {
    const ParallelScan & self = * ((const ParallelScan *) arg );

    const std::pair<size_t,size_t> work = exec.work_range( self.m_work );

    typename Reduce::reference_type update = Reduce::reference( exec.reduce_base() );

    self.m_func.init( update );

    for ( size_t iwork = work.first, work_end = work.second ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update , false );
    }

    // Compile time selection of scan algorithm to support unit testing
    // of both large and small thread count algorithms.
    if ( ! is_same< WorkSpec , ThreadsExecUseScanSmall >::value ) {
      exec.scan_large( self.m_func );
    }
    else {
      exec.scan_small( self.m_func );
    }

    for ( size_t iwork = work.first, work_end = work.second ; iwork < work_end ; ++iwork ) {
      self.m_func( iwork , update , true );
    }

    exec.fan_in();
  }

  ParallelScan( const FunctorType & functor , const size_t nwork )
    : m_func( functor )
    , m_work( nwork )
    {
      ThreadsExec::resize_reduce_scratch( 2 * Reduce::value_size( m_func ) );
      ThreadsExec::start( & ParallelScan::execute , this );
      ThreadsExec::fence();
    }

  inline void wait() {}

  inline ~ParallelScan() { wait(); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template<>
class MultiFunctorParallelReduce< Threads > {
private:

  struct MemberBase {
    virtual void init( Impl::ThreadsExec & ) const = 0 ;
    virtual void exec( Impl::ThreadsExec & ) const = 0 ;
    virtual void fan_in_reduce( Impl::ThreadsExec & ) const = 0 ;
    virtual void output( void * ) const = 0 ;
    virtual ~MemberBase() {}
  };

  template< class FunctorType >
  struct Member : public MemberBase {
    typedef Impl::ReduceAdapter< FunctorType >   Reduce ;
    typedef typename Reduce::pointer_type  pointer_type ;

    const FunctorType  m_func ;
    const size_t       m_work ;

    ~Member() {}

    Member( const FunctorType & func , const size_t work )
      : m_func( func ), m_work( work )
      {
        Impl::ThreadsExec::resize_reduce_scratch( Reduce::value_size( m_func ) );
      }

    void init( Impl::ThreadsExec & exec_ ) const
      { m_func.init( Reduce::reference( exec_.reduce_base() ) ); }

    void exec( Impl::ThreadsExec & exec_ ) const
      {
        typename Reduce::reference_type update = Reduce::reference( exec_.reduce_base() );

        const std::pair<size_t,size_t> work = exec_.work_range( m_work );

        for ( size_t iwork = work.first, work_end = work.second ; iwork < work_end ; ++iwork ) {
          m_func( iwork , update );
        }
      }

    void fan_in_reduce( Impl::ThreadsExec & exec_ ) const
      { exec_.fan_in_reduce( m_func ); }

    void output( void * ptr ) const
      {
        const pointer_type result = (pointer_type) ptr ;
        const pointer_type data   = (pointer_type) Impl::ThreadsExec::root_reduce_scratch();

        Impl::ThreadsExec::fence();

        Reduce::final( m_func , data );

        if ( result ) {
          const unsigned n = Reduce::value_count( m_func );
          for ( unsigned i = 0 ; i < n ; ++i ) { result[i] = data[i]; }
        }
      }
  };

  std::vector< MemberBase * > m_members ;

  static void execute_members( Impl::ThreadsExec & exec_ , const void * arg )
  {
    const MultiFunctorParallelReduce & self = * ((const MultiFunctorParallelReduce *) arg );

    // First functor initializes:

    self.m_members.front()->init( exec_ ); // Initialize thread-local value

    for ( unsigned i = 0 ; i < self.m_members.size() ; ++i ) {
      self.m_members[i]->exec( exec_ );
    }

    // Last functor fan-in reduce:

    self.m_members.back()->fan_in_reduce( exec_ );
  }

public:

  inline
  void execute( void * host_ptr ) const
    {
      if ( ! m_members.empty() ) {
        Impl::ThreadsExec::start( & MultiFunctorParallelReduce::execute_members , this );
        m_members.back()->output( host_ptr );
      }
    }

  inline
  void wait() const {}

  template< class FunctorType >
  void push_back( const size_t work_count , const FunctorType & f )
  {
    MemberBase * const m = new Member< FunctorType >( f , work_count );
    m_members.push_back( m );
  }

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_members.empty() ) {
      delete m_members.back();
      m_members.pop_back();
    }
  }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_THREADS_PARALLEL_HPP */

