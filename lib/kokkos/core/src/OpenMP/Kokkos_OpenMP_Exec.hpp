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

#ifndef KOKKOS_OPENMPEXEC_HPP
#define KOKKOS_OPENMPEXEC_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_OPENMP )

#if !defined(_OPENMP) && !defined(__CUDA_ARCH__)
#error "You enabled Kokkos OpenMP support without enabling OpenMP in the compiler!"
#endif

#include <Kokkos_OpenMP.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>

#include <Kokkos_Atomic.hpp>

#include <Kokkos_UniqueToken.hpp>

#include <iostream>
#include <sstream>
#include <fstream>

#include <omp.h>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos { namespace Impl {

class OpenMPExec;

extern int g_openmp_hardware_max_threads;

extern __thread int t_openmp_hardware_id;
extern __thread OpenMPExec * t_openmp_instance;

//----------------------------------------------------------------------------
/** \brief  Data for OpenMP thread execution */

class OpenMPExec {
public:

  friend class Kokkos::OpenMP ;

  enum { MAX_THREAD_COUNT = 512 };

  void clear_thread_data();

  static void validate_partition( const int nthreads
                                , int & num_partitions
                                , int & partition_size
                                );

private:
  OpenMPExec( int arg_pool_size )
    : m_pool_size{ arg_pool_size }
    , m_level{ omp_get_level() }
    , m_pool()
  {}

  ~OpenMPExec()
  {
    clear_thread_data();
  }

  int m_pool_size;
  int m_level;

  HostThreadTeamData * m_pool[ MAX_THREAD_COUNT ];

public:

  static void verify_is_master( const char * const );

  void resize_thread_data( size_t pool_reduce_bytes
                         , size_t team_reduce_bytes
                         , size_t team_shared_bytes
                         , size_t thread_local_bytes );

  inline
  HostThreadTeamData * get_thread_data() const noexcept
  { return m_pool[ m_level == omp_get_level() ? 0 : omp_get_thread_num() ]; }

  inline
  HostThreadTeamData * get_thread_data( int i ) const noexcept
  { return m_pool[i]; }
};

}} // namespace Kokkos::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline OpenMP::OpenMP() noexcept
{}

inline
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
bool OpenMP::is_initialized() noexcept
#else
bool OpenMP::impl_is_initialized() noexcept
#endif
{ return Impl::t_openmp_instance != nullptr; }

inline
bool OpenMP::in_parallel( OpenMP const& ) noexcept
{
  //t_openmp_instance is only non-null on a master thread
  return   !Impl::t_openmp_instance
         || Impl::t_openmp_instance->m_level < omp_get_level()
         ;
}

inline
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
int OpenMP::thread_pool_size() noexcept
#else
int OpenMP::impl_thread_pool_size() noexcept
#endif
{
  return   OpenMP::in_parallel()
         ? omp_get_num_threads()
         : Impl::t_openmp_instance->m_pool_size
         ;
}

KOKKOS_INLINE_FUNCTION
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
int OpenMP::thread_pool_rank() noexcept
#else
int OpenMP::impl_thread_pool_rank() noexcept
#endif
{
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  return Impl::t_openmp_instance ? 0 : omp_get_thread_num();
#else
  return -1 ;
#endif
}

inline
void OpenMP::impl_static_fence( OpenMP const& instance ) noexcept {}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
inline
void OpenMP::fence( OpenMP const& instance ) noexcept {}
#endif

inline
bool OpenMP::is_asynchronous( OpenMP const& instance ) noexcept
{ return false; }

template <typename F>
void OpenMP::partition_master( F const& f
                             , int num_partitions
                             , int partition_size
                             )
{
  if (omp_get_nested()) {
    using Exec = Impl::OpenMPExec;

    Exec * prev_instance = Impl::t_openmp_instance;

    Exec::validate_partition( prev_instance->m_pool_size, num_partitions, partition_size );

    OpenMP::memory_space space;

    #pragma omp parallel num_threads(num_partitions)
    {
      void * const ptr = space.allocate( sizeof(Exec) );

      Impl::t_openmp_instance = new (ptr) Exec( partition_size );

      size_t pool_reduce_bytes  =   32 * partition_size ;
      size_t team_reduce_bytes  =   32 * partition_size ;
      size_t team_shared_bytes  = 1024 * partition_size ;
      size_t thread_local_bytes = 1024 ;

      Impl::t_openmp_instance->resize_thread_data( pool_reduce_bytes
                                                 , team_reduce_bytes
                                                 , team_shared_bytes
                                                 , thread_local_bytes
                                                 );

      omp_set_num_threads(partition_size);
      f( omp_get_thread_num(), omp_get_num_threads() );

      Impl::t_openmp_instance->~Exec();
      space.deallocate( Impl::t_openmp_instance, sizeof(Exec) );
      Impl::t_openmp_instance = nullptr;
    }

    Impl::t_openmp_instance  = prev_instance;
  }
  else {
    // nested openmp not enabled
    f(0,1);
  }
}


namespace Experimental {

template<>
class MasterLock<OpenMP>
{
public:
  void lock()     { omp_set_lock( &m_lock );   }
  void unlock()   { omp_unset_lock( &m_lock ); }
  bool try_lock() { return static_cast<bool>(omp_test_lock( &m_lock )); }

  MasterLock()  { omp_init_lock( &m_lock ); }
  ~MasterLock() { omp_destroy_lock( &m_lock ); }

  MasterLock( MasterLock const& ) = delete;
  MasterLock( MasterLock && )     = delete;
  MasterLock & operator=( MasterLock const& ) = delete;
  MasterLock & operator=( MasterLock && )     = delete;

private:
  omp_lock_t m_lock;

};

template<>
class UniqueToken< OpenMP, UniqueTokenScope::Instance>
{
public:
  using execution_space = OpenMP;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken( execution_space const& = execution_space() ) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( KOKKOS_ENABLE_DEPRECATED_CODE )
      return Kokkos::OpenMP::thread_pool_size();
#else
      return Kokkos::OpenMP::impl_thread_pool_size();
#endif
#else
      return 0 ;
#endif
    }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const  noexcept
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
#if defined( KOKKOS_ENABLE_DEPRECATED_CODE )
      return Kokkos::OpenMP::thread_pool_rank();
#else
      return Kokkos::OpenMP::impl_thread_pool_rank();
#endif
#else
      return 0 ;
#endif
    }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release( int ) const noexcept {}
};

template<>
class UniqueToken< OpenMP, UniqueTokenScope::Global>
{
public:
  using execution_space = OpenMP;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken( execution_space const& = execution_space() ) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept
    {
      #if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      return Kokkos::Impl::g_openmp_hardware_max_threads ;
      #else
      return 0 ;
      #endif
    }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept
    {
      #if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      return Kokkos::Impl::t_openmp_hardware_id ;
      #else
      return 0 ;
      #endif
    }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release( int ) const noexcept {}
};

} // namespace Experimental

inline
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
int OpenMP::thread_pool_size( int depth )
#else
int OpenMP::impl_thread_pool_size( int depth )
#endif
{
  return depth < 2
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
         ? thread_pool_size()
#else
         ? impl_thread_pool_size()
#endif
         : 1;
}

KOKKOS_INLINE_FUNCTION
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
int OpenMP::hardware_thread_id() noexcept
#else
int OpenMP::impl_hardware_thread_id() noexcept
#endif
{
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  return Impl::t_openmp_hardware_id;
#else
  return -1 ;
#endif
}

inline
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
int OpenMP::max_hardware_threads() noexcept
#else
int OpenMP::impl_max_hardware_threads() noexcept
#endif
{
  return Impl::g_openmp_hardware_max_threads;
}

} // namespace Kokkos

#endif
#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */

