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

#ifndef KOKKOS_QTHREAD_HPP
#define KOKKOS_QTHREAD_HPP

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Layout.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class QthreadExec ;
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Execution space supported by Qthread */
class Qthread {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{
  //! The tag (what type of kokkos_object is this).
  typedef Impl::ExecutionSpaceTag  kokkos_tag ;

  typedef Qthread                  device_type ;
  typedef Qthread                  execution_space ;
  typedef Kokkos::HostSpace        memory_space ;
  typedef Qthread                  scratch_memory_space ;
  typedef memory_space::size_type  size_type ;
  typedef Kokkos::LayoutRight      array_layout ;
  typedef Kokkos::Qthread          host_mirror_device_type ;

  //@}
  /*------------------------------------------------------------------------*/
  /** \brief  Initialization will construct one or more instances */
  static Qthread & instance( int = 0 );

  /** \brief  Set the execution space to a "sleep" state.
   *
   * This function sets the "sleep" state in which it is not ready for work.
   * This may consume less resources than in an "ready" state,
   * but it may also take time to transition to the "ready" state.
   *
   * \return True if enters or is in the "sleep" state.
   *         False if functions are currently executing.
   */
  bool sleep();

  /** \brief  Wake from the sleep state.
   * 
   *  \return True if enters or is in the "ready" state.
   *          False if functions are currently executing.
   */
  static bool wake();

  /** \brief Wait until all dispatched functions to complete.
   * 
   *  The parallel_for or parallel_reduce dispatch of a functor may
   *  return asynchronously, before the functor completes.  This
   *  method does not return until all dispatched functors on this
   *  device have completed.
   */
  static void fence();

  /*------------------------------------------------------------------------*/

  Qthread( Impl::QthreadExec & e ) : m_exec(e) {}

  void * get_shmem( const int ) const ;

  static int team_recommended();
  static int team_max();

  /*------------------------------------------------------------------------*/

  static void initialize( int thread_count );
  static void finalize();

  /** \brief Print configuration information to the given output stream. */
  static void print_configuration( std::ostream & , const bool detail = false );

  int shepherd_size() const ;
  int shepherd_worker_size() const ;

private:

  friend class Impl::QthreadExec ;

  Impl::QthreadExec & m_exec ;

};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <Kokkos_Parallel.hpp>
#include <Qthread/Kokkos_QthreadExec.hpp>
#include <Qthread/Kokkos_Qthread_Parallel.hpp>

#endif /* #define KOKKOS_QTHREAD_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

