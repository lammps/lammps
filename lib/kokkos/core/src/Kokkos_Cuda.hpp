/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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

#ifndef KOKKOS_CUDA_HPP
#define KOKKOS_CUDA_HPP

#include <Kokkos_Macros.hpp>

//----------------------------------------------------------------------------
// If CUDA execution space is enabled then use this header file.

#if defined( KOKKOS_HAVE_CUDA )

#if defined( __CUDACC__ )

#include <cuda.h>

/*  Compiling with a CUDA compiler.
 *
 *  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
 *    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
 *
 *  When generating device code the __CUDA_ARCH__ macro is defined as:
 *    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )
 */
#if ! defined( CUDA_VERSION )
#error "#include <cuda.h> did not define CUDA_VERSION"
#endif

#if ( CUDA_VERSION < 4010 )
#error "Cuda version 4.1 or greater required"
#endif

#if defined( __CUDA_ARCH__ ) && ( __CUDA_ARCH__ < 200 )
/*  Compiling with CUDA compiler for device code. */
#error "Cuda device capability >= 2.0 is required"
#endif

#endif /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#if defined( KOKKOS_HAVE_OPENMP )
#include <Kokkos_OpenMP.hpp>
#elif defined( KOKKOS_HAVE_PTHREAD )
#include <Kokkos_Threads.hpp>
#else
#endif

#include <Kokkos_Serial.hpp>

#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class CudaExec ;
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class Cuda
/// \brief Kokkos device that uses CUDA to run on GPUs.
///
/// A "device" represents a parallel execution model.  It tells Kokkos
/// how to parallelize the execution of kernels in a parallel_for or
/// parallel_reduce.  For example, the Threads device uses Pthreads or
/// C++11 threads on a CPU, the OpenMP device uses the OpenMP language
/// extensions, and the Serial device executes "parallel" kernels
/// sequentially.  The Cuda device uses NVIDIA's CUDA programming
/// model to execute kernels in parallel on GPUs.
class Cuda {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! The tag (what type of kokkos_object is this).
  typedef Impl::ExecutionSpaceTag  kokkos_tag ;
  //! The device type (same as this class).
  typedef Cuda                  device_type ;
  //! This device's execution space.
  typedef Cuda                  execution_space ;
  //! This device's preferred memory space.
  typedef CudaSpace             memory_space ;
  //! The size_type typedef best suited for this device.
  typedef CudaSpace::size_type  size_type ;
  //! This device's preferred array layout.
  typedef LayoutLeft            array_layout ;

  typedef ScratchMemorySpace< Cuda >  scratch_memory_space ;

  //! This device's host mirror type.
#if defined( KOKKOS_HAVE_OPENMP )
  typedef Kokkos::OpenMP       host_mirror_device_type ;
#elif defined( KOKKOS_HAVE_PTHREAD )
  typedef Kokkos::Threads      host_mirror_device_type ;
#else
  typedef Kokkos::Serial       host_mirror_device_type ;
#endif

  //@}
  //! \name Functions that all Kokkos devices must implement.
  //@{

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined( __CUDA_ARCH__ )
    return true;
#else
    return false;
#endif
  }

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
  static void fence();

  //! Free any resources being consumed by the device.
  static void finalize();

  //! Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool detail = false );

  //@}
  //--------------------------------------------------------------------------
  //! \name Device-specific functions
  //@{

  struct SelectDevice {
    int cuda_device_id ;
    SelectDevice() : cuda_device_id(0) {}
    explicit SelectDevice( int id ) : cuda_device_id( id ) {}
  };

  //! Initialize, telling the CUDA run-time library which device to use.
  static void initialize( const SelectDevice = SelectDevice() );
  static void initialize( int device );
  static void initialize( int device , int );

  static int is_initialized();

  /// \brief Cuda device architecture of the selected device.
  ///
  /// This matches the __CUDA_ARCH__ specification.
  static size_type device_arch();

  //! Query device count.
  static size_type detect_device_count();

  /** \brief  Detect the available devices and their architecture
   *          as defined by the __CUDA_ARCH__ specification.
   */
  static std::vector<unsigned> detect_device_arch();

  static unsigned team_max();
  static unsigned team_recommended();

  //@}
  //--------------------------------------------------------------------------
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<>
struct VerifyExecutionCanAccessMemorySpace
  < Kokkos::Cuda::memory_space
  , Kokkos::Cuda::scratch_memory_space
  >
{
  KOKKOS_INLINE_FUNCTION static void verify( void ) { }
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) { }
};

template<>
struct VerifyExecutionCanAccessMemorySpace
  < Kokkos::HostSpace
  , Kokkos::Cuda::scratch_memory_space
  >
{
  inline static void verify( void ) { CudaSpace::access_error(); }
  inline static void verify( const void * p ) { CudaSpace::access_error(p); }
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_View.hpp>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>

//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_HAVE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDA_HPP */



