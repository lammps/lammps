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

#include <iosfwd>
#include <vector>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#else
#ifdef KOKKOS_HAVE_PTHREAD
#include <Kokkos_Threads.hpp>
#else
#include <Kokkos_Serial.hpp>
#endif
#endif
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>

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

  //! The device type (same as this class).
  typedef Cuda                  device_type ;
  //! This device's preferred memory space.
  typedef CudaSpace             memory_space ;
  //! The size_type typedef best suited for this device.
  typedef CudaSpace::size_type  size_type ;
  //! This device's preferred array layout.
  typedef LayoutLeft            array_layout ;
  //! This device's host mirror type.
#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::OpenMP       host_mirror_device_type ;
#else
#ifdef KOKKOS_HAVE_PTHREAD
  typedef Kokkos::Threads       host_mirror_device_type ;
#else
  typedef Kokkos::Serial       host_mirror_device_type ;
#endif
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

  //@}
  //--------------------------------------------------------------------------
#if defined( __CUDA_ARCH__ )
  //! \name Functions for the functor device interface
  //@{

  __device__ inline int league_size() const { return gridDim.x ; }
  __device__ inline int league_rank() const { return blockIdx.x ; }

  __device__ inline int team_size() const { return blockDim.x ; }
  __device__ inline int team_rank() const { return threadIdx.x ; }

  __device__ inline void team_barrier() const { __syncthreads(); }
  __device__ inline unsigned int team_barrier_count(bool value) const
             { return __syncthreads_count(value); }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  __device__ inline Type team_scan( const Type & value );

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename TypeLocal , typename TypeGlobal >
  __device__ inline TypeGlobal team_scan( const TypeLocal & value , TypeGlobal * const global_accum );


  //! Get a pointer to shared memory for this team.
  __device__ inline void * get_shmem( const int size );

  __device__ inline Cuda( Impl::CudaExec & exec ) : m_exec(exec) {}
  __device__ inline Cuda( const Cuda & rhs ) : m_exec(rhs.m_exec) {}

  //@}
  //--------------------------------------------------------------------------

private:

  Impl::CudaExec & m_exec ;

  //--------------------------------------------------------------------------
#else

  int league_size() const ;
  int league_rank() const ;

  int team_size() const ;
  int team_rank() const ;

  void team_barrier() const ;
  unsigned int team_barrier_count(bool) const ;

  template< typename T >
    inline T team_scan(const T& value);

  template< typename TypeLocal , typename TypeGlobal >
    inline TypeGlobal team_scan( const TypeLocal & value , TypeGlobal * const global_accum );

  void * get_shmem( const int size );

  Cuda( Impl::CudaExec & );

#endif

};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief Cuda-specific parallel work configuration */

struct CudaWorkConfig {
  Cuda::size_type  grid[3] ;   //< Grid dimensions
  Cuda::size_type  block[3] ;  //< Block dimensions
  Cuda::size_type  shared ;    //< Shared memory size

  CudaWorkConfig()
  {
    enum { WarpSize = 32 };
    grid[0] = grid[1] = grid[2] = 1 ;
    block[1] = block[2] = 1 ;
    block[0] = 8 * WarpSize ;
    shared = 0 ;
  }
};

template< class FunctorType >
inline
void parallel_for( const CudaWorkConfig & work_config ,
                   const FunctorType    & functor )
{
  Impl::ParallelFor< FunctorType , CudaWorkConfig , Cuda >
    ( work_config , functor );
}

template< class FunctorType , class FinalizeType >
inline
void parallel_reduce( const CudaWorkConfig & work_config ,
                      const FunctorType    & functor ,
                      const FinalizeType   & finalize );

template< class FunctorType >
inline
typename FunctorType::value_type
parallel_reduce( const CudaWorkConfig & work_config ,
                 const FunctorType    & functor );

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_View.hpp>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>

#endif /* #ifndef KOKKOS_CUDA_HPP */

//----------------------------------------------------------------------------


