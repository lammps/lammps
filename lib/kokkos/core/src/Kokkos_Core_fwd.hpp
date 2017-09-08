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

#ifndef KOKKOS_CORE_FWD_HPP
#define KOKKOS_CORE_FWD_HPP

//----------------------------------------------------------------------------
// Kokkos_Macros.hpp does introspection on configuration options
// and compiler environment then sets a collection of #define macros.

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Utilities.hpp>

#include <Kokkos_UniqueToken.hpp>
#include <Kokkos_MasterLock.hpp>

//----------------------------------------------------------------------------
// Have assumed a 64bit build (8byte pointers) throughout the code base.

static_assert( sizeof(void*) == 8
             , "Kokkos assumes 64-bit build; i.e., 8-byte pointers" );

//----------------------------------------------------------------------------

namespace Kokkos {

struct AUTO_t {
  KOKKOS_INLINE_FUNCTION
  constexpr const AUTO_t & operator()() const { return *this; }
};

namespace {
/**\brief Token to indicate that a parameter's value is to be automatically selected */
constexpr AUTO_t AUTO = Kokkos::AUTO_t();
}

struct InvalidType {};

} // namespace Kokkos

//----------------------------------------------------------------------------
// Forward declarations for class inter-relationships

namespace Kokkos {

class HostSpace; ///< Memory space for main process and CPU execution spaces

#ifdef KOKKOS_ENABLE_HBWSPACE
namespace Experimental {
class HBWSpace; /// Memory space for hbw_malloc from memkind (e.g. for KNL processor)
}
#endif

#if defined( KOKKOS_ENABLE_SERIAL )
class Serial;    ///< Execution space main process on CPU.
#endif

#if defined( KOKKOS_ENABLE_QTHREADS )
class Qthreads;  ///< Execution space with Qthreads back-end.
#endif

#if defined( KOKKOS_ENABLE_THREADS )
class Threads;   ///< Execution space with pthreads back-end.
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
class OpenMP;    ///< OpenMP execution space.
#endif

#if defined( KOKKOS_ENABLE_OPENMPTARGET )
namespace Experimental {
class OpenMPTarget;    ///< OpenMPTarget execution space.
class OpenMPTargetSpace;
}
#endif


#if defined( KOKKOS_ENABLE_CUDA )
class CudaSpace;            ///< Memory space on Cuda GPU
class CudaUVMSpace;         ///< Memory space on Cuda GPU with UVM
class CudaHostPinnedSpace;  ///< Memory space on Host accessible to Cuda GPU
class Cuda;                 ///< Execution space for Cuda GPU
#endif

#if defined( KOKKOS_ENABLE_ROCM )
namespace Experimental {
class ROCmSpace ;            ///< Memory space on ROCm GPU
class ROCm ;                 ///< Execution space for ROCm GPU
}
#endif

template<class ExecutionSpace, class MemorySpace>
struct Device;

} // namespace Kokkos

//----------------------------------------------------------------------------
// Set the default execution space.

/// Define Kokkos::DefaultExecutionSpace as per configuration option
/// or chosen from the enabled execution spaces in the following order:
/// Kokkos::Cuda, Kokkos::Experimental::OpenMPTarget, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Serial

namespace Kokkos {

#if   defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA )
  typedef Cuda DefaultExecutionSpace;
#elif defined ( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET )
  typedef Experimental::OpenMPTarget DefaultExecutionSpace ;
#elif defined ( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_ROCM )
  typedef Experimental::ROCm DefaultExecutionSpace ;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef OpenMP DefaultExecutionSpace;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Threads DefaultExecutionSpace;
//#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS )
//  typedef Qthreads DefaultExecutionSpace;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL )
  typedef Serial DefaultExecutionSpace;
#else
#  error "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::Cuda, Kokkos::Experimental::OpenMPTarget, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Qthreads, or Kokkos::Serial."
#endif

#if defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef OpenMP DefaultHostExecutionSpace;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Threads DefaultHostExecutionSpace;
//#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS )
//  typedef Qthreads DefaultHostExecutionSpace;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL )
  typedef Serial DefaultHostExecutionSpace;
#elif defined( KOKKOS_ENABLE_OPENMP )
  typedef OpenMP DefaultHostExecutionSpace;
#elif defined( KOKKOS_ENABLE_THREADS )
  typedef Threads DefaultHostExecutionSpace;
//#elif defined( KOKKOS_ENABLE_QTHREADS )
//  typedef Qthreads DefaultHostExecutionSpace;
#elif defined( KOKKOS_ENABLE_SERIAL )
  typedef Serial DefaultHostExecutionSpace;
#else
#  error "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::OpenMP, Kokkos::Threads, Kokkos::Qthreads, or Kokkos::Serial."
#endif

} // namespace Kokkos

//----------------------------------------------------------------------------
// Detect the active execution space and define its memory space.
// This is used to verify whether a running kernel can access
// a given memory space.

namespace Kokkos {

namespace Impl {

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA ) && defined( KOKKOS_ENABLE_CUDA )
typedef Kokkos::CudaSpace  ActiveExecutionMemorySpace;
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_ROCM_GPU )
typedef Kokkos::HostSpace  ActiveExecutionMemorySpace ;
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
typedef Kokkos::HostSpace  ActiveExecutionMemorySpace;
#else
typedef void ActiveExecutionMemorySpace;
#endif

template< class ActiveSpace, class MemorySpace >
struct VerifyExecutionCanAccessMemorySpace {
  enum {value = 0};
};

template< class Space >
struct VerifyExecutionCanAccessMemorySpace< Space, Space >
{
  enum {value = 1};
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void *) {}
};

} // namespace Impl

} // namespace Kokkos

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA( DATA_SPACE, DATA_PTR ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace, DATA_SPACE >::verify( DATA_PTR )

#define KOKKOS_RESTRICT_EXECUTION_TO_( DATA_SPACE ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace, DATA_SPACE >::verify()

//----------------------------------------------------------------------------

namespace Kokkos {
  void fence();
}

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template< class Functor
        , class Policy
        , class EnableFunctor = void
        , class EnablePolicy = void
        >
struct FunctorPolicyExecutionSpace;

//----------------------------------------------------------------------------
/// \class ParallelFor
/// \brief Implementation of the ParallelFor operator that has a
///   partial specialization for the device.
///
/// This is an implementation detail of parallel_for.  Users should
/// skip this and go directly to the nonmember function parallel_for.
template< class FunctorType, class ExecPolicy, class ExecutionSpace =
          typename Impl::FunctorPolicyExecutionSpace< FunctorType, ExecPolicy >::execution_space
        > class ParallelFor;

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType, class ExecPolicy, class ReducerType = InvalidType, class ExecutionSpace =
          typename Impl::FunctorPolicyExecutionSpace< FunctorType, ExecPolicy >::execution_space
        > class ParallelReduce;

/// \class ParallelScan
/// \brief Implementation detail of parallel_scan.
///
/// This is an implementation detail of parallel_scan.  Users should
/// skip this and go directly to the documentation of the nonmember
/// template function Kokkos::parallel_scan.
template< class FunctorType, class ExecPolicy, class ExecutionSapce =
          typename Impl::FunctorPolicyExecutionSpace< FunctorType, ExecPolicy >::execution_space
        > class ParallelScan;

} // namespace Impl

namespace Experimental {
template<class ScalarType , class Space = HostSpace> struct Sum;
template<class ScalarType , class Space = HostSpace> struct Prod;
template<class ScalarType , class Space = HostSpace> struct Min;
template<class ScalarType , class Space = HostSpace> struct Max;
template<class ScalarType , class Space = HostSpace> struct MinMax;
template<class ScalarType , class Index, class Space = HostSpace> struct MinLoc;
template<class ScalarType , class Index, class Space = HostSpace> struct MaxLoc;
template<class ScalarType , class Index, class Space = HostSpace> struct MinMaxLoc;
template<class ScalarType , class Space = HostSpace> struct BAnd;
template<class ScalarType , class Space = HostSpace> struct BOr;
template<class ScalarType , class Space = HostSpace> struct LAnd;
template<class ScalarType , class Space = HostSpace> struct LOr;
}
} // namespace Kokkos

#endif /* #ifndef KOKKOS_CORE_FWD_HPP */

