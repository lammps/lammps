/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_CORE_FWD_HPP
#define KOKKOS_CORE_FWD_HPP

//----------------------------------------------------------------------------
// Kokkos_Macros.hpp does introspection on configuration options
// and compiler environment then sets a collection of #define macros.

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Utilities.hpp>

#include <Kokkos_MasterLock.hpp>

//----------------------------------------------------------------------------
// Have assumed a 64bit build (8byte pointers) throughout the code base.

static_assert(sizeof(void *) == 8,
              "Kokkos assumes 64-bit build; i.e., 8-byte pointers");

//----------------------------------------------------------------------------

namespace Kokkos {

struct AUTO_t {
  KOKKOS_INLINE_FUNCTION
  constexpr const AUTO_t &operator()() const { return *this; }
};

namespace {
/**\brief Token to indicate that a parameter's value is to be automatically
 * selected */
constexpr AUTO_t AUTO = Kokkos::AUTO_t();
}  // namespace

struct InvalidType {};

}  // namespace Kokkos

//----------------------------------------------------------------------------
// Forward declarations for class inter-relationships

namespace Kokkos {

class HostSpace;  ///< Memory space for main process and CPU execution spaces
class AnonymousSpace;

template <class ExecutionSpace, class MemorySpace>
struct Device;

// forward declare here so that backend initializer calls can use it.
struct InitArguments;

}  // namespace Kokkos

// Include backend forward statements as determined by build options
#include <KokkosCore_Config_FwdBackend.hpp>

//----------------------------------------------------------------------------
// Set the default execution space.

/// Define Kokkos::DefaultExecutionSpace as per configuration option
/// or chosen from the enabled execution spaces in the following order:
/// Kokkos::Cuda, Kokkos::Experimental::OpenMPTarget, Kokkos::OpenMP,
/// Kokkos::Threads, Kokkos::Serial

#if defined(__clang_analyzer__)
#define KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION \
  [[clang::annotate("DefaultExecutionSpace")]]
#define KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION \
  [[clang::annotate("DefaultHostExecutionSpace")]]
#else
#define KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION
#define KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION
#endif

namespace Kokkos {

#if defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION = Cuda;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Experimental::OpenMPTarget;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HIP)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Experimental::HIP;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SYCL)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Experimental::SYCL;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION = OpenMP;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION = Threads;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HPX)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Kokkos::Experimental::HPX;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION = Serial;
#else
#error \
    "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::Cuda, Kokkos::Experimental::HIP, Kokkos::Experimental::SYCL, Kokkos::Experimental::OpenMPTarget, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Experimental::HPX, or Kokkos::Serial."
#endif

#if defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    OpenMP;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Threads;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HPX)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Kokkos::Experimental::HPX;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Serial;
#elif defined(KOKKOS_ENABLE_OPENMP)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    OpenMP;
#elif defined(KOKKOS_ENABLE_THREADS)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Threads;
#elif defined(KOKKOS_ENABLE_HPX)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Kokkos::Experimental::HPX;
#elif defined(KOKKOS_ENABLE_SERIAL)
using DefaultHostExecutionSpace KOKKOS_IMPL_DEFAULT_HOST_EXEC_SPACE_ANNOTATION =
    Serial;
#else
#error \
    "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::OpenMP, Kokkos::Threads, Kokkos::Experimental::HPX, or Kokkos::Serial."
#endif

}  // namespace Kokkos

//----------------------------------------------------------------------------
// Detect the active execution space and define its memory space.
// This is used to verify whether a running kernel can access
// a given memory space.

namespace Kokkos {
namespace Impl {

#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA) && \
    defined(KOKKOS_ENABLE_CUDA)
using ActiveExecutionMemorySpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL)
using ActiveExecutionMemorySpace = Kokkos::Experimental::SYCLDeviceUSMSpace;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU)
using ActiveExecutionMemorySpace = Kokkos::Experimental::HIPSpace;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
using ActiveExecutionMemorySpace = Kokkos::HostSpace;
#else
using ActiveExecutionMemorySpace = void;
#endif

template <typename DstMemorySpace, typename SrcMemorySpace>
struct MemorySpaceAccess;

template <typename DstMemorySpace, typename SrcMemorySpace,
          bool = Kokkos::Impl::MemorySpaceAccess<DstMemorySpace,
                                                 SrcMemorySpace>::accessible>
struct verify_space {
  KOKKOS_FUNCTION static void check() {}
};

template <typename DstMemorySpace, typename SrcMemorySpace>
struct verify_space<DstMemorySpace, SrcMemorySpace, false> {
  KOKKOS_FUNCTION static void check() {
    Kokkos::abort(
        "Kokkos::View ERROR: attempt to access inaccessible memory space");
  };
};

// Base class for exec space initializer factories
class ExecSpaceInitializerBase;

}  // namespace Impl

namespace Experimental {
template <class, class, class, class>
class LogicalMemorySpace;
}

}  // namespace Kokkos

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA(DATA_SPACE, DATA_PTR)        \
  Kokkos::Impl::verify_space<Kokkos::Impl::ActiveExecutionMemorySpace, \
                             DATA_SPACE>::check();

#define KOKKOS_RESTRICT_EXECUTION_TO_(DATA_SPACE)                      \
  Kokkos::Impl::verify_space<Kokkos::Impl::ActiveExecutionMemorySpace, \
                             DATA_SPACE>::check();

//----------------------------------------------------------------------------

namespace Kokkos {
void fence();
}

//----------------------------------------------------------------------------

namespace Kokkos {

template <class DataType, class... Properties>
class View;

namespace Impl {

template <class DstSpace, class SrcSpace,
          class ExecutionSpace = typename DstSpace::execution_space>
struct DeepCopy;

template <class ViewType, class Layout = typename ViewType::array_layout,
          class ExecSpace = typename ViewType::execution_space,
          int Rank = ViewType::Rank, typename iType = int64_t>
struct ViewFill;

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          int Rank, typename iType>
struct ViewCopy;

template <class Functor, class Policy>
struct FunctorPolicyExecutionSpace;

//----------------------------------------------------------------------------
/// \class ParallelFor
/// \brief Implementation of the ParallelFor operator that has a
///   partial specialization for the device.
///
/// This is an implementation detail of parallel_for.  Users should
/// skip this and go directly to the nonmember function parallel_for.
template <class FunctorType, class ExecPolicy,
          class ExecutionSpace = typename Impl::FunctorPolicyExecutionSpace<
              FunctorType, ExecPolicy>::execution_space>
class ParallelFor;

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template <class FunctorType, class ExecPolicy, class ReducerType = InvalidType,
          class ExecutionSpace = typename Impl::FunctorPolicyExecutionSpace<
              FunctorType, ExecPolicy>::execution_space>
class ParallelReduce;

/// \class ParallelScan
/// \brief Implementation detail of parallel_scan.
///
/// This is an implementation detail of parallel_scan.  Users should
/// skip this and go directly to the documentation of the nonmember
/// template function Kokkos::parallel_scan.
template <class FunctorType, class ExecPolicy,
          class ExecutionSapce = typename Impl::FunctorPolicyExecutionSpace<
              FunctorType, ExecPolicy>::execution_space>
class ParallelScan;

template <class FunctorType, class ExecPolicy, class ReturnType = InvalidType,
          class ExecutionSapce = typename Impl::FunctorPolicyExecutionSpace<
              FunctorType, ExecPolicy>::execution_space>
class ParallelScanWithTotal;

}  // namespace Impl

template <class ScalarType, class Space = HostSpace>
struct Sum;
template <class ScalarType, class Space = HostSpace>
struct Prod;
template <class ScalarType, class Space = HostSpace>
struct Min;
template <class ScalarType, class Space = HostSpace>
struct Max;
template <class ScalarType, class Space = HostSpace>
struct MinMax;
template <class ScalarType, class Index, class Space = HostSpace>
struct MinLoc;
template <class ScalarType, class Index, class Space = HostSpace>
struct MaxLoc;
template <class ScalarType, class Index, class Space = HostSpace>
struct MinMaxLoc;
template <class ScalarType, class Space = HostSpace>
struct BAnd;
template <class ScalarType, class Space = HostSpace>
struct BOr;
template <class ScalarType, class Space = HostSpace>
struct LAnd;
template <class ScalarType, class Space = HostSpace>
struct LOr;

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_CORE_FWD_HPP */
