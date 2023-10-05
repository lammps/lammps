//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_CORE_FWD_HPP
#define KOKKOS_CORE_FWD_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE_FWD
#endif

//----------------------------------------------------------------------------
// Kokkos_Macros.hpp does introspection on configuration options
// and compiler environment then sets a collection of #define macros.

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Utilities.hpp>

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
#include <Kokkos_MasterLock.hpp>
#endif

//----------------------------------------------------------------------------
// Have assumed a 64-bit build (8-byte pointers) throughout the code base.
// 32-bit build allowed but unsupported.
#ifdef KOKKOS_IMPL_32BIT
static_assert(sizeof(void *) == 4,
              "Kokkos assumes 64-bit build; i.e., 4-byte pointers");
#else
static_assert(sizeof(void *) == 8,
              "Kokkos assumes 64-bit build; i.e., 8-byte pointers");
#endif
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
// Forward declarations for class interrelationships

namespace Kokkos {

class HostSpace;  ///< Memory space for main process and CPU execution spaces
class AnonymousSpace;

template <class ExecutionSpace, class MemorySpace>
struct Device;

// forward declare here so that backend initializer calls can use it.
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
struct InitArguments;
#endif
class InitializationSettings;

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
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION = HIP;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SYCL)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Experimental::SYCL;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENACC)
using DefaultExecutionSpace KOKKOS_IMPL_DEFAULT_EXEC_SPACE_ANNOTATION =
    Experimental::OpenACC;
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
    "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::Cuda, Kokkos::HIP, Kokkos::Experimental::SYCL, Kokkos::Experimental::OpenMPTarget, Kokkos::Experimental::OpenACC, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Experimental::HPX, or Kokkos::Serial."
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

// check for devices that support sharedSpace
#if defined(KOKKOS_ENABLE_CUDA)
using SharedSpace = CudaUVMSpace;
#define KOKKOS_HAS_SHARED_SPACE
#elif defined(KOKKOS_ENABLE_HIP)
using SharedSpace = HIPManagedSpace;
#define KOKKOS_HAS_SHARED_SPACE
#elif defined(KOKKOS_ENABLE_SYCL)
using SharedSpace = Experimental::SYCLSharedUSMSpace;
#define KOKKOS_HAS_SHARED_SPACE
// if only host compile point to HostSpace
#elif !defined(KOKKOS_ENABLE_OPENACC) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
using SharedSpace               = HostSpace;
#define KOKKOS_HAS_SHARED_SPACE
#endif

inline constexpr bool has_shared_space =
#if defined KOKKOS_HAS_SHARED_SPACE
    true;
#else
    false;
#endif

#if defined(KOKKOS_ENABLE_CUDA)
using SharedHostPinnedSpace = CudaHostPinnedSpace;
#define KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
#elif defined(KOKKOS_ENABLE_HIP)
using SharedHostPinnedSpace = HIPHostPinnedSpace;
#define KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
#elif defined(KOKKOS_ENABLE_SYCL)
    using SharedHostPinnedSpace = Experimental::SYCLHostUSMSpace;
#define KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
#elif !defined(KOKKOS_ENABLE_OPENACC) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
    using SharedHostPinnedSpace = HostSpace;
#define KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
#endif

inline constexpr bool has_shared_host_pinned_space =
#if defined KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
    true;
#else
    false;
#endif

}  // namespace Kokkos

//----------------------------------------------------------------------------
// Detect the active execution space and define its memory space.
// This is used to verify whether a running kernel can access
// a given memory space.

namespace Kokkos {

template <class AccessSpace, class MemorySpace>
struct SpaceAccessibility;

namespace Impl {

// primary template: memory space is accessible, do nothing.
template <class MemorySpace, class AccessSpace,
          bool = SpaceAccessibility<AccessSpace, MemorySpace>::accessible>
struct RuntimeCheckMemoryAccessViolation {
  KOKKOS_FUNCTION RuntimeCheckMemoryAccessViolation(char const *const) {}
};

// explicit specialization: memory access violation will occur, call abort with
// the specified error message.
template <class MemorySpace, class AccessSpace>
struct RuntimeCheckMemoryAccessViolation<MemorySpace, AccessSpace, false> {
  KOKKOS_FUNCTION RuntimeCheckMemoryAccessViolation(char const *const msg) {
    Kokkos::abort(msg);
  }
};

// calls abort with default error message at runtime if memory access violation
// will occur
template <class MemorySpace>
KOKKOS_FUNCTION void runtime_check_memory_access_violation() {
  KOKKOS_IF_ON_HOST((
      RuntimeCheckMemoryAccessViolation<MemorySpace, DefaultHostExecutionSpace>(
          "ERROR: attempt to access inaccessible memory space");))
  KOKKOS_IF_ON_DEVICE(
      (RuntimeCheckMemoryAccessViolation<MemorySpace, DefaultExecutionSpace>(
           "ERROR: attempt to access inaccessible memory space");))
}

// calls abort with specified error message at runtime if memory access
// violation will occur
template <class MemorySpace>
KOKKOS_FUNCTION void runtime_check_memory_access_violation(
    char const *const msg) {
  KOKKOS_IF_ON_HOST((
      (void)RuntimeCheckMemoryAccessViolation<MemorySpace,
                                              DefaultHostExecutionSpace>(msg);))
  KOKKOS_IF_ON_DEVICE((
      (void)
          RuntimeCheckMemoryAccessViolation<MemorySpace, DefaultExecutionSpace>(
              msg);))
}

}  // namespace Impl

namespace Experimental {
template <class, class, class, class>
class LogicalMemorySpace;
}

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
// Getting ICE in Trilinos in Sacado and Intrepid in deep_copy
// See issue https://github.com/kokkos/kokkos/issues/5290
// Simply taking string by value did not resolve the issue
#ifdef KOKKOS_COMPILER_INTEL
void fence();
void fence(const std::string &name);
#else
void fence(const std::string &name = "Kokkos::fence: Unnamed Global Fence");
#endif
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

template <class DataType, class... Properties>
class View;

namespace Impl {

template <class DstSpace, class SrcSpace,
          class ExecutionSpace = typename DstSpace::execution_space,
          class Enable         = void>
struct DeepCopy;

template <class ViewType, class Layout = typename ViewType::array_layout,
          class ExecSpace = typename ViewType::execution_space,
          int Rank = ViewType::rank, typename iType = int64_t>
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
template <typename CombinedFunctorReducerType, typename PolicyType,
          typename ExecutionSpaceType>
class ParallelReduce;

template <typename FunctorType, typename FunctorAnalysisReducerType,
          typename Enable = void>
class CombinedFunctorReducer;

/// \class ParallelScan
/// \brief Implementation detail of parallel_scan.
///
/// This is an implementation detail of parallel_scan.  Users should
/// skip this and go directly to the documentation of the nonmember
/// template function Kokkos::parallel_scan.
template <class FunctorType, class ExecPolicy,
          class ExecutionSpace = typename Impl::FunctorPolicyExecutionSpace<
              FunctorType, ExecPolicy>::execution_space>
class ParallelScan;

template <class FunctorType, class ExecPolicy, class ReturnType = InvalidType,
          class ExecutionSpace = typename Impl::FunctorPolicyExecutionSpace<
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

template <class Scalar, class Index, class Space = HostSpace>
struct MaxFirstLoc;
template <class Scalar, class Index, class ComparatorType,
          class Space = HostSpace>
struct MaxFirstLocCustomComparator;

template <class Scalar, class Index, class Space = HostSpace>
struct MinFirstLoc;
template <class Scalar, class Index, class ComparatorType,
          class Space = HostSpace>
struct MinFirstLocCustomComparator;

template <class Scalar, class Index, class Space = HostSpace>
struct MinMaxFirstLastLoc;
template <class Scalar, class Index, class ComparatorType,
          class Space = HostSpace>
struct MinMaxFirstLastLocCustomComparator;

template <class Index, class Space = HostSpace>
struct FirstLoc;
template <class Index, class Space = HostSpace>
struct LastLoc;
template <class Index, class Space = HostSpace>
struct StdIsPartitioned;
template <class Index, class Space = HostSpace>
struct StdPartitionPoint;
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE_FWD
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE_FWD
#endif
#endif /* #ifndef KOKKOS_CORE_FWD_HPP */
