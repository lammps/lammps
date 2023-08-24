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

#ifndef KOKKOS_IMPL_ERROR_HPP
#define KOKKOS_IMPL_ERROR_HPP

#include <string>
#include <iosfwd>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_abort.hpp>
#endif
#ifdef KOKKOS_ENABLE_HIP
#include <HIP/Kokkos_HIP_Abort.hpp>
#endif
#ifdef KOKKOS_ENABLE_SYCL
#include <SYCL/Kokkos_SYCL_Abort.hpp>
#endif

namespace Kokkos {
namespace Impl {

[[noreturn]] void host_abort(const char *const);

#if defined(KOKKOS_ENABLE_CUDA) && defined(__CUDA_ARCH__)

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
// required to workaround failures in random number generator unit tests with
// pre-volta architectures
#define KOKKOS_IMPL_ABORT_NORETURN
#else
// cuda_abort aborts when building for other platforms than macOS
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#endif

#elif defined(KOKKOS_COMPILER_NVHPC)

#define KOKKOS_IMPL_ABORT_NORETURN

#elif defined(KOKKOS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)
// HIP aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)
// FIXME_SYCL SYCL doesn't abort
#define KOKKOS_IMPL_ABORT_NORETURN
#elif !defined(KOKKOS_ENABLE_OPENMPTARGET) && !defined(KOKKOS_ENABLE_OPENACC)
// Host aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#else
// Everything else does not abort
#define KOKKOS_IMPL_ABORT_NORETURN
#endif

// FIXME_SYCL
// Accomodate host pass for device functions that are not [[noreturn]]
#if defined(KOKKOS_ENABLE_SYCL) || \
    (defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK))
#define KOKKOS_IMPL_ABORT_NORETURN_DEVICE
#else
#define KOKKOS_IMPL_ABORT_NORETURN_DEVICE KOKKOS_IMPL_ABORT_NORETURN
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) ||          \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC)
KOKKOS_IMPL_ABORT_NORETURN_DEVICE inline KOKKOS_IMPL_DEVICE_FUNCTION void
device_abort(const char *const msg) {
#if defined(KOKKOS_ENABLE_CUDA)
  ::Kokkos::Impl::cuda_abort(msg);
#elif defined(KOKKOS_ENABLE_HIP)
  ::Kokkos::Impl::hip_abort(msg);
#elif defined(KOKKOS_ENABLE_SYCL)
  ::Kokkos::Impl::sycl_abort(msg);
#elif defined(KOKKOS_ENABLE_OPENMPTARGET) || defined(KOKKOS_ENABLE_OPENACC)
  printf("%s", msg);  // FIXME_OPENMPTARGET FIXME_OPENACC
#else
#error faulty logic
#endif
}
#endif

[[noreturn]] void throw_runtime_exception(const std::string &msg);

void traceback_callstack(std::ostream &);

std::string human_memory_size(size_t arg_bytes);

}  // namespace Impl

namespace Experimental {

class RawMemoryAllocationFailure : public std::bad_alloc {
 public:
  enum class FailureMode {
    OutOfMemoryError,
    AllocationNotAligned,
    InvalidAllocationSize,
    MaximumCudaUVMAllocationsExceeded,
    Unknown
  };
  enum class AllocationMechanism {
    StdMalloc,
    PosixMemAlign,
    PosixMMap,
    IntelMMAlloc,
    CudaMalloc,
    CudaMallocManaged,
    CudaHostAlloc,
    HIPMalloc,
    HIPHostMalloc,
    HIPMallocManaged,
    SYCLMallocDevice,
    SYCLMallocShared,
    SYCLMallocHost
  };

 private:
  size_t m_attempted_size;
  size_t m_attempted_alignment;
  FailureMode m_failure_mode;
  AllocationMechanism m_mechanism;

 public:
  RawMemoryAllocationFailure(
      size_t arg_attempted_size, size_t arg_attempted_alignment,
      FailureMode arg_failure_mode = FailureMode::OutOfMemoryError,
      AllocationMechanism arg_mechanism =
          AllocationMechanism::StdMalloc) noexcept
      : m_attempted_size(arg_attempted_size),
        m_attempted_alignment(arg_attempted_alignment),
        m_failure_mode(arg_failure_mode),
        m_mechanism(arg_mechanism) {}

  RawMemoryAllocationFailure() noexcept = delete;

  RawMemoryAllocationFailure(RawMemoryAllocationFailure const &) noexcept =
      default;
  RawMemoryAllocationFailure(RawMemoryAllocationFailure &&) noexcept = default;

  RawMemoryAllocationFailure &operator             =(
      RawMemoryAllocationFailure const &) noexcept = default;
  RawMemoryAllocationFailure &operator             =(
      RawMemoryAllocationFailure &&) noexcept = default;

  ~RawMemoryAllocationFailure() noexcept override = default;

  [[nodiscard]] const char *what() const noexcept override {
    if (m_failure_mode == FailureMode::OutOfMemoryError) {
      return "Memory allocation error: out of memory";
    } else if (m_failure_mode == FailureMode::AllocationNotAligned) {
      return "Memory allocation error: allocation result was under-aligned";
    }

    return nullptr;  // unreachable
  }

  [[nodiscard]] size_t attempted_size() const noexcept {
    return m_attempted_size;
  }

  [[nodiscard]] size_t attempted_alignment() const noexcept {
    return m_attempted_alignment;
  }

  [[nodiscard]] AllocationMechanism allocation_mechanism() const noexcept {
    return m_mechanism;
  }

  [[nodiscard]] FailureMode failure_mode() const noexcept {
    return m_failure_mode;
  }

  void print_error_message(std::ostream &o) const;
  [[nodiscard]] std::string get_error_message() const;

  virtual void append_additional_error_information(std::ostream &) const {}
};

}  // end namespace Experimental

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

KOKKOS_IMPL_ABORT_NORETURN KOKKOS_INLINE_FUNCTION void abort(
    const char *const message) {
  KOKKOS_IF_ON_HOST(::Kokkos::Impl::host_abort(message);)
  KOKKOS_IF_ON_DEVICE(::Kokkos::Impl::device_abort(message);)
}

#undef KOKKOS_IMPL_ABORT_NORETURN

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if !defined(NDEBUG) || defined(KOKKOS_ENFORCE_CONTRACTS) || \
    defined(KOKKOS_ENABLE_DEBUG)
#define KOKKOS_EXPECTS(...)                                                    \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Expected precondition `" #__VA_ARGS__                             \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
#define KOKKOS_ENSURES(...)                                                    \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Ensured postcondition `" #__VA_ARGS__                             \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
// some projects already define this for themselves, so don't mess
// them up
#ifndef KOKKOS_ASSERT
#define KOKKOS_ASSERT(...)                                                     \
  {                                                                            \
    if (!bool(__VA_ARGS__)) {                                                  \
      ::Kokkos::abort(                                                         \
          "Kokkos contract violation:\n  "                                     \
          "  Asserted condition `" #__VA_ARGS__                                \
          "` evaluated false.\n"                                               \
          "Error at " KOKKOS_IMPL_TOSTRING(__FILE__) ":" KOKKOS_IMPL_TOSTRING( \
              __LINE__) " \n");                                                \
    }                                                                          \
  }
#endif  // ifndef KOKKOS_ASSERT
#else   // not debug mode
#define KOKKOS_EXPECTS(...)
#define KOKKOS_ENSURES(...)
#ifndef KOKKOS_ASSERT
#define KOKKOS_ASSERT(...)
#endif  // ifndef KOKKOS_ASSERT
#endif  // end debug mode ifdefs

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_ERROR_HPP */
