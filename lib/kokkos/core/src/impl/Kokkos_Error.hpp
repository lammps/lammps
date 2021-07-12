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

#ifndef KOKKOS_ABORT_MESSAGE_BUFFER_SIZE
#define KOKKOS_ABORT_MESSAGE_BUFFER_SIZE 2048
#endif  // ifndef KOKKOS_ABORT_MESSAGE_BUFFER_SIZE

namespace Kokkos {
namespace Impl {

[[noreturn]] void host_abort(const char *const);

void throw_runtime_exception(const std::string &);

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
    SYCLMallocDevice,
    SYCLMallocShared
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

  KOKKOS_ATTRIBUTE_NODISCARD
  const char *what() const noexcept override {
    if (m_failure_mode == FailureMode::OutOfMemoryError) {
      return "Memory allocation error: out of memory";
    } else if (m_failure_mode == FailureMode::AllocationNotAligned) {
      return "Memory allocation error: allocation result was under-aligned";
    }

    return nullptr;  // unreachable
  }

  KOKKOS_ATTRIBUTE_NODISCARD
  size_t attempted_size() const noexcept { return m_attempted_size; }

  KOKKOS_ATTRIBUTE_NODISCARD
  size_t attempted_alignment() const noexcept { return m_attempted_alignment; }

  KOKKOS_ATTRIBUTE_NODISCARD
  AllocationMechanism allocation_mechanism() const noexcept {
    return m_mechanism;
  }

  KOKKOS_ATTRIBUTE_NODISCARD
  FailureMode failure_mode() const noexcept { return m_failure_mode; }

  void print_error_message(std::ostream &o) const;
  KOKKOS_ATTRIBUTE_NODISCARD
  std::string get_error_message() const;

  virtual void append_additional_error_information(std::ostream &) const {}
};

}  // end namespace Experimental

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA) && defined(__CUDA_ARCH__)

#if defined(__APPLE__) || defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
// cuda_abort does not abort when building for macOS.
// required to workaround failures in random number generator unit tests with
// pre-volta architectures
#define KOKKOS_IMPL_ABORT_NORETURN
#else
// cuda_abort aborts when building for other platforms than macOS
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#endif

#elif defined(KOKKOS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)
// HIP aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)
// FIXME_SYCL SYCL doesn't abort
#define KOKKOS_IMPL_ABORT_NORETURN
#elif !defined(KOKKOS_ENABLE_OPENMPTARGET)
// Host aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#else
// Everything else does not abort
#define KOKKOS_IMPL_ABORT_NORETURN
#endif

namespace Kokkos {
KOKKOS_IMPL_ABORT_NORETURN KOKKOS_INLINE_FUNCTION void abort(
    const char *const message) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(__CUDA_ARCH__)
  Kokkos::Impl::cuda_abort(message);
#elif defined(KOKKOS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)
  Kokkos::Impl::hip_abort(message);
#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)
  Kokkos::Impl::sycl_abort(message);
#elif !defined(KOKKOS_ENABLE_OPENMPTARGET)
  Kokkos::Impl::host_abort(message);
#else
  (void)message;  // FIXME_OPENMPTARGET
#endif
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if !defined(NDEBUG) || defined(KOKKOS_ENFORCE_CONTRACTS) || \
    defined(KOKKOS_ENABLE_DEBUG)
#define KOKKOS_EXPECTS(...)                                               \
  {                                                                       \
    if (!bool(__VA_ARGS__)) {                                             \
      ::Kokkos::abort(                                                    \
          "Kokkos contract violation:\n  "                                \
          "  Expected precondition `" #__VA_ARGS__ "` evaluated false."); \
    }                                                                     \
  }
#define KOKKOS_ENSURES(...)                                               \
  {                                                                       \
    if (!bool(__VA_ARGS__)) {                                             \
      ::Kokkos::abort(                                                    \
          "Kokkos contract violation:\n  "                                \
          "  Ensured postcondition `" #__VA_ARGS__ "` evaluated false."); \
    }                                                                     \
  }
// some projects already define this for themselves, so don't mess them up
#ifndef KOKKOS_ASSERT
#define KOKKOS_ASSERT(...)                                             \
  {                                                                    \
    if (!bool(__VA_ARGS__)) {                                          \
      ::Kokkos::abort(                                                 \
          "Kokkos contract violation:\n  "                             \
          "  Asserted condition `" #__VA_ARGS__ "` evaluated false."); \
    }                                                                  \
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
