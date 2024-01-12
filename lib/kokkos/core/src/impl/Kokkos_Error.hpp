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
#include <Kokkos_Abort.hpp>
#include <Kokkos_Assert.hpp>

namespace Kokkos {
namespace Impl {

[[noreturn]] void throw_runtime_exception(const std::string &msg);

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
    PosixMemAlign KOKKOS_DEPRECATED,
    PosixMMap KOKKOS_DEPRECATED,
    IntelMMAlloc KOKKOS_DEPRECATED,
#endif
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

#endif /* #ifndef KOKKOS_IMPL_ERROR_HPP */
