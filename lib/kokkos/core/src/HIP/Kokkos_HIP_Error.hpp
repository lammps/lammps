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

#ifndef KOKKOS_HIP_ERROR_HPP
#define KOKKOS_HIP_ERROR_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>

#include <hip/hip_runtime.h>

#include <ostream>

namespace Kokkos {
namespace Impl {

void hip_internal_error_throw(hipError_t e, const char* name,
                              const char* file = nullptr, const int line = 0);

inline void hip_internal_safe_call(hipError_t e, const char* name,
                                   const char* file = nullptr,
                                   const int line   = 0) {
  if (hipSuccess != e) {
    hip_internal_error_throw(e, name, file, line);
  }
}

}  // namespace Impl
}  // namespace Kokkos

#define KOKKOS_IMPL_HIP_SAFE_CALL(call) \
  Kokkos::Impl::hip_internal_safe_call(call, #call, __FILE__, __LINE__)

namespace Kokkos {
namespace Experimental {

class HIPRawMemoryAllocationFailure : public RawMemoryAllocationFailure {
 private:
  hipError_t m_error_code = hipSuccess;

  static FailureMode get_failure_mode(hipError_t error_code) {
    switch (error_code) {
      case hipErrorMemoryAllocation: return FailureMode::OutOfMemoryError;
      case hipErrorInvalidValue: return FailureMode::InvalidAllocationSize;
      default: return FailureMode::Unknown;
    }
  }

 public:
  HIPRawMemoryAllocationFailure(size_t arg_attempted_size,
                                hipError_t arg_error_code,
                                AllocationMechanism arg_mechanism) noexcept
      : RawMemoryAllocationFailure(
            arg_attempted_size, /* HIPSpace doesn't handle alignment? */ 1,
            get_failure_mode(arg_error_code), arg_mechanism),
        m_error_code(arg_error_code) {}

  void append_additional_error_information(std::ostream& o) const override {
    if (m_error_code != hipSuccess) {
      o << "  The HIP allocation returned the error code \""
        << hipGetErrorName(m_error_code) << "\".";
    }
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
