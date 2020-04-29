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

#ifndef KOKKOS_HIP_ERROR_HPP
#define KOKKOS_HIP_ERROR_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>

#include <hip/hip_runtime.h>

#include <ostream>

namespace Kokkos {
namespace Impl {

void hip_internal_error_throw(hipError_t e, const char* name,
                              const char* file = NULL, const int line = 0);

inline void hip_internal_safe_call(hipError_t e, const char* name,
                                   const char* file = NULL,
                                   const int line   = 0) {
  if (hipSuccess != e) {
    hip_internal_error_throw(e, name, file, line);
  }
}

}  // namespace Impl
}  // namespace Kokkos

#define HIP_SAFE_CALL(call) \
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
      o << "  The HIP allocation returned the error code \"\""
        << hipGetErrorName(m_error_code) << "\".";
    }
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
