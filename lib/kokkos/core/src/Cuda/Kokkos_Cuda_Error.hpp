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

#ifndef KOKKOS_CUDA_ERROR_HPP
#define KOKKOS_CUDA_ERROR_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <impl/Kokkos_Error.hpp>

#include <iostream>

namespace Kokkos {
namespace Impl {

void cuda_device_synchronize();

void cuda_internal_error_throw(cudaError e, const char* name,
                               const char* file = NULL, const int line = 0);

inline void cuda_internal_safe_call(cudaError e, const char* name,
                                    const char* file = NULL,
                                    const int line   = 0) {
  if (cudaSuccess != e) {
    cuda_internal_error_throw(e, name, file, line);
  }
}

#define CUDA_SAFE_CALL(call) \
  Kokkos::Impl::cuda_internal_safe_call(call, #call, __FILE__, __LINE__)

}  // namespace Impl

namespace Experimental {

class CudaRawMemoryAllocationFailure : public RawMemoryAllocationFailure {
 private:
  using base_t = RawMemoryAllocationFailure;

  cudaError_t m_error_code = cudaSuccess;

  static FailureMode get_failure_mode(cudaError_t error_code) {
    switch (error_code) {
      case cudaErrorMemoryAllocation: return FailureMode::OutOfMemoryError;
      case cudaErrorInvalidValue: return FailureMode::InvalidAllocationSize;
      // TODO handle cudaErrorNotSupported for cudaMallocManaged
      default: return FailureMode::Unknown;
    }
  }

 public:
  // using base_t::base_t;
  // would trigger
  //
  // error: cannot determine the exception specification of the default
  // constructor due to a circular dependency
  //
  // using NVCC 9.1 and gcc 7.4
  CudaRawMemoryAllocationFailure(
      size_t arg_attempted_size, size_t arg_attempted_alignment,
      FailureMode arg_failure_mode = FailureMode::OutOfMemoryError,
      AllocationMechanism arg_mechanism =
          AllocationMechanism::StdMalloc) noexcept
      : base_t(arg_attempted_size, arg_attempted_alignment, arg_failure_mode,
               arg_mechanism) {}

  CudaRawMemoryAllocationFailure(size_t arg_attempted_size,
                                 cudaError_t arg_error_code,
                                 AllocationMechanism arg_mechanism) noexcept
      : base_t(arg_attempted_size, /* CudaSpace doesn't handle alignment? */ 1,
               get_failure_mode(arg_error_code), arg_mechanism),
        m_error_code(arg_error_code) {}

  void append_additional_error_information(std::ostream& o) const override {
    if (m_error_code != cudaSuccess) {
      o << "  The Cuda allocation returned the error code \"\""
        << cudaGetErrorName(m_error_code) << "\".";
    }
  }
};

}  // end namespace Experimental

}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_CUDA
#endif  // KOKKOS_CUDA_ERROR_HPP
