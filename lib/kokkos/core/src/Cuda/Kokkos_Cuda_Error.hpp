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

#ifndef KOKKOS_CUDA_ERROR_HPP
#define KOKKOS_CUDA_ERROR_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Profiling.hpp>

namespace Kokkos {
namespace Impl {

void cuda_device_synchronize(const std::string& name);
void cuda_stream_synchronize(const cudaStream_t stream,
                             const std::string& name);

[[noreturn]] void cuda_internal_error_throw(cudaError e, const char* name,
                                            const char* file = nullptr,
                                            const int line   = 0);

#ifndef KOKKOS_COMPILER_NVHPC
[[noreturn]]
#endif
             void cuda_internal_error_abort(cudaError e, const char* name,
                                            const char* file = nullptr,
                                            const int line   = 0);

inline void cuda_internal_safe_call(cudaError e, const char* name,
                                    const char* file = nullptr,
                                    const int line   = 0) {
  // 1. Success -> normal continuation.
  // 2. Error codes for which, to continue using CUDA, the process must be
  //    terminated and relaunched -> call abort on the host-side.
  // 3. Any other error code -> throw a runtime error.
  switch (e) {
    case cudaSuccess: break;
    case cudaErrorIllegalAddress:
    case cudaErrorAssert:
    case cudaErrorHardwareStackError:
    case cudaErrorIllegalInstruction:
    case cudaErrorMisalignedAddress:
    case cudaErrorInvalidAddressSpace:
    case cudaErrorInvalidPc:
    case cudaErrorLaunchFailure:
      cuda_internal_error_abort(e, name, file, line);
      break;
    default: cuda_internal_error_throw(e, name, file, line); break;
  }
}

#define KOKKOS_IMPL_CUDA_SAFE_CALL(call) \
  Kokkos::Impl::cuda_internal_safe_call(call, #call, __FILE__, __LINE__)

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_CUDA
#endif  // KOKKOS_CUDA_ERROR_HPP
