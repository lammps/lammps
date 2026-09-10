// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_ERROR_HPP
#define KOKKOS_HIP_ERROR_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>

#include <hip/hip_runtime.h>

namespace Kokkos {
namespace Impl {

void hip_internal_error_throw(hipError_t e, const char* name,
                              const char* file = nullptr, const int line = 0);

void hip_internal_error_abort(hipError_t e, const char* name,
                              const char* file = nullptr, const int line = 0);

inline void hip_internal_safe_call(hipError_t e, const char* name,
                                   const char* file = nullptr,
                                   const int line   = 0) {
  // 1. Success -> normal continuation.
  // 2. Error codes for which, to continue using HIP, the process must be
  //    terminated and relaunched -> call abort on the host-side.
  // 3. Any other error code -> throw a runtime error.
  switch (e) {
    case hipSuccess: break;
    case hipErrorInvalidValue:
    case hipErrorOutOfMemory:
    case hipErrorInitializationError:
    case hipErrorDeinitialized:
    case hipErrorInvalidConfiguration:
    case hipErrorInvalidSymbol:
    case hipErrorInvalidDevicePointer:
    case hipErrorInvalidMemcpyDirection:
    case hipErrorInsufficientDriver:
    case hipErrorMissingConfiguration:
    case hipErrorPriorLaunchFailure:
    case hipErrorInvalidDeviceFunction:
    case hipErrorNoDevice:
    case hipErrorInvalidDevice:
    case hipErrorInvalidContext:
    case hipErrorNoBinaryForGpu:
    case hipErrorInvalidSource:
    case hipErrorIllegalState:
    case hipErrorNotFound:
    case hipErrorIllegalAddress:
    case hipErrorLaunchOutOfResources:
    case hipErrorLaunchTimeOut:
    case hipErrorAssert:
    case hipErrorLaunchFailure:
    case hipErrorNotSupported:
    case hipErrorStreamCaptureUnsupported:
    case hipErrorCapturedEvent:
    case hipErrorGraphExecUpdateFailure:
    case hipErrorUnknown: hip_internal_error_abort(e, name, file, line); break;
    default: hip_internal_error_throw(e, name, file, line);
  }
}

}  // namespace Impl
}  // namespace Kokkos

#define KOKKOS_IMPL_HIP_SAFE_CALL(call) \
  Kokkos::Impl::hip_internal_safe_call(call, #call, __FILE__, __LINE__)

#endif
