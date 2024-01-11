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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <cstring>
#include <cstdlib>

#include <ostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <impl/Kokkos_Error.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void throw_runtime_exception(const std::string &msg) {
  throw std::runtime_error(msg);
}

std::string human_memory_size(size_t arg_bytes) {
  double bytes   = arg_bytes;
  const double K = 1024;
  const double M = K * 1024;
  const double G = M * 1024;

  std::ostringstream out;
  if (bytes < K) {
    out << std::setprecision(4) << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << std::setprecision(4) << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << std::setprecision(4) << bytes << " M";
  } else {
    bytes /= G;
    out << std::setprecision(4) << bytes << " G";
  }
  return out.str();
}

}  // namespace Impl

void Experimental::RawMemoryAllocationFailure::print_error_message(
    std::ostream &o) const {
  o << "Allocation of size " << Impl::human_memory_size(m_attempted_size);
  o << " failed";
  switch (m_failure_mode) {
    case FailureMode::OutOfMemoryError:
      o << ", likely due to insufficient memory.";
      break;
    case FailureMode::AllocationNotAligned:
      o << " because the allocation was improperly aligned.";
      break;
    case FailureMode::InvalidAllocationSize:
      o << " because the requested allocation size is not a valid size for the"
           " requested allocation mechanism (it's probably too large).";
      break;
    // TODO move this to the subclass for Cuda-related things
    case FailureMode::MaximumCudaUVMAllocationsExceeded:
      o << " because the maximum Cuda UVM allocations was exceeded.";
      break;
    case FailureMode::Unknown: o << " because of an unknown error."; break;
  }
  o << "  (The allocation mechanism was ";
  switch (m_mechanism) {
    case AllocationMechanism::StdMalloc: o << "standard malloc()."; break;
    case AllocationMechanism::CudaMalloc: o << "cudaMalloc()."; break;
    case AllocationMechanism::CudaMallocManaged:
      o << "cudaMallocManaged().";
      break;
    case AllocationMechanism::CudaHostAlloc: o << "cudaHostAlloc()."; break;
    case AllocationMechanism::HIPMalloc: o << "hipMalloc()."; break;
    case AllocationMechanism::HIPHostMalloc: o << "hipHostMalloc()."; break;
    case AllocationMechanism::HIPMallocManaged:
      o << "hipMallocManaged().";
      break;
    case AllocationMechanism::SYCLMallocDevice:
      o << "sycl::malloc_device().";
      break;
    case AllocationMechanism::SYCLMallocShared:
      o << "sycl::malloc_shared().";
      break;
    case AllocationMechanism::SYCLMallocHost:
      o << "sycl::malloc_host().";
      break;
    default: o << "unsupported.";
  }
  append_additional_error_information(o);
  o << ")" << std::endl;
}

std::string Experimental::RawMemoryAllocationFailure::get_error_message()
    const {
  std::ostringstream out;
  print_error_message(out);
  return out.str();
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

#ifdef KOKKOS_ENABLE_CUDA
namespace Experimental {

void CudaRawMemoryAllocationFailure::append_additional_error_information(
    std::ostream &o) const {
  if (m_error_code != cudaSuccess) {
    o << "  The Cuda allocation returned the error code \""
      << cudaGetErrorName(m_error_code) << "\".";
  }
}

}  // end namespace Experimental
#endif

}  // namespace Kokkos
