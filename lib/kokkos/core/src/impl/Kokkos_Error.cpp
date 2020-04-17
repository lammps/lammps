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

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <ostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <impl/Kokkos_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void host_abort(const char *const message) {
  fwrite(message, 1, strlen(message), stderr);
  fflush(stderr);
  ::abort();
}

void throw_runtime_exception(const std::string &msg) {
  std::ostringstream o;
  o << msg;
  traceback_callstack(o);
  throw std::runtime_error(o.str());
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
    case AllocationMechanism::PosixMemAlign: o << "posix_memalign()."; break;
    case AllocationMechanism::PosixMMap: o << "POSIX mmap()."; break;
    case AllocationMechanism::IntelMMAlloc:
      o << "the Intel _mm_malloc() intrinsic.";
      break;
    case AllocationMechanism::CudaMalloc: o << "cudaMalloc()."; break;
    case AllocationMechanism::CudaMallocManaged:
      o << "cudaMallocManaged().";
      break;
    case AllocationMechanism::CudaHostAlloc: o << "cudaHostAlloc()."; break;
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
namespace Impl {

void traceback_callstack(std::ostream &msg) {
  msg << std::endl << "Traceback functionality not available" << std::endl;
}

}  // namespace Impl
}  // namespace Kokkos
