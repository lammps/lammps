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

#ifndef KOKKOS_SYCLDEEPCOPY_HPP
#define KOKKOS_SYCLDEEPCOPY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_SYCL.hpp>

#ifdef KOKKOS_ENABLE_SYCL

namespace Kokkos {
namespace Impl {

template <>
struct DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                Kokkos::Experimental::SYCLDeviceUSMSpace,
                Kokkos::Experimental::SYCL> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::SYCL&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
                Kokkos::Experimental::SYCL> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::SYCL&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
                Kokkos::Experimental::SYCL> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::SYCL&, void* dst, const void* src,
           size_t);
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                Kokkos::Experimental::SYCLDeviceUSMSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                   Kokkos::Experimental::SYCLDeviceUSMSpace,
                   Kokkos::Experimental::SYCL>(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
             Kokkos::Experimental::SYCLDeviceUSMSpace,
             Kokkos::Experimental::SYCL>(Kokkos::Experimental::SYCL(), dst, src,
                                         n);
    Kokkos::Experimental::SYCL().fence();
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
                   Kokkos::Experimental::SYCL>(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
             Kokkos::Experimental::SYCL>(Kokkos::Experimental::SYCL(), dst, src,
                                         n);
    Kokkos::Experimental::SYCL().fence();
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
                   Kokkos::Experimental::SYCL>(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
             Kokkos::Experimental::SYCL>(Kokkos::Experimental::SYCL(), dst, src,
                                         n);
    Kokkos::Experimental::SYCL().fence();
  }
};

}  // namespace Impl
}  // namespace Kokkos
#endif
#endif
