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

#include <vector>

#ifdef KOKKOS_ENABLE_SYCL

namespace Kokkos {
namespace Impl {

template <class DT, class... DP>
struct ZeroMemset<Kokkos::Experimental::SYCL, DT, DP...> {
  ZeroMemset(const Kokkos::Experimental::SYCL& exec_space,
             const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    auto event = exec_space.impl_internal_space_instance()->m_queue->memset(
        dst.data(), 0,
        dst.size() * sizeof(typename View<DT, DP...>::value_type));
    exec_space.impl_internal_space_instance()
        ->m_queue->ext_oneapi_submit_barrier(std::vector<sycl::event>{event});
  }

  ZeroMemset(const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    Experimental::Impl::SYCLInternal::singleton().m_queue->memset(
        dst.data(), 0,
        dst.size() * sizeof(typename View<DT, DP...>::value_type));
  }
};

void DeepCopySYCL(void* dst, const void* src, size_t n);
void DeepCopyAsyncSYCL(const Kokkos::Experimental::SYCL& instance, void* dst,
                       const void* src, size_t n);
void DeepCopyAsyncSYCL(void* dst, const void* src, size_t n);

template <class MemSpace>
struct DeepCopy<MemSpace, HostSpace, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace>
struct DeepCopy<HostSpace, MemSpace, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2>
struct DeepCopy<MemSpace1, MemSpace2, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace1>::value &&
                                 is_sycl_type_space<MemSpace2>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2, class ExecutionSpace>
struct DeepCopy<
    MemSpace1, MemSpace2, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace1>::value &&
        is_sycl_type_space<MemSpace2>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace1::name() + "Space, " +
        MemSpace2::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<
    MemSpace, HostSpace, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace::name() +
        "Space, HostSpace, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<
    HostSpace, MemSpace, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<HostSpace, ") + MemSpace::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

}  // namespace Impl
}  // namespace Kokkos
#endif
#endif
