// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_ZEROMEMSET_HPP
#define KOKKOS_SYCL_ZEROMEMSET_HPP

#include <impl/Kokkos_ZeroMemset_fwd.hpp>
#include <SYCL/Kokkos_SYCL.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct ZeroMemset<Kokkos::SYCL> {
  ZeroMemset(const Kokkos::SYCL& exec_space, void* dst, size_t cnt) {
    auto event = exec_space.sycl_queue().memset(dst, 0, cnt);
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
    exec_space.sycl_queue().ext_oneapi_submit_barrier(
        std::vector<sycl::event>{event});
#endif
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_SYCL_ZEROMEMSET_HPP)
