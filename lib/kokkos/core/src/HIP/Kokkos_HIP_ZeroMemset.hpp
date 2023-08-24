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
#ifndef KOKKOS_HIP_ZEROMEMSET_HPP
#define KOKKOS_HIP_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos {
namespace Impl {

template <class T, class... P>
struct ZeroMemset<HIP, View<T, P...>> {
  ZeroMemset(const HIP& exec_space, const View<T, P...>& dst,
             typename View<T, P...>::const_value_type&) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipMemsetAsync(
        dst.data(), 0, dst.size() * sizeof(typename View<T, P...>::value_type),
        exec_space.hip_stream()));
  }

  ZeroMemset(const View<T, P...>& dst,
             typename View<T, P...>::const_value_type&) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipMemset(dst.data(), 0,
                  dst.size() * sizeof(typename View<T, P...>::value_type)));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_HIP_ZEROMEMSET_HPP)
