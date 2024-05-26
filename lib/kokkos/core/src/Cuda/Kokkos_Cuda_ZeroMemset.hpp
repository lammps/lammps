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
#ifndef KOKKOS_CUDA_ZEROMEMSET_HPP
#define KOKKOS_CUDA_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <Cuda/Kokkos_Cuda.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos {
namespace Impl {

template <class T, class... P>
struct ZeroMemset<Kokkos::Cuda, View<T, P...>> {
  ZeroMemset(const Kokkos::Cuda& exec_space_instance,
             const View<T, P...>& dst) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (exec_space_instance.impl_internal_space_instance()
             ->cuda_memset_async_wrapper(
                 dst.data(), 0,
                 dst.size() * sizeof(typename View<T, P...>::value_type))));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_CUDA_ZEROMEMSET_HPP)
