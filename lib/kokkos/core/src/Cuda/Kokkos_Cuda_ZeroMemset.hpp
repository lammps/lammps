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

template <>
struct ZeroMemset<Kokkos::Cuda> {
  ZeroMemset(const Kokkos::Cuda& exec_space_instance, void* dst, size_t cnt) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (exec_space_instance.impl_internal_space_instance()
             ->cuda_memset_async_wrapper(dst, 0, cnt)));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_CUDA_ZEROMEMSET_HPP)
