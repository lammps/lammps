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

#ifndef KOKKOS_COPY_OPS_FOR_BINSORT_IMPL_HPP_
#define KOKKOS_COPY_OPS_FOR_BINSORT_IMPL_HPP_

#include <Kokkos_Macros.hpp>
#include <cstddef>

namespace Kokkos {
namespace Impl {

template <class DstViewType, class SrcViewType, int Rank = DstViewType::rank>
struct CopyOp;

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 1> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    dst(i_dst) = src(i_src);
  }
};

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 2> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    for (int j = 0; j < (int)dst.extent(1); j++) dst(i_dst, j) = src(i_src, j);
  }
};

template <class DstViewType, class SrcViewType>
struct CopyOp<DstViewType, SrcViewType, 3> {
  KOKKOS_INLINE_FUNCTION
  static void copy(DstViewType const& dst, size_t i_dst, SrcViewType const& src,
                   size_t i_src) {
    for (int j = 0; j < dst.extent(1); j++)
      for (int k = 0; k < dst.extent(2); k++)
        dst(i_dst, j, k) = src(i_src, j, k);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
