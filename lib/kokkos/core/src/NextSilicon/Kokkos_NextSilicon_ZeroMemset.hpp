// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICONSPACE_ZEROMEMSET_HPP
#define KOKKOS_NEXTSILICONSPACE_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <NextSilicon/Kokkos_NextSiliconSpace.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

#include <nextapi/memory.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct ZeroMemset<Kokkos::Experimental::NextSilicon> {
  ZeroMemset(const Kokkos::Experimental::NextSilicon& /*exec_space*/, void* dst,
             size_t cnt) {
    uint8_t pattern = 0;
    nextapi_memory_fill(dst, &pattern, sizeof(pattern), cnt);
  }
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_NEXTSILICONSPACE_ZEROMEMSET_HPP
