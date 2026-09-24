// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_PRINTF_HPP
#define KOKKOS_NEXTSILICON_PRINTF_HPP

#include <Kokkos_Macros.hpp>

#include <cstdio>

namespace Kokkos {
namespace Impl {

template <typename... Args>
void nextsilicon_printf(const char* format, Args... args) {
  // FIXME_NEXTSILICON: CS-515 tracks printing from device
  if (!__next_is_in_handed_off_code()) {
    if constexpr (sizeof...(Args) == 0)
      ::printf("%s", format);
    else
      ::printf(format, args...);
  }
}

}  // namespace Impl
}  // namespace Kokkos

#endif
