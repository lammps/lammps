// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SERIAL_ZEROMEMSET_HPP
#define KOKKOS_SERIAL_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>
#include <Serial/Kokkos_Serial.hpp>

#include <type_traits>
#include <cstring>

namespace Kokkos {
namespace Impl {

template <>
struct ZeroMemset<Serial> {
  ZeroMemset(const Serial&, void* dst, size_t cnt) { std::memset(dst, 0, cnt); }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_SERIAL_ZEROMEMSET_HPP)
