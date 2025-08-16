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

#ifndef KOKKOS_SERIAL_ZEROMEMSET_HPP
#define KOKKOS_SERIAL_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>
#include <Serial/Kokkos_Serial.hpp>

#include <type_traits>
#include <cstring>

namespace Kokkos {
namespace Impl {

// We only need to provide a specialization for Serial if there is a host
// parallel execution space since the specialization for
// DefaultHostExecutionSpace is defined elsewhere.
struct DummyExecutionSpace;
template <>
struct ZeroMemset<
    std::conditional_t<!std::is_same_v<Serial, DefaultHostExecutionSpace>,
                       Serial, DummyExecutionSpace>> {
  ZeroMemset(const Serial&, void* dst, size_t cnt) { std::memset(dst, 0, cnt); }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_SERIAL_ZEROMEMSET_HPP)
