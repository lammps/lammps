// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_SERIAL_HPP
#define KOKKOS_DECLARE_SERIAL_HPP

#if defined(KOKKOS_ENABLE_SERIAL)
#include <Serial/Kokkos_Serial.hpp>
#include <Serial/Kokkos_Serial_MDRangePolicy.hpp>
#include <Serial/Kokkos_Serial_ZeroMemset.hpp>
#endif

#endif
