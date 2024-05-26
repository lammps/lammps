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

#ifndef KOKKOS_DECLARE_SERIAL_HPP
#define KOKKOS_DECLARE_SERIAL_HPP

#if defined(KOKKOS_ENABLE_SERIAL)
#include <Serial/Kokkos_Serial.hpp>
#include <Serial/Kokkos_Serial_MDRangePolicy.hpp>
#include <Serial/Kokkos_Serial_ZeroMemset.hpp>
#endif

#endif
