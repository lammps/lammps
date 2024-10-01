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

#ifndef KOKKOS_DESUL_ATOMICS_CONFIG_HPP
#define KOKKOS_DESUL_ATOMICS_CONFIG_HPP

#if defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL)
#define DESUL_CUDA_ARCH_IS_PRE_PASCAL
#endif

#if defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL) || \
    defined(KOKKOS_ARCH_PASCAL)
#define DESUL_CUDA_ARCH_IS_PRE_VOLTA
#endif

#endif
