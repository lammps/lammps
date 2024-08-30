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

/// \file Kokkos_Vectorization.hpp
/// \brief Declaration and definition of Kokkos::Vectorization interface.
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_VECTORIZATION_HPP
#define KOKKOS_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Vectorization.hpp>
#elif defined(KOKKOS_ENABLE_HIP)
#include <HIP/Kokkos_HIP_Vectorization.hpp>
#endif

#endif
