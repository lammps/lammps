// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
