// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HALF_HPP_
#define KOKKOS_HALF_HPP_
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_HALF
#endif

#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_Half_Conversion.hpp>
#endif
#ifdef KOKKOS_ENABLE_HIP
#include <HIP/Kokkos_HIP_Half_Conversion.hpp>
#endif
#ifdef KOKKOS_ENABLE_SYCL
#include <SYCL/Kokkos_SYCL_Half_Conversion.hpp>
#endif
#include <impl/Kokkos_Half_NumericTraits.hpp>
#include <impl/Kokkos_Half_ReductionIdentity.hpp>
#include <impl/Kokkos_Half_MathematicalFunctions.hpp>

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_HALF
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_HALF
#endif
#endif  // KOKKOS_HALF_HPP_
