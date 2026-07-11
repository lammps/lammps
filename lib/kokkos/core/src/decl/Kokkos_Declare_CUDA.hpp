// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_DECLARE_CUDA_HPP
#define KOKKOS_DECLARE_CUDA_HPP

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_Half_Impl_Type.hpp>
#include <Cuda/Kokkos_Cuda_Half_MathematicalFunctions.hpp>
#include <Cuda/Kokkos_Cuda_Half_Conversion.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_MDRange.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_Range.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_Team.hpp>
#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>
#include <Cuda/Kokkos_Cuda_View.hpp>
#include <Cuda/Kokkos_Cuda_Team.hpp>
#include <Cuda/Kokkos_Cuda_MDRangePolicy.hpp>
#include <Cuda/Kokkos_Cuda_UniqueToken.hpp>
#include <Cuda/Kokkos_Cuda_ZeroMemset.hpp>
#endif

#endif
