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

#ifndef KOKKOS_DECLARE_CUDA_HPP
#define KOKKOS_DECLARE_CUDA_HPP

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_Half_Impl_Type.hpp>
#include <Cuda/Kokkos_Cuda_Half_Conversion.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_MDRange.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_Range.hpp>
#include <Cuda/Kokkos_Cuda_Parallel_Team.hpp>
#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_Instance.hpp>
#include <Cuda/Kokkos_Cuda_View.hpp>
#include <Cuda/Kokkos_Cuda_Team.hpp>
#include <Cuda/Kokkos_Cuda_Task.hpp>
#include <Cuda/Kokkos_Cuda_MDRangePolicy.hpp>
#include <Cuda/Kokkos_Cuda_UniqueToken.hpp>
#include <Cuda/Kokkos_Cuda_ZeroMemset.hpp>
#endif

#endif
