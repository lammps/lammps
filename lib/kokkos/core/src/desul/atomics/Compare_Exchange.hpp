/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_HPP_

#include "desul/atomics/Macros.hpp"

#include "desul/atomics/Compare_Exchange_ScopeCaller.hpp"

#ifdef DESUL_HAVE_GCC_ATOMICS
#include "desul/atomics/Compare_Exchange_GCC.hpp"
#endif
#ifdef DESUL_HAVE_MSVC_ATOMICS
#include "desul/atomics/Compare_Exchange_MSVC.hpp"
#endif
#ifdef DESUL_HAVE_SERIAL_ATOMICS
#include "desul/atomics/Compare_Exchange_Serial.hpp"
#endif
#ifdef DESUL_HAVE_CUDA_ATOMICS
#include "desul/atomics/Compare_Exchange_CUDA.hpp"
#endif
#ifdef DESUL_HAVE_HIP_ATOMICS
#include "desul/atomics/Compare_Exchange_HIP.hpp"
#endif
#ifdef DESUL_HAVE_OPENMP_ATOMICS
#include "desul/atomics/Compare_Exchange_OpenMP.hpp"
#endif
#ifdef DESUL_HAVE_SYCL_ATOMICS
#include "desul/atomics/Compare_Exchange_SYCL.hpp"
#endif
#endif
