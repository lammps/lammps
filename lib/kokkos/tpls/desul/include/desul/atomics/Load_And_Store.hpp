/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_HPP_

#include <desul/atomics/Macros.hpp>

#ifdef DESUL_HAVE_GCC_ATOMICS
#include <desul/atomics/Load_And_Store_GCC.hpp>
#endif
#ifdef DESUL_HAVE_MSVC_ATOMICS
#include <desul/atomics/Load_And_Store_MSVC.hpp>
#endif
#ifdef DESUL_HAVE_CUDA_ATOMICS
#include <desul/atomics/Load_And_Store_CUDA.hpp>
#endif
#ifdef DESUL_HAVE_HIP_ATOMICS
#include <desul/atomics/Load_And_Store_HIP.hpp>
#endif
#ifdef DESUL_HAVE_OPENMP_ATOMICS
#include <desul/atomics/Load_And_Store_OpenMP.hpp>
#endif
#ifdef DESUL_HAVE_OPENACC_ATOMICS
#include <desul/atomics/Load_And_Store_OpenACC.hpp>
#endif
#ifdef DESUL_HAVE_SYCL_ATOMICS
#include <desul/atomics/Load_And_Store_SYCL.hpp>
#endif

#include <desul/atomics/Load_And_Store_ScopeCaller.hpp>

#endif
