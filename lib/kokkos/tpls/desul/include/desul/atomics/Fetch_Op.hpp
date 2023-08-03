/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_FETCH_OP_HPP_
#define DESUL_ATOMICS_FETCH_OP_HPP_

#include <desul/atomics/Macros.hpp>

#ifdef DESUL_HAVE_GCC_ATOMICS
#include <desul/atomics/Fetch_Op_GCC.hpp>
#endif
#ifdef DESUL_HAVE_CUDA_ATOMICS
#include <desul/atomics/Fetch_Op_CUDA.hpp>
#endif
#ifdef DESUL_HAVE_HIP_ATOMICS
#include <desul/atomics/Fetch_Op_HIP.hpp>
#endif
#ifdef DESUL_HAVE_OPENMP_ATOMICS
#include <desul/atomics/Fetch_Op_OpenMP.hpp>
#endif
#ifdef DESUL_HAVE_SYCL_ATOMICS
#include <desul/atomics/Fetch_Op_SYCL.hpp>
#endif

#include <desul/atomics/Fetch_Op_ScopeCaller.hpp>

// Must come last
#include <desul/atomics/Fetch_Op_Generic.hpp>

#endif
