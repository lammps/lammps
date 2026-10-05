/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_OPENMP_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_OPENMP_HPP_

#include <desul/atomics/Compare_Exchange_OpenMP.hpp>
#include <desul/atomics/Lock_Free_Types_OpenMP.hpp>
#include <desul/atomics/Macros.hpp>

namespace desul {
namespace Impl {

DESUL_IMPL_ATOMIC_LOAD_AND_STORE_WITH_CAS(DESUL_IMPL_HOST_FUNCTION, host)

}  // namespace Impl
}  // namespace desul

#endif
