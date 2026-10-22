/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_MSVC_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_MSVC_HPP_

#include <desul/atomics/Compare_Exchange_MSVC.hpp>
#include <desul/atomics/Lock_Free_Types_MSVC.hpp>
#include <desul/atomics/Macros.hpp>

namespace desul {
namespace Impl {

// MSVC does only provide _InterlockedCompareExchange but no load and store intrinsics.
// If we can require c++20 we can switch to std::atomic_ref

DESUL_IMPL_ATOMIC_LOAD_AND_STORE_WITH_CAS(DESUL_IMPL_HOST_FUNCTION, host)

}  // namespace Impl
}  // namespace desul

#endif
