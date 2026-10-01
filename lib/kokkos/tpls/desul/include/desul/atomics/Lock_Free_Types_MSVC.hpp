/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_LOCK_FREE_TYPES_MSVC_HPP_
#define DESUL_LOCK_FREE_TYPES_MSVC_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

template <class T>
inline constexpr bool host_atomic_always_lock_free<T, void> = (sizeof(T) == 1) ||
                                                              (sizeof(T) == 2) ||
                                                              (sizeof(T) == 4) ||
                                                              (sizeof(T) == 8) ||
                                                              (sizeof(T) == 16);

}
}  // namespace desul
#endif
