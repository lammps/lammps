/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_LOCK_FREE_TYPES_GCC_HPP_
#define DESUL_LOCK_FREE_TYPES_GCC_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

template <class T>
inline constexpr bool host_atomic_always_lock_free<T, void> =
#ifndef DESUL_HAVE_LIBATOMIC
    ((sizeof(T) == 1 && alignof(T) == 1) || (sizeof(T) == 2 && alignof(T) == 2) ||
     (sizeof(T) == 4 && alignof(T) == 4) || (sizeof(T) == 8 && alignof(T) == 8)
#ifdef DESUL_HAVE_16BYTE_LOCK_FREE_ATOMICS_HOST
     || (sizeof(T) == 16 && alignof(T) == 16)
#endif
         ) &&
#endif
    std::is_trivially_copyable<T>::value;

}  // namespace Impl
}  // namespace desul
#endif
