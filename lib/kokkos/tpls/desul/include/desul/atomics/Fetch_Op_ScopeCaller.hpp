/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_FETCH_OP_SCOPECALLER_HPP_
#define DESUL_ATOMICS_FETCH_OP_SCOPECALLER_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Macros.hpp>

namespace desul {
namespace Impl {

#define DESUL_IMPL_ATOMIC_FETCH_OPER(ANNOTATION, HOST_OR_DEVICE) \
  template <class Oper, class T, class MemoryOrder>              \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch_oper(               \
      const Oper& op,                                            \
      T* const dest,                                             \
      dont_deduce_this_parameter_t<const T> val,                 \
      MemoryOrder /*order*/,                                     \
      MemoryScopeCaller /*scope*/) {                             \
    T oldval = *dest;                                            \
    *dest = op.apply(oldval, val);                               \
    return oldval;                                               \
  }                                                              \
                                                                 \
  template <class Oper, class T, class MemoryOrder>              \
  ANNOTATION T HOST_OR_DEVICE##_atomic_oper_fetch(               \
      const Oper& op,                                            \
      T* const dest,                                             \
      dont_deduce_this_parameter_t<const T> val,                 \
      MemoryOrder /*order*/,                                     \
      MemoryScopeCaller /*scope*/) {                             \
    T oldval = *dest;                                            \
    T newval = op.apply(oldval, val);                            \
    *dest = newval;                                              \
    return newval;                                               \
  }

DESUL_IMPL_ATOMIC_FETCH_OPER(DESUL_IMPL_HOST_FUNCTION, host)
DESUL_IMPL_ATOMIC_FETCH_OPER(DESUL_IMPL_DEVICE_FUNCTION, device)

#undef DESUL_IMPL_ATOMIC_FETCH_OPER

}  // namespace Impl
}  // namespace desul

// FIXME consider implementing directly atomic_fetch_##OP and atomic_##OP##_fetch or
// dropping this placeholder

#endif
