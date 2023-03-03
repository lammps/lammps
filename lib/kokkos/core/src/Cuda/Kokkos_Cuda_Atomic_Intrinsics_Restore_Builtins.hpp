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

#ifdef KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED

#undef __atomic_load_n
#undef __atomic_load
#undef __atomic_store_n
#undef __atomic_store
#undef __atomic_exchange_n
#undef __atomic_exchange
#undef __atomic_compare_exchange_n
#undef __atomic_compare_exchange
#undef __atomic_fetch_add
#undef __atomic_fetch_sub
#undef __atomic_fetch_and
#undef __atomic_fetch_xor
#undef __atomic_fetch_or
#undef __atomic_test_and_set
#undef __atomic_clear
#undef __atomic_always_lock_free
#undef __atomic_is_lock_free
#undef __atomic_thread_fence
#undef __atomic_signal_fence

#undef KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED

#endif  // KOKKOS_SIMT_ATOMIC_BUILTIN_REPLACEMENTS_DEFINED
