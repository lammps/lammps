/* 
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_OPENMP40_HPP_
#define DESUL_ATOMICS_OPENMP40_HPP_
#include<type_traits>

namespace desul {
namespace Impl {
  template<class MEMORY_ORDER_TMP, class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_pre_capture_flush(MEMORY_ORDER_TMP, MEMORY_SCOPE_TMP) {}
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_pre_capture_flush(MemoryOrderAcquire, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderAcquire(), MEMORY_SCOPE_TMP());
  }
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_pre_capture_flush(MemoryOrderAcqRel, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderAcqRel(), MEMORY_SCOPE_TMP());
  }
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_pre_capture_flush(MemoryOrderSeqCst, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderSeqCst(), MEMORY_SCOPE_TMP());
  }

  template<class MEMORY_ORDER_TMP, class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_post_capture_flush(MEMORY_ORDER_TMP, MEMORY_SCOPE_TMP) {}
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_post_capture_flush(MemoryOrderRelease, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderRelease(), MEMORY_SCOPE_TMP());
  }
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_post_capture_flush(MemoryOrderAcqRel, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderAcqRel(), MEMORY_SCOPE_TMP());
  }
  template<class MEMORY_SCOPE_TMP>
  void openmp_maybe_call_post_capture_flush(MemoryOrderSeqCst, MEMORY_SCOPE_TMP) {
    atomic_thread_fence(MemoryOrderSeqCst(), MEMORY_SCOPE_TMP());
  }

  template<class T>
  struct is_openmp_atomic_type_t {
    static constexpr bool value = std::is_arithmetic<T>::value;
  };
  template<class T>
  constexpr bool is_openmp_atomic_type_v = is_openmp_atomic_type_t<T>::value;
}
}

namespace desul {
// Can't use a macro approach to get all definitions since the ops include #pragma omp
// So gonna use multiple inclusion of the same code snippet here.

// Can't do Node level atomics this way with OpenMP Target, but we could 
// have a define which says whether or not Device level IS node level (e.g. for pure CPU node)

#define MEMORY_ORDER MemoryOrderRelaxed
// #define MEMORY_SCOPE MemoryScopeNode
// #include<desul/atomics/openmp/OpenMP_40_op.inc>
// #undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeDevice
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeCore
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#undef MEMORY_ORDER

#define MEMORY_ORDER MemoryOrderAcqRel
// #define MEMORY_SCOPE MemoryScopeNode
// #include<desul/atomics/openmp/OpenMP_40_op.inc>
// #undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeDevice
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeCore
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#undef MEMORY_ORDER

#define MEMORY_ORDER MemoryOrderSeqCst
// #define MEMORY_SCOPE MemoryScopeNode
// #include<desul/atomics/openmp/OpenMP_40_op.inc>
// #undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeDevice
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#define MEMORY_SCOPE MemoryScopeCore
#include<desul/atomics/openmp/OpenMP_40_op.inc>
#undef MEMORY_SCOPE
#undef MEMORY_ORDER
}  // namespace desul
#endif
