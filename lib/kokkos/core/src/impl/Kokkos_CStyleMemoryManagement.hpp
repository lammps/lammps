// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_C_STYLE_MEMORY_MANAGEMENT_HPP
#define KOKKOS_C_STYLE_MEMORY_MANAGEMENT_HPP

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_C_STYLE_MEMORY_MANAGEMENT
#endif

#include <Kokkos_Abort.hpp>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_InitializeFinalize.hpp>

#include <KokkosCore_Config_DeclareBackend.hpp>  // FIXME

#include <sstream>

namespace Kokkos::Impl {

inline void check_init_final([[maybe_unused]] char const* func_name) {
// FIXME_THREADS: Checking for calls to kokkos_malloc, kokkos_realloc,
// kokkos_free before initialize or after finalize is currently disabled
// for the Threads backend. Refer issue #7944.
#if !defined(KOKKOS_ENABLE_THREADS)
  if (is_finalized()) {
    std::stringstream ss;
    ss << "Kokkos ERROR: attempting to perform C-style memory management "
          "via ";
    ss << func_name << "() **after** Kokkos::finalize() was called\n";
    Kokkos::abort(ss.str().c_str());
  } else if (!is_initialized()) {
    std::stringstream ss;
    ss << "Kokkos ERROR: attempting to perform C-style memory management "
          "via ";
    ss << func_name << "() **before** Kokkos::initialize() was called\n";
    Kokkos::abort(ss.str().c_str());
  }
#endif
}

}  // namespace Kokkos::Impl

namespace Kokkos {

/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const std::string& arg_alloc_label,
                           const size_t arg_alloc_size) {
  Impl::check_init_final("kokkos_malloc");
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), arg_alloc_label, arg_alloc_size);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const size_t arg_alloc_size) {
  Impl::check_init_final("kokkos_malloc");
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), "no-label", arg_alloc_size);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void kokkos_free(void* arg_alloc) {
  Impl::check_init_final("kokkos_free");
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::deallocate_tracked(
      arg_alloc);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_realloc(void* arg_alloc, const size_t arg_alloc_size) {
  Impl::check_init_final("kokkos_realloc");
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::reallocate_tracked(
      arg_alloc, arg_alloc_size);
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_C_STYLE_MEMORY_MANAGEMENT
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_C_STYLE_MEMORY_MANAGEMENT
#endif

#endif
