// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <nextapi/memory.h>

#include <NextSilicon/Kokkos_NextSilicon.hpp>
#include <NextSilicon/Kokkos_NextSiliconSpace.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

#include <cstdlib>

namespace Kokkos {
namespace Experimental {

namespace {

// Guidance from NextSilicon is that typically we want each allocation backed by
// a few pages:
// - fewer/larger pages improves exception-based page migration and
// - fewer/larger pages reduces the number of required TLB entries, but...
// - fewer/larger pages increases the amount of over-allocation
// We choose to aim for at least 4 pages per allocation to have few but not too
// few
constexpr size_t pick_desired_page_size(size_t size) {
  // Maverick-2 supported page sizes in descending order
  constexpr size_t page_sizes[] = {
      16ULL << 30,   //  16 GiB
      4ULL << 30,    //   4 GiB
      1ULL << 30,    //   1 GiB
      256ULL << 20,  // 256 MiB
      64ULL << 20,   //  64 MiB
      16ULL << 20,   //  16 MiB
      4ULL << 20,    //   4 MiB
      1ULL << 20,    //   1 MiB
      256ULL << 10,  // 256 KiB
      64ULL << 10,   //  64 KiB
      16ULL << 10,   //  16 KiB
      4ULL << 10,    //   4 KiB
  };

  constexpr size_t min_pages = 4;

  // number of pages of size ps needed to cover allocation sz
  auto pages_needed = [](size_t sz, size_t ps) -> size_t {
    return (sz + ps - 1) / ps;
  };

  // largest page size that still provides min_pages
  for (auto page_size : page_sizes) {
    if (pages_needed(size, page_size) >= min_pages) return page_size;
  }

  // allow small allocations to share a page
  return 64;
}

static_assert(pick_desired_page_size(0) == 64);
static_assert(pick_desired_page_size(7) == 64);
static_assert(pick_desired_page_size(8'000'000) == 1 << 20);
static_assert(pick_desired_page_size(4'000'000'000) == 1 << 30);
static_assert(pick_desired_page_size(120'000'000'000) == 16ULL << 30);

}  //  namespace

void *NextSiliconSharedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void *NextSiliconSharedSpace::allocate(const char *arg_label,
                                       const size_t arg_alloc_size,
                                       const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

void *NextSiliconSharedSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  // The NextSilicon UVM migration runtime chooses a page size based on the
  // alignment of the allocation. Choose an alignment that corresponds to a
  // reasonable page size for the allocation size to improve the likelihood of
  // good performance from the UVM migration runtime.
  size_t alignment = pick_desired_page_size(arg_alloc_size);
  void *ptr        = std::aligned_alloc(alignment, arg_alloc_size);

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void NextSiliconSharedSpace::deallocate(void *const arg_alloc_ptr,
                                        const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size,
             /*arg_logical_size=*/0);
}

void NextSiliconSharedSpace::deallocate(const char *arg_label,
                                        void *const arg_alloc_ptr,
                                        const size_t arg_alloc_size,
                                        const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}

void NextSiliconSharedSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }

  std::free(arg_alloc_ptr);
}

}  // namespace Experimental
}  // namespace Kokkos
