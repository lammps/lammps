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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_VIEW_CHECKING_HPP
#define KOKKOS_VIEW_CHECKING_HPP

#include <cstring>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <View/MDSpan/Kokkos_MDSpan_Header.hpp>

namespace Kokkos::Impl {
// primary template: memory space is accessible, do nothing.
template <class MemorySpace, class AccessSpace,
          bool = SpaceAccessibility<AccessSpace, MemorySpace>::accessible>
struct RuntimeCheckBasicViewMemoryAccessViolation {
  KOKKOS_FUNCTION RuntimeCheckBasicViewMemoryAccessViolation(
      Kokkos::Impl::SharedAllocationTracker const &) {}
};

// explicit specialization: memory access violation will occur, call abort with
// the specified error message.
template <class MemorySpace, class AccessSpace>
struct RuntimeCheckBasicViewMemoryAccessViolation<MemorySpace, AccessSpace,
                                                  false> {
  KOKKOS_FUNCTION RuntimeCheckBasicViewMemoryAccessViolation(
      Kokkos::Impl::SharedAllocationTracker const &tracker) {
    char err[256] =
        "Kokkos::View ERROR: attempt to access inaccessible memory space "
        "(label=\"";

    KOKKOS_IF_ON_HOST(({
      if (tracker.has_record()) {
        strncat(err, tracker.template get_label<void>().c_str(), 128);
      } else {
        strcat(err, "**UNMANAGED**");
      }
    }))

    KOKKOS_IF_ON_DEVICE(({
      strcat(err, "**UNAVAILABLE**");
      (void)tracker;
    }))

    strcat(err, "\")");

    Kokkos::abort(err);
  }
};

template <class MemorySpace>
KOKKOS_FUNCTION void runtime_check_memory_access_violation(
    SharedAllocationTracker const &track) {
  KOKKOS_IF_ON_HOST(((void)RuntimeCheckBasicViewMemoryAccessViolation<
                         MemorySpace, DefaultHostExecutionSpace>(track);))
  KOKKOS_IF_ON_DEVICE(((void)RuntimeCheckBasicViewMemoryAccessViolation<
                           MemorySpace, DefaultExecutionSpace>(track);))
}

template <class IndexType, std::size_t... Extents, class... Indices,
          std::size_t... Enumerate>
KOKKOS_FUNCTION bool within_range(
    Kokkos::extents<IndexType, Extents...> const &exts,
    std::index_sequence<Enumerate...>, Indices... indices) {
  // FIXME[CUDA11]: This is written so weirdly to avoid warnings with CUDA 11
  // Without the workaround this could just be written as:
  // return ((indices < exts.extent(Enumerate)) && ...) &&
  //     ((std::is_unsigned_v<decltype(indices)> ||
  //       (indices >= static_cast<decltype(indices)>(0))) &&
  //      ...);
  [[maybe_unused]] auto check_index_min = [](auto idx) {
    if constexpr (!std::is_unsigned_v<decltype(idx)>) {
      return idx >= static_cast<decltype(idx)>(0);
    }

    return true;
  };

  return ((static_cast<size_t>(indices) <
           static_cast<size_t>(exts.extent(Enumerate))) &&
          ...) &&
         (check_index_min(indices) && ...);
}

template <class... Indices>
KOKKOS_FUNCTION constexpr char *append_formatted_multidimensional_index(
    char *dest, Indices... indices) {
  char *d = dest;
  strcat(d, "[");
  (
      [&] {
        d += strlen(d);
        to_chars_i(d,
                   d + 20,  // 20 digits ought to be enough
                   indices);
        strcat(d, ",");
      }(),
      ...);
  d[strlen(d) - 1] = ']';  // overwrite trailing comma
  return dest;
}

template <class IndexType, size_t... Extents, std::size_t... Enumerate>
KOKKOS_FUNCTION void print_extents(
    char *dest, Kokkos::extents<IndexType, Extents...> const &exts,
    std::index_sequence<Enumerate...>) {
  append_formatted_multidimensional_index(dest, exts.extent(Enumerate)...);
}

template <class ExtentsType, class... IndexTypes>
KOKKOS_INLINE_FUNCTION void view_verify_operator_bounds(
    SharedAllocationTracker const &tracker, const ExtentsType &exts,
    [[maybe_unused]] const void *data, IndexTypes... idx) {
  using idx_t = typename ExtentsType::index_type;
  if (!within_range(exts, std::make_index_sequence<sizeof...(IndexTypes)>(),
                    idx...)) {
    char err[256] = "";
    strcat(err, "Kokkos::View ERROR: out of bounds access");
    strcat(err, " label=(\"");
    KOKKOS_IF_ON_HOST(
        if (tracker.has_record()) {
          strncat(err, tracker.template get_label<void>().c_str(), 128);
        } else { strcat(err, "**UNMANAGED**"); })
    KOKKOS_IF_ON_DEVICE([&] {
      if (!tracker.has_record()) {
        strcat(err, "**UNMANAGED**");
        return;
      }
      SharedAllocationHeader const *const header =
          SharedAllocationHeader::get_header(data);
      char const *const label = header->label();
      strcat(err, label);
    }();)
    strcat(err, "\") with indices ");
    append_formatted_multidimensional_index(err, static_cast<idx_t>(idx)...);
    strcat(err, " but extents ");
    print_extents(err, exts, std::make_index_sequence<sizeof...(IndexTypes)>());
    Kokkos::abort(err);
  }
}
}  // namespace Kokkos::Impl

#endif  // KOKKOS_VIEW_CHECKING_HPP
