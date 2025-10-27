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
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_MEMORYTRAITS_HPP
#define KOKKOS_MEMORYTRAITS_HPP

#include <impl/Kokkos_Traits.hpp>
#include <Kokkos_BitManipulation.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Memory access traits for views, an extension point.
 *
 *  These traits should be orthogonal.  If there are dependencies then
 *  the MemoryTraits template must detect and enforce dependencies.
 *
 *  A zero value is the default for a View, indicating that none of
 *  these traits are present.
 */
enum MemoryTraitsFlags {
  Unmanaged    = 0x01,
  RandomAccess = 0x02,
  Atomic       = 0x04,
  Restrict     = 0x08,
  Aligned      = 0x10
};

template <unsigned T = 0>
struct MemoryTraits {
  //! Tag this class as a kokkos memory traits:
  using memory_traits = MemoryTraits<T>;

  static constexpr unsigned impl_value = T;

  static constexpr bool is_unmanaged =
      (unsigned(0) != (T & unsigned(Kokkos::Unmanaged)));
  static constexpr bool is_random_access =
      (unsigned(0) != (T & unsigned(Kokkos::RandomAccess)));
  static constexpr bool is_atomic =
      (unsigned(0) != (T & unsigned(Kokkos::Atomic)));
  static constexpr bool is_restrict =
      (unsigned(0) != (T & unsigned(Kokkos::Restrict)));
  static constexpr bool is_aligned =
      (unsigned(0) != (T & unsigned(Kokkos::Aligned)));
};

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
using MemoryManaged KOKKOS_DEPRECATED = Kokkos::MemoryTraits<>;
#endif
using MemoryUnmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
using MemoryRandomAccess =
    Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>;

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief Memory alignment settings
 *
 *  Sets global value for memory alignment.  Must be a power of two!
 *  Enable compatibility of views from different devices with static stride.
 *  Use compiler flag to enable overwrites.
 */
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
static constexpr unsigned MEMORY_ALIGNMENT = KOKKOS_IMPL_MEMORY_ALIGNMENT;
static constexpr unsigned MEMORY_ALIGNMENT_THRESHOLD =
    KOKKOS_IMPL_MEMORY_ALIGNMENT_THRESHOLD;
#else
static constexpr unsigned MEMORY_ALIGNMENT           = 64;
static constexpr unsigned MEMORY_ALIGNMENT_THRESHOLD = 1;
#endif
static_assert(has_single_bit(MEMORY_ALIGNMENT),
              "MEMORY_ALIGNMENT must be a power of 2");

// ------------------------------------------------------------------ //
//  this identifies the default memory trait
//
template <typename Tp>
struct is_default_memory_trait : std::false_type {};

template <>
struct is_default_memory_trait<Kokkos::MemoryTraits<>> : std::true_type {};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_MEMORYTRAITS_HPP */
