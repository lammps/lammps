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
#ifndef KOKKOS_ANONYMOUSSPACE_HPP
#define KOKKOS_ANONYMOUSSPACE_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>
#include <cstddef>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class AnonymousSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space    = AnonymousSpace;
  using execution_space = Kokkos::DefaultExecutionSpace;
  using size_type       = size_t;

  //! This memory space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  /**\brief  Default memory space instance */
  AnonymousSpace()                          = default;
  AnonymousSpace(AnonymousSpace &&rhs)      = default;
  AnonymousSpace(const AnonymousSpace &rhs) = default;
  AnonymousSpace &operator=(AnonymousSpace &&) = default;
  AnonymousSpace &operator=(const AnonymousSpace &) = default;
  ~AnonymousSpace()                                 = default;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char *name() { return "Anonymous"; }
};

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <typename OtherSpace>
struct MemorySpaceAccess<Kokkos::AnonymousSpace, OtherSpace> {
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <typename OtherSpace>
struct MemorySpaceAccess<OtherSpace, Kokkos::AnonymousSpace> {
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::AnonymousSpace, Kokkos::AnonymousSpace> {
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

}  // namespace Impl

}  // namespace Kokkos

#endif  // #define KOKKOS_ANONYMOUSSPACE_HPP
