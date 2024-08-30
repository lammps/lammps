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

#ifndef KOKKOS_SERIAL_UNIQUE_TOKEN_HPP
#define KOKKOS_SERIAL_UNIQUE_TOKEN_HPP

#include <Kokkos_UniqueToken.hpp>

namespace Kokkos {
namespace Experimental {

template <>
class UniqueToken<Serial, UniqueTokenScope::Instance> {
 public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief create object size for requested size on given instance
  ///
  /// It is the users responsibility to only acquire size tokens concurrently
  UniqueToken(size_type, execution_space const& = execution_space()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

template <>
class UniqueToken<Serial, UniqueTokenScope::Global> {
 public:
  using execution_space = Serial;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return 1; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept { return 0; }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
