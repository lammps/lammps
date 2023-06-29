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

#ifndef KOKKOS_OPENMPTARGET_UNIQUE_TOKEN_HPP
#define KOKKOS_OPENMPTARGET_UNIQUE_TOKEN_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_OPENMPTARGET

#include <OpenMPTarget/Kokkos_OpenMPTargetSpace.hpp>
#include <Kokkos_UniqueToken.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_ConcurrentBitset.hpp>

namespace Kokkos {
namespace Experimental {

// both global and instance Unique Tokens are implemented in the same way
template <>
class UniqueToken<OpenMPTarget, UniqueTokenScope::Global> {
 protected:
  uint32_t volatile* m_buffer;
  uint32_t m_count;

 public:
  using execution_space = OpenMPTarget;
  using size_type       = int32_t;

  explicit UniqueToken(execution_space const& = execution_space());

  KOKKOS_DEFAULTED_FUNCTION
  UniqueToken(const UniqueToken&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  UniqueToken(UniqueToken&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  UniqueToken& operator=(const UniqueToken&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  UniqueToken& operator=(UniqueToken&&) = default;

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type size() const noexcept { return m_count; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type acquire() const {
    const Kokkos::pair<int, int> result =
        Kokkos::Impl::concurrent_bitset::acquire_bounded(
            m_buffer, m_count, Kokkos::Impl::clock_tic() % m_count);

    if (result.first < 0) {
      Kokkos::abort(
          "UniqueToken<OpenMPTarget> failure to acquire tokens, no tokens "
          "available");
    }

    return result.first;
  }

  /// \brief release an acquired value
  KOKKOS_INLINE_FUNCTION
  void release(size_type i) const noexcept {
    Kokkos::Impl::concurrent_bitset::release(m_buffer, i);
  }
};

template <>
class UniqueToken<OpenMPTarget, UniqueTokenScope::Instance>
    : public UniqueToken<OpenMPTarget, UniqueTokenScope::Global> {
 private:
  Kokkos::View<uint32_t*, ::Kokkos::Experimental::OpenMPTargetSpace>
      m_buffer_view;

 public:
  explicit UniqueToken(execution_space const& arg = execution_space())
      : UniqueToken<OpenMPTarget, UniqueTokenScope::Global>(arg) {}

  UniqueToken(size_type max_size, execution_space const& = execution_space())
      : m_buffer_view(
            "Kokkos::UniqueToken::m_buffer_view",
            ::Kokkos::Impl::concurrent_bitset::buffer_bound(max_size)) {
    m_buffer = m_buffer_view.data();
    m_count  = max_size;
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_OPENMPTARGET
#endif  // KOKKOS_OPENMPTARGET_UNIQUE_TOKEN_HPP
