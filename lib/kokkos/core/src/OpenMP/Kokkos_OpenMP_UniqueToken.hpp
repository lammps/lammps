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

#ifndef KOKKOS_OPENMP_UNIQUE_TOKEN_HPP
#define KOKKOS_OPENMP_UNIQUE_TOKEN_HPP

#include <Kokkos_UniqueToken.hpp>

namespace Kokkos::Experimental {
template <>
class UniqueToken<OpenMP, UniqueTokenScope::Instance> {
 public:
  using execution_space = OpenMP;
  using size_type       = int;

 private:
  using buffer_type = Kokkos::View<uint32_t*, Kokkos::HostSpace>;
  execution_space m_exec;
  size_type m_count;
  buffer_type m_buffer_view;
  uint32_t volatile* m_buffer;

 public:
  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& exec = execution_space()) noexcept
      : m_exec(exec),
        m_count(m_exec.impl_thread_pool_size()),
        m_buffer_view(buffer_type()),
        m_buffer(nullptr) {}

  UniqueToken(size_type max_size,
              execution_space const& exec = execution_space())
      : m_exec(exec),
        m_count(max_size),
        m_buffer_view("UniqueToken::m_buffer_view",
                      ::Kokkos::Impl::concurrent_bitset::buffer_bound(m_count)),
        m_buffer(m_buffer_view.data()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept {
    KOKKOS_IF_ON_HOST((return m_count;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST(
        (if (m_count >= m_exec.impl_thread_pool_size()) return m_exec
             .impl_thread_pool_rank();
         const ::Kokkos::pair<int, int> result =
             ::Kokkos::Impl::concurrent_bitset::acquire_bounded(
                 m_buffer, m_count, ::Kokkos::Impl::clock_tic() % m_count);

         if (result.first < 0) {
           ::Kokkos::abort(
               "UniqueToken<OpenMP> failure to acquire tokens, no tokens "
               "available");
         }

         return result.first;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int i) const noexcept {
    KOKKOS_IF_ON_HOST((if (m_count < m_exec.impl_thread_pool_size()) {
      ::Kokkos::Impl::concurrent_bitset::release(m_buffer, i);
    }))

    KOKKOS_IF_ON_DEVICE(((void)i;))
  }
};

template <>
class UniqueToken<OpenMP, UniqueTokenScope::Global> {
 public:
  using execution_space = OpenMP;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept {
    KOKKOS_IF_ON_HOST(
        (return Kokkos::Impl::OpenMPInternal::max_hardware_threads();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief acquire value such that 0 <= value < size()
  // FIXME this is wrong when using nested parallelism. In that case multiple
  // threads have the same thread ID.
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((return omp_get_thread_num();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};
}  // namespace Kokkos::Experimental

#endif
