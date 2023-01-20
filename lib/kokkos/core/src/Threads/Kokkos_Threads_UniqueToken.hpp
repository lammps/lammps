/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_THREADS_UNIQUETOKEN_HPP
#define KOKKOS_THREADS_UNIQUETOKEN_HPP

#include <Kokkos_UniqueToken.hpp>

namespace Kokkos {
namespace Experimental {

template <>
class UniqueToken<Threads, UniqueTokenScope::Instance> {
 private:
  using buffer_type = Kokkos::View<uint32_t *, Kokkos::HostSpace>;
  int m_count;
  buffer_type m_buffer_view;
  uint32_t volatile *m_buffer;

 public:
  using execution_space = Threads;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const & = execution_space()) noexcept
      : m_count(::Kokkos::Threads::impl_thread_pool_size()),
        m_buffer_view(buffer_type()),
        m_buffer(nullptr) {}

  UniqueToken(size_type max_size, execution_space const & = execution_space())
      : m_count(max_size > ::Kokkos::Threads::impl_thread_pool_size()
                    ? ::Kokkos::Threads::impl_thread_pool_size()
                    : max_size),
        m_buffer_view(
            max_size > ::Kokkos::Threads::impl_thread_pool_size()
                ? buffer_type()
                : buffer_type("UniqueToken::m_buffer_view",
                              ::Kokkos::Impl::concurrent_bitset::buffer_bound(
                                  m_count))),
        m_buffer(m_buffer_view.data()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept { return m_count; }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((
        if (m_buffer == nullptr) {
          return Threads::impl_thread_pool_rank();
        } else {
          const ::Kokkos::pair<int, int> result =
              ::Kokkos::Impl::concurrent_bitset::acquire_bounded(
                  m_buffer, m_count, ::Kokkos::Impl::clock_tic() % m_count);

          if (result.first < 0) {
            ::Kokkos::abort(
                "UniqueToken<Threads> failure to acquire tokens, no tokens "
                "available");
          }
          return result.first;
        }))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int i) const noexcept {
    KOKKOS_IF_ON_HOST((if (m_buffer != nullptr) {
      ::Kokkos::Impl::concurrent_bitset::release(m_buffer, i);
    }))

    KOKKOS_IF_ON_DEVICE(((void)i;))
  }
};

template <>
class UniqueToken<Threads, UniqueTokenScope::Global> {
 public:
  using execution_space = Threads;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const & = execution_space()) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept {
    KOKKOS_IF_ON_HOST((return Threads::impl_thread_pool_size();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((return Threads::impl_thread_pool_rank();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
