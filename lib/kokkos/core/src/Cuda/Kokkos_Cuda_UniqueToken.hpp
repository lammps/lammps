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

#ifndef KOKKOS_CUDA_UNIQUE_TOKEN_HPP
#define KOKKOS_CUDA_UNIQUE_TOKEN_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_UniqueToken.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_ConcurrentBitset.hpp>

namespace Kokkos {
namespace Experimental {

// both global and instance Unique Tokens are implemented in the same way
template <>
class UniqueToken<Cuda, UniqueTokenScope::Global> {
 private:
  uint32_t volatile* m_buffer;
  uint32_t m_count;

 public:
  using execution_space = Cuda;
  using size_type       = int32_t;

#if defined(KOKKOS_ENABLE_DEPRECATED_CODE)
  explicit UniqueToken(execution_space const&);

  KOKKOS_INLINE_FUNCTION
  UniqueToken() : m_buffer(0), m_count(0) {}
#else
  explicit UniqueToken(execution_space const& = execution_space());
#endif

#ifdef KOKKOS_CUDA_9_DEFAULTED_BUG_WORKAROUND
  KOKKOS_INLINE_FUNCTION
  UniqueToken(const UniqueToken& rhs)
      : m_buffer(rhs.m_buffer), m_count(rhs.m_count) {}

  KOKKOS_INLINE_FUNCTION
  UniqueToken(UniqueToken&& rhs)
      : m_buffer(std::move(rhs.m_buffer)), m_count(std::move(rhs.m_count)) {}

  KOKKOS_INLINE_FUNCTION
  UniqueToken& operator=(const UniqueToken& rhs) {
    m_buffer = rhs.m_buffer;
    m_count  = rhs.m_count;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  UniqueToken& operator=(UniqueToken&& rhs) {
    m_buffer = std::move(rhs.m_buffer);
    m_count  = std::move(rhs.m_count);
    return *this;
  }
#else
  KOKKOS_INLINE_FUNCTION
  UniqueToken(const UniqueToken&) = default;

  KOKKOS_INLINE_FUNCTION
  UniqueToken(UniqueToken&&) = default;

  KOKKOS_INLINE_FUNCTION
  UniqueToken& operator=(const UniqueToken&) = default;

  KOKKOS_INLINE_FUNCTION
  UniqueToken& operator=(UniqueToken&&) = default;
#endif

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
          "UniqueToken<Cuda> failure to release tokens, no tokens available");
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
class UniqueToken<Cuda, UniqueTokenScope::Instance>
    : public UniqueToken<Cuda, UniqueTokenScope::Global> {
 public:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  explicit UniqueToken(execution_space const& arg)
      : UniqueToken<Cuda, UniqueTokenScope::Global>(arg) {}
#else
  explicit UniqueToken(execution_space const& arg = execution_space())
      : UniqueToken<Cuda, UniqueTokenScope::Global>(arg) {}
#endif
};

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_CUDA
#endif  // KOKKOS_CUDA_UNIQUE_TOKEN_HPP
