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

#ifndef KOKKOS_HIP_UNIQUE_TOKEN_HPP
#define KOKKOS_HIP_UNIQUE_TOKEN_HPP

#include <Kokkos_HIP_Space.hpp>
#include <Kokkos_UniqueToken.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

namespace Kokkos {

namespace Impl {
Kokkos::View<uint32_t*, Kokkos::Experimental::HIPSpace>
hip_global_unique_token_locks(bool deallocate = false);
}

namespace Experimental {

// both global and instance Unique Tokens are implemented in the same way
// the global version has one shared static lock array underneath
// but it can't be a static member variable since we need to acces it on device
// and we share the implementation with the instance version
template <>
class UniqueToken<HIP, UniqueTokenScope::Global> {
 protected:
  View<uint32_t*, HIPSpace> m_locks;

 public:
  using execution_space = HIP;
  using size_type       = int32_t;

  explicit UniqueToken(execution_space const& = HIP())
      : m_locks(Kokkos::Impl::hip_global_unique_token_locks()) {}

 protected:
  // These are constructors for the Instance version
  UniqueToken(size_type max_size) {
    m_locks = Kokkos::View<uint32_t*, HIPSpace>("Kokkos::UniqueToken::m_locks",
                                                max_size);
  }
  UniqueToken(size_type max_size, execution_space const& exec) {
    m_locks = Kokkos::View<uint32_t*, HIPSpace>(
        Kokkos::view_alloc(exec, "Kokkos::UniqueToken::m_locks"), max_size);
  }

 public:
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
  size_type size() const noexcept { return m_locks.extent(0); }

 private:
  __device__ size_type impl_acquire() const {
    int idx = blockIdx.x * (blockDim.x * blockDim.y) +
              threadIdx.y * blockDim.x + threadIdx.x;
    idx                            = idx % size();
    unsigned long long active      = __ballot(1);
    unsigned long long done_active = 0;
    bool done                      = false;
    while (active != done_active) {
      if (!done) {
        // Using m_locks(idx) fails self containment test of Kokkos_HIP.hpp
        // That failure stems from the space access verification because the
        // Host execution space is not defined
        if (Kokkos::atomic_compare_exchange(m_locks.data() + idx, 0, 1) == 0) {
          done = true;
        } else {
          idx += blockDim.y * blockDim.x + 1;
          idx = idx % size();
        }
      }
      done_active = __ballot(done ? 1 : 0);
    }

// Make sure that all writes in the previous lock owner are visible to me
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
    desul::atomic_thread_fence(desul::MemoryOrderAcquire(),
                               desul::MemoryScopeDevice());
#else
    Kokkos::memory_fence();
#endif
    return idx;
  }

 public:
  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type acquire() const {
    KOKKOS_IF_ON_DEVICE(return impl_acquire();)
    KOKKOS_IF_ON_HOST(return 0;)
  }

  /// \brief release an acquired value
  KOKKOS_INLINE_FUNCTION
  void release(size_type idx) const noexcept {
// Make sure my writes are visible to the next lock owner
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
    desul::atomic_thread_fence(desul::MemoryOrderRelease(),
                               desul::MemoryScopeDevice());
#else
    Kokkos::memory_fence();
#endif
    (void)Kokkos::atomic_exchange(m_locks.data() + idx, 0);
  }
};

template <>
class UniqueToken<HIP, UniqueTokenScope::Instance>
    : public UniqueToken<HIP, UniqueTokenScope::Global> {
 public:
  // The instance version will forward to protected constructor which creates
  // a lock array per instance
  UniqueToken()
      : UniqueToken<HIP, UniqueTokenScope::Global>(
            Kokkos::Experimental::HIP().concurrency()) {}
  explicit UniqueToken(execution_space const& arg)
      : UniqueToken<HIP, UniqueTokenScope::Global>(
            Kokkos::Experimental::HIP().concurrency(), arg) {}
  explicit UniqueToken(size_type max_size)
      : UniqueToken<HIP, UniqueTokenScope::Global>(max_size) {}
  UniqueToken(size_type max_size, execution_space const& arg)
      : UniqueToken<HIP, UniqueTokenScope::Global>(max_size, arg) {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_HIP_UNIQUE_TOKEN_HPP
