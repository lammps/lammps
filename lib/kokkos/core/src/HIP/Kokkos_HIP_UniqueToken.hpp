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

#ifndef KOKKOS_HIP_UNIQUE_TOKEN_HPP
#define KOKKOS_HIP_UNIQUE_TOKEN_HPP

#include <HIP/Kokkos_HIP_Space.hpp>
#include <Kokkos_UniqueToken.hpp>

namespace Kokkos {

namespace Impl {
Kokkos::View<uint32_t*, HIPSpace> hip_global_unique_token_locks(
    bool deallocate = false);
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
    desul::atomic_thread_fence(desul::MemoryOrderAcquire(),
                               desul::MemoryScopeDevice());
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
    desul::atomic_thread_fence(desul::MemoryOrderRelease(),
                               desul::MemoryScopeDevice());
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
      : UniqueToken<HIP, UniqueTokenScope::Global>(HIP().concurrency()) {}
  explicit UniqueToken(execution_space const& arg)
      : UniqueToken<HIP, UniqueTokenScope::Global>(HIP().concurrency(), arg) {}
  explicit UniqueToken(size_type max_size)
      : UniqueToken<HIP, UniqueTokenScope::Global>(max_size) {}
  UniqueToken(size_type max_size, execution_space const& arg)
      : UniqueToken<HIP, UniqueTokenScope::Global>(max_size, arg) {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_HIP_UNIQUE_TOKEN_HPP
