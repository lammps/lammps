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

#ifndef KOKKOS_SYCL_UNIQUE_TOKEN_HPP
#define KOKKOS_SYCL_UNIQUE_TOKEN_HPP

#include <impl/Kokkos_ConcurrentBitset.hpp>
#include <SYCL/Kokkos_SYCL_Space.hpp>
#include <Kokkos_UniqueToken.hpp>

namespace Kokkos {
namespace Experimental {

namespace Impl {
Kokkos::View<uint32_t*, SYCLDeviceUSMSpace> sycl_global_unique_token_locks(
    bool deallocate = false);
}

// both global and instance Unique Tokens are implemented in the same way
// the global version has one shared static lock array underneath
// but it can't be a static member variable since we need to acces it on device
// and we share the implementation with the instance version
template <>
class UniqueToken<SYCL, UniqueTokenScope::Global> {
  Kokkos::View<uint32_t*, SYCLDeviceUSMSpace> m_locks;

 public:
  using execution_space = SYCL;
  using size_type       = int32_t;

  explicit UniqueToken(execution_space const& = execution_space())
      : m_locks(Impl::sycl_global_unique_token_locks()) {}

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

 protected:
  // Constructors for the Instance version
  UniqueToken(size_type max_size)
      : m_locks(Kokkos::View<uint32_t*, SYCLDeviceUSMSpace>(
            "Kokkos::UniqueToken::m_locks", max_size)) {}

  UniqueToken(size_type max_size, execution_space const& arg)
      : m_locks(Kokkos::View<uint32_t*, SYCLDeviceUSMSpace>(
            Kokkos::view_alloc(arg, "Kokkos::UniqueToken::m_locks"),
            max_size)) {}

 private:
  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  size_type impl_acquire() const {
    auto item = sycl::ext::oneapi::experimental::this_nd_item<3>();
    std::size_t threadIdx[3] = {item.get_local_id(2), item.get_local_id(1),
                                item.get_local_id(0)};
    std::size_t blockIdx[3]  = {item.get_group(2), item.get_group(1),
                               item.get_group(0)};
    std::size_t blockDim[3] = {item.get_local_range(2), item.get_local_range(1),
                               item.get_local_range(0)};

    int idx = blockIdx[0] * (blockDim[0] * blockDim[1]) +
              threadIdx[1] * blockDim[0] + threadIdx[0];
    idx %= size();

    while (Kokkos::atomic_compare_exchange(&m_locks(idx), 0, 1) == 1) {
      idx += blockDim[1] * blockDim[0] + 1;
      idx %= size();
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
    (void)Kokkos::atomic_exchange(&m_locks(idx), 0);
  }
};

template <>
class UniqueToken<SYCL, UniqueTokenScope::Instance>
    : public UniqueToken<SYCL, UniqueTokenScope::Global> {
 public:
  UniqueToken()
      : UniqueToken<SYCL, UniqueTokenScope::Global>(
            Kokkos::Experimental::SYCL().concurrency()) {}

  explicit UniqueToken(execution_space const& arg)
      : UniqueToken<SYCL, UniqueTokenScope::Global>(
            Kokkos::Experimental::SYCL().concurrency(), arg) {}

  explicit UniqueToken(size_type max_size)
      : UniqueToken<SYCL, UniqueTokenScope::Global>(max_size) {}

  UniqueToken(size_type max_size, execution_space const& arg)
      : UniqueToken<SYCL, UniqueTokenScope::Global>(max_size, arg) {}
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
