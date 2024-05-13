/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_THREAD_FENCE_MSVC_HPP_
#define DESUL_ATOMICS_THREAD_FENCE_MSVC_HPP_

#include <atomic>
#include <desul/atomics/Adapt_CXX.hpp>
#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

template <class MemoryOrder, class MemoryScope>
void host_atomic_thread_fence(MemoryOrder, MemoryScope) {
  std::atomic_thread_fence(Impl::CXXMemoryOrder<MemoryOrder>::value);
}

}  // namespace Impl
}  // namespace desul

#endif
