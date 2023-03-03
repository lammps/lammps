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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_HPX
#include <Kokkos_HPX.hpp>

#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <hpx/local/condition_variable.hpp>
#include <hpx/local/init.hpp>
#include <hpx/local/thread.hpp>
#include <hpx/local/mutex.hpp>

#include <atomic>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

namespace Kokkos {
namespace Experimental {

bool HPX::m_hpx_initialized = false;
#if defined(KOKKOS_ENABLE_HPX_ASYNC_DISPATCH)
std::atomic<uint32_t> HPX::m_next_instance_id{HPX::impl_default_instance_id() +
                                              1};
uint32_t HPX::m_active_parallel_region_count{0};
hpx::spinlock HPX::m_active_parallel_region_count_mutex;
hpx::condition_variable_any HPX::m_active_parallel_region_count_cond;
HPX::instance_data HPX::m_default_instance_data;
#else
Kokkos::Impl::thread_buffer HPX::m_default_buffer;
#endif

int HPX::concurrency() {
  hpx::runtime *rt = hpx::get_runtime_ptr();
  if (rt == nullptr) {
    return hpx::threads::hardware_concurrency();
  } else {
    if (hpx::threads::get_self_ptr() == nullptr) {
      return hpx::resource::get_thread_pool(0).get_os_thread_count();
    } else {
      return hpx::this_thread::get_pool()->get_os_thread_count();
    }
  }
}

void HPX::impl_initialize(InitializationSettings const &settings) {
  hpx::runtime *rt = hpx::get_runtime_ptr();
  if (rt == nullptr) {
    hpx::local::init_params i;
    i.cfg = {
#ifdef KOKKOS_ENABLE_DEBUG
        "--hpx:attach-debugger=exception",
#endif
    };
    if (settings.has_num_threads()) {
      i.cfg.emplace_back("hpx.os_threads=" +
                         std::to_string(settings.get_num_threads()));
    }
    int argc_hpx     = 1;
    char name[]      = "kokkos_hpx";
    char *argv_hpx[] = {name, nullptr};
    hpx::local::start(nullptr, argc_hpx, argv_hpx, i);

    m_hpx_initialized = true;
  }
}

bool HPX::impl_is_initialized() noexcept {
  hpx::runtime *rt = hpx::get_runtime_ptr();
  return rt != nullptr;
}

void HPX::impl_finalize() {
  if (m_hpx_initialized) {
    hpx::runtime *rt = hpx::get_runtime_ptr();
    if (rt != nullptr) {
      hpx::apply([]() { hpx::local::finalize(); });
      hpx::local::stop();
    } else {
      Kokkos::abort(
          "Kokkos::Experimental::HPX::impl_finalize: Kokkos started "
          "HPX but something else already stopped HPX\n");
    }
  }
}

}  // namespace Experimental

namespace Impl {

int g_hpx_space_factory_initialized =
    initialize_space_factory<Kokkos::Experimental::HPX>("060_HPX");

}  // namespace Impl

}  // namespace Kokkos

#else
void KOKKOS_CORE_SRC_IMPL_HPX_PREVENT_LINK_ERROR() {}
#endif  //#ifdef KOKKOS_ENABLE_HPX
