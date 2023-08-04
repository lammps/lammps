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
#include <HPX/Kokkos_HPX.hpp>

#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <hpx/local/condition_variable.hpp>
#include <hpx/local/init.hpp>
#include <hpx/local/runtime.hpp>
#include <hpx/local/thread.hpp>
#include <hpx/local/mutex.hpp>
#include <hpx/version.hpp>

#include <atomic>
#include <chrono>
#include <memory>
#include <ostream>
#include <string>
#include <type_traits>

namespace Kokkos {
namespace Impl {
void hpx_thread_buffer::resize(const std::size_t num_threads,
                               const std::size_t size_per_thread,
                               const std::size_t extra_space) noexcept {
  m_num_threads     = num_threads;
  m_size_per_thread = size_per_thread;
  m_extra_space     = extra_space;

  pad_to_cache_line(m_size_per_thread);

  std::size_t size_total_new =
      m_num_threads * m_size_per_thread + m_extra_space;

  if (m_size_total < size_total_new) {
    // Don't use make_unique here as it value-initializes the elements of the
    // array, which we have no use for, and can be very slow for large arrays.
    m_data       = std::unique_ptr<char[]>(new char[size_total_new]);
    m_size_total = size_total_new;
  }
}

void *hpx_thread_buffer::get(std::size_t thread_num) const noexcept {
  KOKKOS_EXPECTS(thread_num < m_num_threads);
  if (!m_data) {
    return nullptr;
  }
  return &m_data[thread_num * m_size_per_thread];
}

void *hpx_thread_buffer::get_extra_space() const noexcept {
  KOKKOS_EXPECTS(m_extra_space > 0);
  if (!m_data) {
    return nullptr;
  }
  return &m_data[m_num_threads * m_size_per_thread];
}
}  // namespace Impl

namespace Experimental {

bool HPX::m_hpx_initialized = false;
std::atomic<uint32_t> HPX::m_next_instance_id{HPX::impl_default_instance_id() +
                                              1};
uint32_t HPX::m_active_parallel_region_count{0};
hpx::spinlock HPX::m_active_parallel_region_count_mutex;
hpx::condition_variable_any HPX::m_active_parallel_region_count_cond;
HPX::instance_data HPX::m_default_instance_data;

void HPX::print_configuration(std::ostream &os, const bool) const {
  os << "Host Parallel Execution Space\n";
  os << "  KOKKOS_ENABLE_HPX: yes\n";
  os << "HPX Options:\n";
#if defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH)
  os << "  KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH: yes\n";
#else
  os << "  KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH: no\n";
#endif
  os << "\nHPX Runtime Configuration:\n";
  os << "Worker threads: " << hpx::get_num_worker_threads() << '\n';
  os << hpx::complete_version() << '\n';
  os << hpx::configuration_string() << '\n';
}

bool &HPX::impl_get_in_parallel() noexcept {
  static thread_local bool in_parallel = false;
  return in_parallel;
}

HPX::impl_in_parallel_scope::impl_in_parallel_scope() noexcept {
  KOKKOS_EXPECTS(!impl_get_in_parallel());
  impl_get_in_parallel() = true;
}

HPX::impl_in_parallel_scope::~impl_in_parallel_scope() noexcept {
  KOKKOS_EXPECTS(impl_get_in_parallel());
  impl_get_in_parallel() = false;
}

HPX::impl_not_in_parallel_scope::impl_not_in_parallel_scope() noexcept {
  KOKKOS_EXPECTS(impl_get_in_parallel());
  impl_get_in_parallel() = false;
}

HPX::impl_not_in_parallel_scope::~impl_not_in_parallel_scope() noexcept {
  KOKKOS_EXPECTS(!impl_get_in_parallel());
  impl_get_in_parallel() = true;
}

void HPX::impl_decrement_active_parallel_region_count() {
  std::unique_lock<hpx::spinlock> l(m_active_parallel_region_count_mutex);
  if (--m_active_parallel_region_count == 0) {
    l.unlock();
    m_active_parallel_region_count_cond.notify_all();
  };
}

void HPX::impl_increment_active_parallel_region_count() {
  std::unique_lock<hpx::spinlock> l(m_active_parallel_region_count_mutex);
  ++m_active_parallel_region_count;
}

void HPX::impl_instance_fence_locked(const std::string &name) const {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HPX>(
      name,
      Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{
          impl_instance_id()},
      [&]() {
        auto &s = impl_get_sender();

        hpx::this_thread::experimental::sync_wait(std::move(s));
        s = hpx::execution::experimental::unique_any_sender(
            hpx::execution::experimental::just());
      });
}

void HPX::impl_instance_fence(const std::string &name) const {
  std::lock_guard<hpx::spinlock> l(impl_get_sender_mutex());
  impl_instance_fence_locked(name);
}

void HPX::impl_static_fence(const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HPX>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      [&]() {
        auto &s = HPX().impl_get_sender();

        std::unique_lock<hpx::spinlock> l(HPX().impl_get_sender_mutex());

        // This is a loose fence. Any work scheduled before this will be waited
        // for, but work scheduled while waiting may also be waited for.
        {
          std::unique_lock<hpx::spinlock> l_count(
              m_active_parallel_region_count_mutex);
          m_active_parallel_region_count_cond.wait(
              l_count, [&]() { return m_active_parallel_region_count == 0; });
        }

        hpx::this_thread::experimental::sync_wait(std::move(s));
        s = hpx::execution::experimental::unique_any_sender(
            hpx::execution::experimental::just());
      });
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int HPX::concurrency() {
#else
int HPX::concurrency() const {
#endif
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
#if HPX_VERSION_FULL >= 0x010900
      hpx::post([]() { hpx::local::finalize(); });
#else
      hpx::apply([]() { hpx::local::finalize(); });
#endif
      hpx::local::stop();
    } else {
      Kokkos::abort(
          "Kokkos::Experimental::HPX::impl_finalize: Kokkos started "
          "HPX but something else already stopped HPX\n");
    }
  }
}

int HPX::impl_thread_pool_size() noexcept {
  hpx::runtime *rt = hpx::get_runtime_ptr();
  if (rt == nullptr) {
    return 0;
  } else {
    if (hpx::threads::get_self_ptr() == nullptr) {
      return hpx::resource::get_thread_pool(0).get_os_thread_count();
    } else {
      return hpx::this_thread::get_pool()->get_os_thread_count();
    }
  }
}

int HPX::impl_thread_pool_rank() noexcept {
  hpx::runtime *rt = hpx::get_runtime_ptr();
  if (rt == nullptr) {
    return 0;
  } else {
    if (hpx::threads::get_self_ptr() == nullptr) {
      return 0;
    } else {
      return hpx::this_thread::get_pool()->get_pool_index();
    }
  }
}

int HPX::impl_thread_pool_size(int depth) {
  if (depth == 0) {
    return impl_thread_pool_size();
  } else {
    return 1;
  }
}

template void HPX::impl_bulk_plain_erased<int>(
    bool, bool, std::function<void(int)> &&, int const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_plain_erased<unsigned int>(
    bool, bool, std::function<void(unsigned int)> &&, unsigned int const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_plain_erased<long>(
    bool, bool, std::function<void(long)> &&, long const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_plain_erased<std::size_t>(
    bool, bool, std::function<void(std::size_t)> &&, std::size_t const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_setup_finalize_erased<int>(
    bool, bool, std::function<void(int)> &&, std::function<void()> &&,
    std::function<void()> &&, int const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_setup_finalize_erased<unsigned int>(
    bool, bool, std::function<void(unsigned int)> &&, std::function<void()> &&,
    std::function<void()> &&, unsigned int const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_setup_finalize_erased<long>(
    bool, bool, std::function<void(long)> &&, std::function<void()> &&,
    std::function<void()> &&, long const,
    hpx::threads::thread_stacksize stacksize) const;

template void HPX::impl_bulk_setup_finalize_erased<std::size_t>(
    bool, bool, std::function<void(std::size_t)> &&, std::function<void()> &&,
    std::function<void()> &&, std::size_t const,
    hpx::threads::thread_stacksize stacksize) const;
}  // namespace Experimental

namespace Impl {
int g_hpx_space_factory_initialized =
    initialize_space_factory<Kokkos::Experimental::HPX>("060_HPX");

}  // namespace Impl

}  // namespace Kokkos

#else
void KOKKOS_CORE_SRC_IMPL_HPX_PREVENT_LINK_ERROR() {}
#endif  //#ifdef KOKKOS_ENABLE_HPX
