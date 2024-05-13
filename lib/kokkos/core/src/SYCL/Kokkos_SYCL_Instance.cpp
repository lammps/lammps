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

#include <Kokkos_Core.hpp>  //kokkos_malloc

#include <impl/Kokkos_CheckedIntegerOps.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

namespace {

// FIXME_SYCL Should be a multiple of the maximum subgroup size.
static constexpr auto sizeScratchGrain =
    sizeof(Kokkos::Experimental::SYCL::size_type[32]);

std::size_t scratch_count(const std::size_t size) {
  return (size + sizeScratchGrain - 1) / sizeScratchGrain;
}

}  // namespace

std::vector<std::optional<sycl::queue>*> SYCLInternal::all_queues;
std::mutex SYCLInternal::mutex;

Kokkos::View<uint32_t*, SYCLDeviceUSMSpace> sycl_global_unique_token_locks(
    bool deallocate) {
  static Kokkos::View<uint32_t*, SYCLDeviceUSMSpace> locks =
      Kokkos::View<uint32_t*, SYCLDeviceUSMSpace>();
  if (!deallocate && locks.extent(0) == 0)
    locks = Kokkos::View<uint32_t*, SYCLDeviceUSMSpace>(
        "Kokkos::UniqueToken<SYCL>::m_locks", SYCL().concurrency());
  if (deallocate) locks = Kokkos::View<uint32_t*, SYCLDeviceUSMSpace>();
  return locks;
}

SYCLInternal::~SYCLInternal() {
  if (!was_finalized || m_scratchSpace || m_scratchHost || m_scratchFlags) {
    std::cerr << "Kokkos::Experimental::SYCL ERROR: Failed to call "
                 "Kokkos::Experimental::SYCL::finalize()"
              << std::endl;
    std::cerr.flush();
  }
}

int SYCLInternal::verify_is_initialized(const char* const label) const {
  if (!is_initialized()) {
    Kokkos::abort((std::string("Kokkos::Experimental::SYCL::") + label +
                   " : ERROR device not initialized\n")
                      .c_str());
  }
  return is_initialized();
}
SYCLInternal& SYCLInternal::singleton() {
  static SYCLInternal self;
  return self;
}

void SYCLInternal::initialize(const sycl::device& d) {
  auto exception_handler = [](sycl::exception_list exceptions) {
    bool asynchronous_error = false;
    for (std::exception_ptr const& e : exceptions) {
      try {
        std::rethrow_exception(e);
      } catch (sycl::exception const& e) {
        std::cerr << e.what() << '\n';
        asynchronous_error = true;
      }
    }
    if (asynchronous_error)
      Kokkos::Impl::throw_runtime_exception(
          "There was an asynchronous SYCL error!\n");
  };
#ifdef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
  initialize(
      sycl::queue{d, exception_handler, sycl::property::queue::in_order()});
#else
  initialize(sycl::queue{d, exception_handler});
#endif
}

// FIXME_SYCL
void SYCLInternal::initialize(const sycl::queue& q) {
  KOKKOS_EXPECTS(!is_initialized());

#define KOKKOS_IMPL_CHECK_SYCL_BACKEND_SUPPORT(BACKEND, REQUIRED)            \
  if (BACKEND != REQUIRED)                                                   \
  Kokkos::abort(                                                             \
      "The SYCL execution space instance was initialized with an "           \
      "unsupported backend type! For this GPU architecture, only " #REQUIRED \
      " is supported.")
#if defined(KOKKOS_ARCH_INTEL_GPU)
  KOKKOS_IMPL_CHECK_SYCL_BACKEND_SUPPORT(q.get_backend(),
                                         sycl::backend::ext_oneapi_level_zero);
#elif defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  KOKKOS_IMPL_CHECK_SYCL_BACKEND_SUPPORT(q.get_backend(),
                                         sycl::backend::ext_oneapi_cuda);
#elif defined(KOKKOS_ARCH_AMD_GPU)
  KOKKOS_IMPL_CHECK_SYCL_BACKEND_SUPPORT(q.get_backend(),
                                         sycl::backend::ext_oneapi_hip);
#endif

  if (was_finalized)
    Kokkos::abort("Calling SYCL::initialize after SYCL::finalize is illegal\n");

  m_queue = q;
  // guard pushing to all_queues
  {
    std::scoped_lock lock(mutex);
    all_queues.push_back(&m_queue);
  }
  const sycl::device& d = m_queue->get_device();

  m_maxWorkgroupSize =
      d.template get_info<sycl::info::device::max_work_group_size>();
  // FIXME_SYCL this should give the correct value for NVIDIA GPUs
  m_maxConcurrency =
      m_maxWorkgroupSize * 2 *
      d.template get_info<sycl::info::device::max_compute_units>();

  m_maxShmemPerBlock =
      d.template get_info<sycl::info::device::local_mem_size>();

  for (auto& usm_mem : m_indirectKernelMem) {
    usm_mem.reset(*m_queue, m_instance_id);
  }

#ifdef KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED
  // Init the array for used for arbitrarily sized atomics
  if (this == &singleton()) {
    desul::Impl::init_lock_arrays();
    desul::Impl::init_lock_arrays_sycl(*m_queue);
  }
#endif
}

int SYCLInternal::acquire_team_scratch_space() {
  // Grab the next scratch memory allocation. We must make sure that the last
  // kernel using the allocation has completed, so we wait for the event that
  // was registered with that kernel.
  int current_team_scratch = desul::atomic_fetch_inc_mod(
      &m_current_team_scratch, m_n_team_scratch - 1,
      desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice());

  m_team_scratch_event[current_team_scratch].wait_and_throw();

  return current_team_scratch;
}

sycl::device_ptr<void> SYCLInternal::resize_team_scratch_space(
    int scratch_pool_id, std::int64_t bytes, bool force_shrink) {
  // Multiple ParallelFor/Reduce Teams can call this function at the same time
  // and invalidate the m_team_scratch_ptr. We use a pool to avoid any race
  // condition.
  if (m_team_scratch_current_size[scratch_pool_id] == 0) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_malloc<Experimental::SYCLDeviceUSMSpace>(
            "Kokkos::Experimental::SYCLDeviceUSMSpace::TeamScratchMemory",
            m_team_scratch_current_size[scratch_pool_id]);
  }
  if ((bytes > m_team_scratch_current_size[scratch_pool_id]) ||
      ((bytes < m_team_scratch_current_size[scratch_pool_id]) &&
       (force_shrink))) {
    m_team_scratch_current_size[scratch_pool_id] = bytes;
    m_team_scratch_ptr[scratch_pool_id] =
        Kokkos::kokkos_realloc<Experimental::SYCLDeviceUSMSpace>(
            m_team_scratch_ptr[scratch_pool_id],
            m_team_scratch_current_size[scratch_pool_id]);
  }
  return m_team_scratch_ptr[scratch_pool_id];
}

void SYCLInternal::register_team_scratch_event(int scratch_pool_id,
                                               sycl::event event) {
  m_team_scratch_event[scratch_pool_id] = event;
}

uint32_t SYCLInternal::impl_get_instance_id() const { return m_instance_id; }

void SYCLInternal::finalize() {
  SYCLInternal::fence(*m_queue,
                      "Kokkos::SYCLInternal::finalize: fence on finalization",
                      m_instance_id);
  was_finalized = true;

  // The global_unique_token_locks array is static and should only be
  // deallocated once by the defualt instance
  if (this == &singleton()) {
    Impl::sycl_global_unique_token_locks(true);
#ifdef KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED
    desul::Impl::finalize_lock_arrays();
    desul::Impl::finalize_lock_arrays_sycl(*m_queue);
#endif
  }

  auto device_mem_space = SYCLDeviceUSMSpace(*m_queue);
  auto host_mem_space   = SYCLHostUSMSpace(*m_queue);
  if (nullptr != m_scratchSpace)
    device_mem_space.deallocate(m_scratchSpace,
                                m_scratchSpaceCount * sizeScratchGrain);
  if (nullptr != m_scratchHost)
    host_mem_space.deallocate(m_scratchHost,
                              m_scratchHostCount * sizeScratchGrain);
  if (nullptr != m_scratchFlags)
    device_mem_space.deallocate(m_scratchFlags,
                                m_scratchFlagsCount * sizeScratchGrain);
  m_syclDev           = -1;
  m_scratchSpaceCount = 0;
  m_scratchSpace      = nullptr;
  m_scratchHostCount  = 0;
  m_scratchHost       = nullptr;
  m_scratchFlagsCount = 0;
  m_scratchFlags      = nullptr;

  for (int i = 0; i < m_n_team_scratch; ++i) {
    if (m_team_scratch_current_size[i] > 0) {
      Kokkos::kokkos_free<Kokkos::Experimental::SYCLDeviceUSMSpace>(
          m_team_scratch_ptr[i]);
      m_team_scratch_current_size[i] = 0;
      m_team_scratch_ptr[i]          = nullptr;
    }
  }

  for (auto& usm_mem : m_indirectKernelMem) usm_mem.reset();
  // guard erasing from all_queues
  {
    std::scoped_lock lock(mutex);
    all_queues.erase(std::find(all_queues.begin(), all_queues.end(), &m_queue));
  }
  m_queue.reset();
}

sycl::device_ptr<void> SYCLInternal::scratch_space(const std::size_t size) {
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount < scratch_count(size)) {
    auto mem_space = Kokkos::Experimental::SYCLDeviceUSMSpace(*m_queue);

    if (nullptr != m_scratchSpace)
      mem_space.deallocate(m_scratchSpace,
                           m_scratchSpaceCount * sizeScratchGrain);

    m_scratchSpaceCount = scratch_count(size);

    std::size_t alloc_size = Kokkos::Impl::multiply_overflow_abort(
        m_scratchSpaceCount, sizeScratchGrain);
    m_scratchSpace = static_cast<size_type*>(mem_space.allocate(
        "Kokkos::Experimental::SYCL::InternalScratchSpace", alloc_size));
  }

  return m_scratchSpace;
}

sycl::host_ptr<void> SYCLInternal::scratch_host(const std::size_t size) {
  if (verify_is_initialized("scratch_unified") &&
      m_scratchHostCount < scratch_count(size)) {
    auto mem_space = Kokkos::Experimental::SYCLHostUSMSpace(*m_queue);

    if (nullptr != m_scratchHost)
      mem_space.deallocate(m_scratchHost,
                           m_scratchHostCount * sizeScratchGrain);

    m_scratchHostCount = scratch_count(size);

    std::size_t alloc_size = Kokkos::Impl::multiply_overflow_abort(
        m_scratchHostCount, sizeScratchGrain);
    m_scratchHost = static_cast<size_type*>(mem_space.allocate(
        "Kokkos::Experimental::SYCL::InternalScratchHost", alloc_size));
  }

  return m_scratchHost;
}

sycl::device_ptr<void> SYCLInternal::scratch_flags(const std::size_t size) {
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount < scratch_count(size)) {
    auto mem_space = Kokkos::Experimental::SYCLDeviceUSMSpace(*m_queue);

    if (nullptr != m_scratchFlags)
      mem_space.deallocate(m_scratchFlags,
                           m_scratchFlagsCount * sizeScratchGrain);

    m_scratchFlagsCount = scratch_count(size);

    std::size_t alloc_size = Kokkos::Impl::multiply_overflow_abort(
        m_scratchFlagsCount, sizeScratchGrain);
    m_scratchFlags = static_cast<size_type*>(mem_space.allocate(
        "Kokkos::Experimental::SYCL::InternalScratchFlags", alloc_size));

    // We only zero-initialize the allocation when we actually allocate.
    // It's the responsibility of the features using scratch_flags,
    // namely parallel_reduce and parallel_scan, to reset the used values to 0.
    auto memset_event = m_queue->memset(m_scratchFlags, 0,
                                        m_scratchFlagsCount * sizeScratchGrain);
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
    m_queue->ext_oneapi_submit_barrier(std::vector{memset_event});
#endif
  }

  return m_scratchFlags;
}

template <typename WAT>
void SYCLInternal::fence_helper(WAT& wat, const std::string& name,
                                uint32_t instance_id) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::SYCL>(
      name, Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{instance_id},
      [&]() {
        try {
          wat.wait_and_throw();
        } catch (sycl::exception const& e) {
          Kokkos::Impl::throw_runtime_exception(
              std::string("There was a synchronous SYCL error:\n") += e.what());
        }
      });
}
template void SYCLInternal::fence_helper<sycl::queue>(sycl::queue&,
                                                      const std::string&,
                                                      uint32_t);
template void SYCLInternal::fence_helper<sycl::event>(sycl::event&,
                                                      const std::string&,
                                                      uint32_t);

// This function cycles through a pool of USM allocations for functors
SYCLInternal::IndirectKernelMem& SYCLInternal::get_indirect_kernel_mem() {
  // Thread safety: atomically increment round robin variable
  // NB: atomic_fetch_inc_mod returns values in range [0-N], not
  // [0-N) as might be expected.
  size_t next_pool = desul::atomic_fetch_inc_mod(
      &m_pool_next, m_usm_pool_size - 1, desul::MemoryOrderRelaxed(),
      desul::MemoryScopeDevice());
  return m_indirectKernelMem[next_pool];
}

template <sycl::usm::alloc Kind>
size_t SYCLInternal::USMObjectMem<Kind>::reserve(size_t n) {
  assert(m_q);

  if (m_capacity < n) {
    AllocationSpace alloc_space(*m_q);
    if (m_data) alloc_space.deallocate(m_data, m_capacity);

    m_data =
        alloc_space.allocate("Kokkos::Experimental::SYCL::USMObjectMem", n);

    if constexpr (sycl::usm::alloc::device == Kind)
      m_staging.reset(new char[n]);
    m_capacity = n;
  }

  return m_capacity;
}

template <sycl::usm::alloc Kind>
void SYCLInternal::USMObjectMem<Kind>::reset() {
  if (m_data) {
    // This implies a fence since this class is not copyable
    // and deallocating implies a fence across all registered queues.
    AllocationSpace alloc_space(*m_q);
    alloc_space.deallocate(m_data, m_capacity);

    m_capacity = 0;
    m_data     = nullptr;
  }
  m_q.reset();
}

int SYCLInternal::m_syclDev;

template class SYCLInternal::USMObjectMem<sycl::usm::alloc::shared>;
template class SYCLInternal::USMObjectMem<sycl::usm::alloc::device>;
template class SYCLInternal::USMObjectMem<sycl::usm::alloc::host>;

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos
