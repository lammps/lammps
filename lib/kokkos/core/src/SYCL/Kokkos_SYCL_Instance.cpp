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

namespace Kokkos {
namespace Experimental {
namespace Impl {

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
  if (!was_finalized || m_scratchSpace || m_scratchFlags) {
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
  // FIXME SYCL+Cuda Apparently, submit_barrier doesn't quite work as expected
  // for oneAPI 2023.0.0 on NVIDIA GPUs.
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  initialize(
      sycl::queue{d, exception_handler, sycl::property::queue::in_order()});
#else
  initialize(sycl::queue{d, exception_handler});
#endif
}

// FIXME_SYCL
void SYCLInternal::initialize(const sycl::queue& q) {
  if (was_finalized)
    Kokkos::abort("Calling SYCL::initialize after SYCL::finalize is illegal\n");

  if (is_initialized()) return;

  if (!HostSpace::execution_space::impl_is_initialized()) {
    const std::string msg(
        "SYCL::initialize ERROR : HostSpace::execution_space is not "
        "initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  const bool ok_init = nullptr == m_scratchSpace || nullptr == m_scratchFlags;
  const bool ok_dev  = true;
  if (ok_init && ok_dev) {
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

  } else {
    std::ostringstream msg;
    msg << "Kokkos::Experimental::SYCL::initialize(...) FAILED";

    if (!ok_init) {
      msg << " : Already initialized";
    }
    Kokkos::Impl::throw_runtime_exception(msg.str());
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

  using RecordSYCL = Kokkos::Impl::SharedAllocationRecord<SYCLDeviceUSMSpace>;
  if (nullptr != m_scratchSpace)
    RecordSYCL::decrement(RecordSYCL::get_record(m_scratchSpace));
  if (nullptr != m_scratchFlags)
    RecordSYCL::decrement(RecordSYCL::get_record(m_scratchFlags));
  m_syclDev           = -1;
  m_scratchSpaceCount = 0;
  m_scratchSpace      = nullptr;
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
  const size_type sizeScratchGrain =
      sizeof(Kokkos::Experimental::SYCL::size_type);
  if (verify_is_initialized("scratch_space") &&
      m_scratchSpaceCount * sizeScratchGrain < size) {
    m_scratchSpaceCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record = Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::SYCLDeviceUSMSpace, void>;

    if (nullptr != m_scratchSpace)
      Record::decrement(Record::get_record(m_scratchSpace));

    Record* const r =
        Record::allocate(Kokkos::Experimental::SYCLDeviceUSMSpace(*m_queue),
                         "Kokkos::Experimental::SYCL::InternalScratchSpace",
                         (sizeScratchGrain * m_scratchSpaceCount));

    Record::increment(r);

    m_scratchSpace = reinterpret_cast<size_type*>(r->data());
  }

  return m_scratchSpace;
}

sycl::device_ptr<void> SYCLInternal::scratch_flags(const std::size_t size) {
  const size_type sizeScratchGrain =
      sizeof(Kokkos::Experimental::SYCL::size_type);
  if (verify_is_initialized("scratch_flags") &&
      m_scratchFlagsCount * sizeScratchGrain < size) {
    m_scratchFlagsCount = (size + sizeScratchGrain - 1) / sizeScratchGrain;

    using Record = Kokkos::Impl::SharedAllocationRecord<
        Kokkos::Experimental::SYCLDeviceUSMSpace, void>;

    if (nullptr != m_scratchFlags)
      Record::decrement(Record::get_record(m_scratchFlags));

    Record* const r =
        Record::allocate(Kokkos::Experimental::SYCLDeviceUSMSpace(*m_queue),
                         "Kokkos::Experimental::SYCL::InternalScratchFlags",
                         (sizeScratchGrain * m_scratchFlagsCount));

    Record::increment(r);

    m_scratchFlags = reinterpret_cast<size_type*>(r->data());
  }
  auto memset_event = m_queue->memset(m_scratchFlags, 0,
                                      m_scratchFlagsCount * sizeScratchGrain);
  m_queue->ext_oneapi_submit_barrier(std::vector{memset_event});

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
    using Record = Kokkos::Impl::SharedAllocationRecord<AllocationSpace, void>;
    // First free what we have (in case malloc can reuse it)
    if (m_data) Record::decrement(Record::get_record(m_data));

    Record* const r = Record::allocate(
        AllocationSpace(*m_q), "Kokkos::Experimental::SYCL::USMObjectMem", n);
    Record::increment(r);

    m_data = r->data();
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
    using Record = Kokkos::Impl::SharedAllocationRecord<AllocationSpace, void>;
    Record::decrement(Record::get_record(m_data));

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
