// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_PARALLELFOR_RANGE_HPP
#define KOKKOS_NEXTSILICON_PARALLELFOR_RANGE_HPP

#include <NextSilicon/Kokkos_NextSilicon_Instance.hpp>
#include <NextSilicon/Kokkos_NextSilicon.hpp>
#include <NextSilicon/Kokkos_NextSilicon_Intrinsics.hpp>
#include <NextSilicon/Kokkos_NextSilicon_ThreadSpaceGuard.hpp>

#include <Kokkos_Parallel.hpp>

#include <nextapi/parallelism.hpp>

#include <cstddef>
#include <mutex>

template <class Functor, class... Traits>
class Kokkos::Impl::ParallelFor<Functor, Kokkos::RangePolicy<Traits...>,
                                Kokkos::Experimental::NextSilicon> {
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;
  using IndexType = typename Policy::index_type;

  Functor m_functor;
  Policy m_policy;

 public:
  ParallelFor(Functor const& functor, Policy const& policy)
      : m_functor(functor), m_policy(policy) {}

  void execute() const {
    // Acquire the device for potential handoff before kernel execution begins
    const std::lock_guard<std::mutex> device_lock =
        this->m_policy.space().impl_internal_space_instance()->lock_device();

    // Clone the driver to prevent the stack from getting migrated to device.
    auto internal_instance = m_policy.space().impl_internal_space_instance();
    auto cloned_driver     = internal_instance->clone_driver(*this);
    cloned_driver->execute_internal();
  }

 private:
  __attribute__((noinline)) void execute_internal() const {
    const IndexType begin = m_policy.begin();
    const IndexType end   = m_policy.end();

    if (end <= begin) return;

    const IndexType chunk_size = m_policy.chunk_size();
    nextapi::parallel_for(begin, end, {.chunk_size = chunk_size},
                          parallel_function, &m_functor);
  }

 private:
  /// Executes a single parallel iteration `index` for `functor`.
  static void parallel_function(IndexType index,
                                Functor const* __restrict functor) {
    // Set the thread's semantic execution space to device
    //   See Kokkos_NextSilicon_ThreadSpaceGuard.hpp for a detailed description
    //   of the mechanism and rationale.
    NextSiliconThreadSpaceGuard thread_guard{};

    // Communicate to the compiler that the functor is a immutable and thread
    // invariant for the duration of the microtask.
    Kokkos::Experimental::Impl::
        __ns_immutable_thread_invariant_parameter_struct(functor);

    // Invokes the functor itself. Expected to inline (if compiler visible) the
    // functor body while passing the extra this pointer.
    (*functor)(index);
  }
};

#endif  // KOKKOS_NEXTSILICON_PARALLELFOR_RANGE_HPP
