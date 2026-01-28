// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)

// FIXME_OPENMPTARGET - macro for workaround implementation in UniqueToken
// constructor. undef'ed at the end
#define KOKKOS_IMPL_OPENMPTARGET_WORKAROUND

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <OpenMPTarget/Kokkos_OpenMPTarget.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_UniqueToken.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Instance.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <sstream>

namespace Kokkos {
namespace Experimental {
namespace Impl {
uint32_t OpenMPTargetInternal::impl_get_instance_id() const noexcept {
  return m_instance_id;
}

void OpenMPTargetInternal::fence(openmp_fence_is_static is_static) {
  fence(
      "Kokkos::Experimental::Impl::OpenMPTargetInternal::fence: Unnamed "
      "Internal Fence",
      is_static);
}
void OpenMPTargetInternal::fence(const std::string& name,
                                 openmp_fence_is_static is_static) {
  if (is_static == openmp_fence_is_static::no) {
    Kokkos::Tools::Experimental::Impl::profile_fence_event<
        Kokkos::Experimental::OpenMPTarget>(
        name,
        Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{
            impl_get_instance_id()},
        [&]() {});
  } else {
    Kokkos::Tools::Experimental::Impl::profile_fence_event<
        Kokkos::Experimental::OpenMPTarget>(
        name,
        Kokkos::Tools::Experimental::SpecialSynchronizationCases::
            GlobalDeviceSynchronization,
        [&]() {});
  }
}
int OpenMPTargetInternal::concurrency() const {
  int max_threads_sm = 2048;
  int max_threads    = max_threads_sm * 80;
#if defined(KOKKOS_ARCH_AMPERE86)
  max_threads = max_threads_sm * 84;
#elif defined(KOKKOS_ARCH_AMPERE87)
  max_threads_sm = 1024;
  max_threads    = max_threads_sm * 32;  // Orin Nano cores
#elif defined(KOKKOS_ARCH_AMPERE80)
  return max_threads_sm * 108;
#elif defined(KOKKOS_ARCH_VOLTA72)
  return max_threads_sm * 84;
#elif defined(KOKKOS_ARCH_VOLTA70)
  return max_threads_sm * 80;
#elif defined(KOKKOS_ARCH_PASCAL60) || defined(KOKKOS_ARCH_PASCAL61)
  return max_threads_sm * 60;
#endif

  return max_threads;
}
const char* OpenMPTargetInternal::name() { return "OpenMPTarget"; }
void OpenMPTargetInternal::print_configuration(std::ostream& os,
                                               bool /*verbose*/) const {
  // FIXME_OPENMPTARGET
  os << "Using OpenMPTarget\n";
}

void OpenMPTargetInternal::impl_finalize() {
  if (m_uniquetoken_ptr != nullptr)
    Kokkos::kokkos_free<Kokkos::Experimental::OpenMPTargetSpace>(
        m_uniquetoken_ptr);
}

void OpenMPTargetInternal::impl_initialize() {
  // FIXME_OPENMPTARGET:  Only fix the number of teams for NVIDIA architectures
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  omp_set_num_teams(512);
#endif
}

OpenMPTargetInternal* OpenMPTargetInternal::impl_singleton() {
  static OpenMPTargetInternal self;
  return &self;
}

void OpenMPTargetInternal::verify_is_process(const char* const label) {
  // Fails if the current task is in a parallel region or is not on the host.
  if (omp_in_parallel() && (!omp_is_initial_device())) {
    std::string msg(label);
    msg.append(" ERROR: in parallel or on device");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void OpenMPTargetInternal::clear_scratch() {
  Kokkos::Experimental::OpenMPTargetSpace space;
  space.deallocate(m_scratch_ptr, m_scratch_size);
  m_scratch_ptr  = nullptr;
  m_scratch_size = 0;
}

void* OpenMPTargetInternal::get_scratch_ptr() { return m_scratch_ptr; }

void OpenMPTargetInternal::resize_scratch(int64_t team_size,
                                          int64_t shmem_size_L0,
                                          int64_t shmem_size_L1,
                                          int64_t league_size) {
  Kokkos::Experimental::OpenMPTargetSpace space;
  // Level-0 scratch when using clang/17 and higher comes from their OpenMP
  // extension, `ompx_dyn_cgroup_mem`.
#if defined(KOKKOS_IMPL_OPENMPTARGET_LLVM_EXTENSIONS)
  shmem_size_L0 = 0;
#endif
  const int64_t shmem_size =
      shmem_size_L0 + shmem_size_L1;  // L0 + L1 scratch memory per team.
  const int64_t padding = shmem_size * 10 / 100;  // Padding per team.

  // Maximum active teams possible.
  // The number should not exceed the maximum in-flight teams possible or the
  // league_size.
  int max_active_teams =
      std::min(OpenMPTargetInternal::concurrency() / team_size, league_size);

  // max_active_teams is the number of active teams on the given hardware.
  // We set the number of teams to be twice the number of max_active_teams for
  // the compiler to pick the right number in its case.
  omp_set_num_teams(max_active_teams * 2);

  // Total amount of scratch memory allocated is depenedent
  // on the maximum number of in-flight teams possible.
  int64_t total_size =
      (shmem_size +
       ::Kokkos::Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE + padding) *
      max_active_teams * 2;

  if (total_size > m_scratch_size) {
    space.deallocate(m_scratch_ptr, m_scratch_size);
    m_scratch_size = total_size;
    m_scratch_ptr  = space.allocate(total_size);
  }
}

}  // namespace Impl

OpenMPTarget::OpenMPTarget()
    : m_space_instance(Impl::OpenMPTargetInternal::impl_singleton()) {}

const char* OpenMPTarget::name() {
  return Impl::OpenMPTargetInternal::impl_singleton()->name();
}
void OpenMPTarget::print_configuration(std::ostream& os, bool verbose) const {
  os << "OpenMPTarget Execution Space:\n";
  os << "  KOKKOS_ENABLE_OPENMPTARGET: yes\n";

  os << "\nOpenMPTarget Runtime Configuration:\n";

  m_space_instance->print_configuration(os, verbose);
}

uint32_t OpenMPTarget::impl_instance_id() const noexcept {
  return m_space_instance->impl_get_instance_id();
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int OpenMPTarget::concurrency() {
  return Impl::OpenMPTargetInternal::impl_singleton()->concurrency();
}
#else
int OpenMPTarget::concurrency() const {
  return m_space_instance->concurrency();
}
#endif

void OpenMPTarget::fence(const std::string& name) {
  Impl::OpenMPTargetInternal::impl_singleton()->fence(name);
}

void OpenMPTarget::impl_static_fence(const std::string& name) {
  Impl::OpenMPTargetInternal::impl_singleton()->fence(
      name, Kokkos::Experimental::Impl::openmp_fence_is_static::yes);
}

void OpenMPTarget::impl_initialize(InitializationSettings const& settings) {
  using Kokkos::Impl::get_visible_devices;
  std::vector<int> const& visible_devices = get_visible_devices();
  using Kokkos::Impl::get_gpu;
  const int device_num = get_gpu(settings).value_or(visible_devices[0]);
  omp_set_default_device(device_num);

  Impl::OpenMPTargetInternal::impl_singleton()->impl_initialize();
}
void OpenMPTarget::impl_finalize() {
  Impl::OpenMPTargetInternal::impl_singleton()->impl_finalize();
}
}  // Namespace Experimental

namespace Impl {
int g_openmptarget_space_factory_initialized =
    Kokkos::Impl::initialize_space_factory<Experimental::OpenMPTarget>(
        "160_OpenMPTarget");

}  // namespace Impl
}  // Namespace Kokkos

namespace Kokkos {
namespace Experimental {

UniqueToken<Kokkos::Experimental::OpenMPTarget,
            Kokkos::Experimental::UniqueTokenScope::Global>::
    UniqueToken(Kokkos::Experimental::OpenMPTarget const& space) {
#ifdef KOKKOS_IMPL_OPENMPTARGET_WORKAROUND
  uint32_t* ptr = space.impl_internal_space_instance()->m_uniquetoken_ptr;
  int count     = Kokkos::Experimental::OpenMPTarget().concurrency();
  if (ptr == nullptr) {
    int size = count * sizeof(uint32_t);
    ptr      = static_cast<uint32_t*>(
        Kokkos::kokkos_malloc<Kokkos::Experimental::OpenMPTargetSpace>(
            "Kokkos::OpenMPTarget::m_uniquetoken_ptr", size));
    std::vector<uint32_t> h_buf(count, 0);
    if (0 < size)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(ptr, h_buf.data(), size, 0,
                                                   0, omp_get_default_device(),
                                                   omp_get_initial_device()));

    space.impl_internal_space_instance()->m_uniquetoken_ptr = ptr;
  }
#else
// FIXME_OPENMPTARGET - 2 versions of non-working implementations to fill `ptr`
// with 0's
// Version 1 - Creating a target region and filling the
// pointer Error - CUDA error: named symbol not found
#pragma omp target teams distribute parallel for is_device_ptr(ptr) \
    map(to : size)
  for (int i = 0; i < count; ++i) ptr[i] = 0;

  // Version 2 : Allocating a view on the device and filling it with a scalar
  // value of 0.
  Kokkos::View<uint32_t*, Kokkos::Experimental::OpenMPTargetSpace> ptr_view(
      ptr, count);
  Kokkos::deep_copy(ptr_view, 0);
#endif
  m_buffer = ptr;
  m_count  = count;
}
}  // namespace Experimental
}  // namespace Kokkos

#undef KOKKOS_IMPL_OPENMPTARGET_WORKAROUND
#endif  // defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)
