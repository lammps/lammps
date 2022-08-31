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

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_HIP.hpp>
#include <Kokkos_HIP_Space.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <atomic>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {

namespace {
hipStream_t get_deep_copy_stream() {
  static hipStream_t s = nullptr;
  if (s == nullptr) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&s));
  }
  return s;
}
}  // namespace

void DeepCopyHIP(void* dst, void const* src, size_t n) {
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyAsync(dst, src, n, hipMemcpyDefault));
}

void DeepCopyAsyncHIP(const Kokkos::Experimental::HIP& instance, void* dst,
                      void const* src, size_t n) {
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

void DeepCopyAsyncHIP(void* dst, void const* src, size_t n) {
  hipStream_t s = get_deep_copy_stream();
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyAsync(dst, src, n, hipMemcpyDefault, s));
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HIP>(
      "Kokkos::Impl::DeepCopyAsyncHIP: Post Deep Copy Fence on Deep-Copy "
      "stream",
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          DeepCopyResourceSynchronization,
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(s)); });
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
namespace Kokkos {

KOKKOS_DEPRECATED void Experimental::HIPSpace::access_error() {
  const std::string msg(
      "Kokkos::Experimental::HIPSpace::access_error attempt to execute "
      "Experimental::HIP function from non-HIP space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

KOKKOS_DEPRECATED void Experimental::HIPSpace::access_error(const void* const) {
  const std::string msg(
      "Kokkos::Experimental::HIPSpace::access_error attempt to execute "
      "Experimental::HIP function from non-HIP space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

}  // namespace Kokkos
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

HIPSpace::HIPSpace() : m_device(HIP().hip_device()) {}

HIPHostPinnedSpace::HIPHostPinnedSpace() {}

void* HIPSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPSpace::allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPSpace::impl_allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  auto const error_code = hipMalloc(&ptr, arg_alloc_size);
  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    throw HIPRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        RawMemoryAllocationFailure::AllocationMechanism::HIPMalloc);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void* HIPHostPinnedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPHostPinnedSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPHostPinnedSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  auto const error_code = hipHostMalloc(&ptr, arg_alloc_size);
  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    throw HIPRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        RawMemoryAllocationFailure::AllocationMechanism::HIPHostMalloc);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}
void HIPSpace::deallocate(void* const arg_alloc_ptr,
                          const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void HIPSpace::deallocate(const char* arg_label, void* const arg_alloc_ptr,
                          const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
}

void HIPHostPinnedSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HIPHostPinnedSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPHostPinnedSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(arg_alloc_ptr));
}

}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::s_root_record;

SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::HIPHostPinnedSpace, void>::s_root_record;
#endif

SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                       void>::~SharedAllocationRecord() {
  auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, (alloc_size - sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                       void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::HIPSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {

  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // Copy to device memory
  Kokkos::Experimental::HIP exec;
  Kokkos::Impl::DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace>(
      exec, RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
  exec.fence(
      "SharedAllocationRecord<Kokkos::Experimental::HIPSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  // Fill in the Header information, directly accessible via host pinned memory
  this->base_t::_fill_host_accessible_header_info(*RecordBase::m_alloc_ptr,
                                                  arg_label);
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {
int get_gpu(const InitArguments& args);
}
namespace Experimental {

int HIP::concurrency() {
  auto const& prop = hip_device_prop();
  return prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount;
}
int HIP::impl_is_initialized() {
  return Impl::HIPInternal::singleton().is_initialized();
}

void HIP::impl_initialize(const HIP::SelectDevice config) {
  Impl::HIPInternal::singleton().initialize(config.hip_device_id);
}

void HIP::impl_finalize() { Impl::HIPInternal::singleton().finalize(); }

HIP::HIP()
    : m_space_instance(&Impl::HIPInternal::singleton(),
                       [](Impl::HIPInternal*) {}) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
}

HIP::HIP(hipStream_t const stream, bool manage_stream)
    : m_space_instance(new Impl::HIPInternal, [](Impl::HIPInternal* ptr) {
        ptr->finalize();
        delete ptr;
      }) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
  m_space_instance->initialize(Impl::HIPInternal::singleton().m_hipDev, stream,
                               manage_stream);
}

void HIP::print_configuration(std::ostream& s, const bool) {
  Impl::HIPInternal::singleton().print_configuration(s);
}

uint32_t HIP::impl_instance_id() const noexcept {
  return m_space_instance->impl_get_instance_id();
}
void HIP::impl_static_fence(const std::string& name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::HIP>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize()); });
}
void HIP::impl_static_fence() {
  impl_static_fence("Kokkos::HIP::impl_static_fence: Unnamed Static Fence");
}

void HIP::fence(const std::string& name) const {
  m_space_instance->fence(name);
}
void HIP::fence() const {
  fence("Kokkos::HIP::fence(): Unnamed Instance Fence");
}

hipStream_t HIP::hip_stream() const { return m_space_instance->m_stream; }

int HIP::hip_device() const { return impl_internal_space_instance()->m_hipDev; }

hipDeviceProp_t const& HIP::hip_device_prop() {
  return Impl::HIPInternal::singleton().m_deviceProp;
}

const char* HIP::name() { return "HIP"; }

}  // namespace Experimental

namespace Impl {

int g_hip_space_factory_initialized =
    initialize_space_factory<HIPSpaceInitializer>("150_HIP");

void HIPSpaceInitializer::initialize(const InitArguments& args) {
  int use_gpu = Impl::get_gpu(args);

  if (std::is_same<Kokkos::Experimental::HIP,
                   Kokkos::DefaultExecutionSpace>::value ||
      0 < use_gpu) {
    if (use_gpu > -1) {
      Kokkos::Experimental::HIP::impl_initialize(
          Kokkos::Experimental::HIP::SelectDevice(use_gpu));
    } else {
      Kokkos::Experimental::HIP::impl_initialize();
    }
  }
}

void HIPSpaceInitializer::finalize(const bool all_spaces) {
  if (std::is_same<Kokkos::Experimental::HIP,
                   Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
    if (Kokkos::Experimental::HIP::impl_is_initialized())
      Kokkos::Experimental::HIP::impl_finalize();
  }
}

void HIPSpaceInitializer::fence() {
  Kokkos::Experimental::HIP::impl_static_fence();
}
void HIPSpaceInitializer::fence(const std::string& name) {
  Kokkos::Experimental::HIP::impl_static_fence(name);
}

void HIPSpaceInitializer::print_configuration(std::ostream& msg,
                                              const bool detail) {
  msg << "Devices:" << std::endl;
  msg << "  KOKKOS_ENABLE_HIP: ";
  msg << "yes" << std::endl;

  msg << "HIP Options:" << std::endl;
  msg << "  KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE: ";
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "\nRuntime Configuration:" << std::endl;
  Experimental::HIP::print_configuration(msg, detail);
}

}  // namespace Impl

#ifdef KOKKOS_ENABLE_CXX14
namespace Tools {
namespace Experimental {
constexpr DeviceType DeviceTypeTraits<Kokkos::Experimental::HIP>::id;
}
}  // namespace Tools
#endif

}  // namespace Kokkos

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

namespace Kokkos {
namespace Impl {

// To avoid additional compilation cost for something that's (mostly?) not
// performance sensitive, we explicity instantiate these CRTP base classes here,
// where we have access to the associated *_timpl.hpp header files.
template class HostInaccessibleSharedAllocationRecordCommon<
    Kokkos::Experimental::HIPSpace>;
template class SharedAllocationRecordCommon<Kokkos::Experimental::HIPSpace>;
template class SharedAllocationRecordCommon<
    Kokkos::Experimental::HIPHostPinnedSpace>;

}  // end namespace Impl
}  // end namespace Kokkos

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
