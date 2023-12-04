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

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <Kokkos_Core.hpp>
#include <Cuda/Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_CudaSpace.hpp>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <atomic>

//#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>

#include <impl/Kokkos_Tools.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

cudaStream_t Kokkos::Impl::cuda_get_deep_copy_stream() {
  static cudaStream_t s = nullptr;
  if (s == nullptr) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (CudaInternal::singleton().cuda_stream_create_wrapper(&s)));
  }
  return s;
}

const std::unique_ptr<Kokkos::Cuda> &Kokkos::Impl::cuda_get_deep_copy_space(
    bool initialize) {
  static std::unique_ptr<Cuda> space = nullptr;
  if (!space && initialize)
    space = std::make_unique<Cuda>(Kokkos::Impl::cuda_get_deep_copy_stream());
  return space;
}

namespace Kokkos {
namespace Impl {

namespace {

static std::atomic<int> num_uvm_allocations(0);

}  // namespace

void DeepCopyCuda(void *dst, const void *src, size_t n) {
  KOKKOS_IMPL_CUDA_SAFE_CALL((CudaInternal::singleton().cuda_memcpy_wrapper(
      dst, src, n, cudaMemcpyDefault)));
}

void DeepCopyAsyncCuda(const Cuda &instance, void *dst, const void *src,
                       size_t n) {
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (instance.impl_internal_space_instance()->cuda_memcpy_async_wrapper(
          dst, src, n, cudaMemcpyDefault)));
}

void DeepCopyAsyncCuda(void *dst, const void *src, size_t n) {
  cudaStream_t s = cuda_get_deep_copy_stream();
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      (CudaInternal::singleton().cuda_memcpy_async_wrapper(
          dst, src, n, cudaMemcpyDefault, s)));
  Impl::cuda_stream_synchronize(
      s,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          DeepCopyResourceSynchronization,
      "Kokkos::Impl::DeepCopyAsyncCuda: Deep Copy Stream Sync");
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
bool CudaUVMSpace::available() { return true; }
#endif

/*--------------------------------------------------------------------------*/

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
// The purpose of the following variable is to allow a state-based choice
// for pinning UVM allocations to the CPU. For now this is considered
// an experimental debugging capability - with the potential to work around
// some CUDA issues.
bool CudaUVMSpace::kokkos_impl_cuda_pin_uvm_to_host_v = false;

bool CudaUVMSpace::cuda_pin_uvm_to_host() {
  return CudaUVMSpace::kokkos_impl_cuda_pin_uvm_to_host_v;
}
void CudaUVMSpace::cuda_set_pin_uvm_to_host(bool val) {
  CudaUVMSpace::kokkos_impl_cuda_pin_uvm_to_host_v = val;
}
#endif
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
bool kokkos_impl_cuda_pin_uvm_to_host() {
  return Kokkos::CudaUVMSpace::cuda_pin_uvm_to_host();
}

void kokkos_impl_cuda_set_pin_uvm_to_host(bool val) {
  Kokkos::CudaUVMSpace::cuda_set_pin_uvm_to_host(val);
}
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

CudaSpace::CudaSpace() : m_device(Kokkos::Cuda().cuda_device()) {}

CudaUVMSpace::CudaUVMSpace() : m_device(Kokkos::Cuda().cuda_device()) {}

CudaHostPinnedSpace::CudaHostPinnedSpace() {}

size_t memory_threshold_g = 40000;  // 40 kB

//==============================================================================
// <editor-fold desc="allocate()"> {{{1

void *CudaSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void *CudaSpace::allocate(const Cuda &exec_space, const char *arg_label,
                          const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  return impl_allocate(exec_space, arg_label, arg_alloc_size, arg_logical_size);
}
void *CudaSpace::allocate(const char *arg_label, const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

namespace {
void *impl_allocate_common(const Cuda &exec_space, const char *arg_label,
                           const size_t arg_alloc_size,
                           const size_t arg_logical_size,
                           const Kokkos::Tools::SpaceHandle arg_handle,
                           bool exec_space_provided) {
  void *ptr = nullptr;

#ifndef CUDART_VERSION
#error CUDART_VERSION undefined!
#elif (defined(KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC) && CUDART_VERSION >= 11020)
  cudaError_t error_code;
  if (arg_alloc_size >= memory_threshold_g) {
    if (exec_space_provided) {
      error_code =
          exec_space.impl_internal_space_instance()->cuda_malloc_async_wrapper(
              &ptr, arg_alloc_size);
      exec_space.fence("Kokkos::Cuda: backend fence after async malloc");
    } else {
      error_code = Impl::CudaInternal::singleton().cuda_malloc_async_wrapper(
          &ptr, arg_alloc_size);
      Impl::cuda_device_synchronize(
          "Kokkos::Cuda: backend fence after async malloc");
    }
  } else {
    error_code =
        (exec_space_provided
             ? exec_space.impl_internal_space_instance()->cuda_malloc_wrapper(
                   &ptr, arg_alloc_size)
             : Impl::CudaInternal::singleton().cuda_malloc_wrapper(
                   &ptr, arg_alloc_size));
  }
#else
  cudaError_t error_code;
  if (exec_space_provided) {
    error_code = exec_space.impl_internal_space_instance()->cuda_malloc_wrapper(
        &ptr, arg_alloc_size);
  } else {
    error_code = Impl::CudaInternal::singleton().cuda_malloc_wrapper(
        &ptr, arg_alloc_size);
  }
#endif
  if (error_code != cudaSuccess) {  // TODO tag as unlikely branch
    // This is the only way to clear the last error, which
    // we should do here since we're turning it into an
    // exception here
    exec_space.impl_internal_space_instance()->cuda_get_last_error_wrapper();
    throw Experimental::CudaRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        Experimental::RawMemoryAllocationFailure::AllocationMechanism::
            CudaMalloc);
  }

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}
}  // namespace

void *CudaSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  return impl_allocate_common(Kokkos::Cuda{}, arg_label, arg_alloc_size,
                              arg_logical_size, arg_handle, false);
}

void *CudaSpace::impl_allocate(
    const Cuda &exec_space, const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  return impl_allocate_common(exec_space, arg_label, arg_alloc_size,
                              arg_logical_size, arg_handle, true);
}

void *CudaUVMSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void *CudaUVMSpace::allocate(const char *arg_label, const size_t arg_alloc_size,
                             const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void *CudaUVMSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void *ptr = nullptr;

  Cuda::impl_static_fence(
      "Kokkos::CudaUVMSpace::impl_allocate: Pre UVM Allocation");
  if (arg_alloc_size > 0) {
    Kokkos::Impl::num_uvm_allocations++;

    auto error_code =
        Impl::CudaInternal::singleton().cuda_malloc_managed_wrapper(
            &ptr, arg_alloc_size, cudaMemAttachGlobal);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
    if (Kokkos::CudaUVMSpace::cuda_pin_uvm_to_host())
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (Impl::CudaInternal::singleton().cuda_mem_advise_wrapper(
              ptr, arg_alloc_size, cudaMemAdviseSetPreferredLocation,
              cudaCpuDeviceId)));
#endif

    if (error_code != cudaSuccess) {  // TODO tag as unlikely branch
      // This is the only way to clear the last error, which
      // we should do here since we're turning it into an
      // exception here
      Impl::CudaInternal::singleton().cuda_get_last_error_wrapper();
      throw Experimental::CudaRawMemoryAllocationFailure(
          arg_alloc_size, error_code,
          Experimental::RawMemoryAllocationFailure::AllocationMechanism::
              CudaMallocManaged);
    }
  }
  Cuda::impl_static_fence(
      "Kokkos::CudaUVMSpace::impl_allocate: Post UVM Allocation");
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}
void *CudaHostPinnedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void *CudaHostPinnedSpace::allocate(const char *arg_label,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void *CudaHostPinnedSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void *ptr = nullptr;

  auto error_code = Impl::CudaInternal::singleton().cuda_host_alloc_wrapper(
      &ptr, arg_alloc_size, cudaHostAllocDefault);
  if (error_code != cudaSuccess) {  // TODO tag as unlikely branch
    // This is the only way to clear the last error, which
    // we should do here since we're turning it into an
    // exception here
    Impl::CudaInternal::singleton().cuda_get_last_error_wrapper();
    throw Experimental::CudaRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        Experimental::RawMemoryAllocationFailure::AllocationMechanism::
            CudaHostAlloc);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}

// </editor-fold> end allocate() }}}1
//==============================================================================
void CudaSpace::deallocate(void *const arg_alloc_ptr,
                           const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void CudaSpace::deallocate(const char *arg_label, void *const arg_alloc_ptr,
                           const size_t arg_alloc_size,
                           const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void CudaSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  try {
#ifndef CUDART_VERSION
#error CUDART_VERSION undefined!
#elif (defined(KOKKOS_ENABLE_IMPL_CUDA_MALLOC_ASYNC) && CUDART_VERSION >= 11020)
    if (arg_alloc_size >= memory_threshold_g) {
      Impl::cuda_device_synchronize(
          "Kokkos::Cuda: backend fence before async free");
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (Impl::CudaInternal::singleton().cuda_free_async_wrapper(
              arg_alloc_ptr)));
      Impl::cuda_device_synchronize(
          "Kokkos::Cuda: backend fence after async free");
    } else {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (Impl::CudaInternal::singleton().cuda_free_wrapper(arg_alloc_ptr)));
    }
#else
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (Impl::CudaInternal::singleton().cuda_free_wrapper(arg_alloc_ptr)));
#endif
  } catch (...) {
  }
}
void CudaUVMSpace::deallocate(void *const arg_alloc_ptr,
                              const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void CudaUVMSpace::deallocate(const char *arg_label, void *const arg_alloc_ptr,
                              const size_t arg_alloc_size

                              ,
                              const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void CudaUVMSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  Cuda::impl_static_fence(
      "Kokkos::CudaUVMSpace::impl_deallocate: Pre UVM Deallocation");
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  try {
    if (arg_alloc_ptr != nullptr) {
      Kokkos::Impl::num_uvm_allocations--;
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          (Impl::CudaInternal::singleton().cuda_free_wrapper(arg_alloc_ptr)));
    }
  } catch (...) {
  }
  Cuda::impl_static_fence(
      "Kokkos::CudaUVMSpace::impl_deallocate: Post UVM Deallocation");
}

void CudaHostPinnedSpace::deallocate(void *const arg_alloc_ptr,
                                     const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void CudaHostPinnedSpace::deallocate(const char *arg_label,
                                     void *const arg_alloc_ptr,
                                     const size_t arg_alloc_size,
                                     const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}

void CudaHostPinnedSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  try {
    KOKKOS_IMPL_CUDA_SAFE_CALL((
        Impl::CudaInternal::singleton().cuda_free_host_wrapper(arg_alloc_ptr)));
  } catch (...) {
  }
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::CudaSpace, void>::s_root_record;

SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::CudaUVMSpace, void>::s_root_record;

SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::CudaHostPinnedSpace, void>::s_root_record;
#endif

//==============================================================================
// <editor-fold desc="SharedAllocationRecord destructors"> {{{1

SharedAllocationRecord<Kokkos::CudaSpace, void>::~SharedAllocationRecord() {
  auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, (alloc_size - sizeof(SharedAllocationHeader)));
}

void SharedAllocationRecord<Kokkos::CudaSpace, void>::deep_copy_header_no_exec(
    void *ptr, const void *header) {
  Kokkos::Cuda exec;
  Kokkos::Impl::DeepCopy<CudaSpace, HostSpace>(exec, ptr, header,
                                               sizeof(SharedAllocationHeader));
  exec.fence(
      "SharedAllocationRecord<Kokkos::CudaSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

SharedAllocationRecord<Kokkos::CudaUVMSpace, void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size,
                     (SharedAllocationRecord<void, void>::m_alloc_size -
                      sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<Kokkos::CudaHostPinnedSpace,
                       void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size,
                     (SharedAllocationRecord<void, void>::m_alloc_size -
                      sizeof(SharedAllocationHeader)));
}

// </editor-fold> end SharedAllocationRecord destructors }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="SharedAllocationRecord constructors"> {{{1

SharedAllocationRecord<Kokkos::CudaSpace, void>::SharedAllocationRecord(
    const Kokkos::CudaSpace &arg_space, const std::string &arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::CudaSpace, void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {

  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // Copy to device memory
  Kokkos::Cuda exec;
  Kokkos::Impl::DeepCopy<CudaSpace, HostSpace>(
      exec, RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
  exec.fence(
      "SharedAllocationRecord<Kokkos::CudaSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

SharedAllocationRecord<Kokkos::CudaSpace, void>::SharedAllocationRecord(
    const Kokkos::Cuda &arg_exec_space, const Kokkos::CudaSpace &arg_space,
    const std::string &arg_label, const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::CudaSpace, void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_exec_space, arg_space,
                                               arg_label, arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {

  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // Copy to device memory
  Kokkos::Impl::DeepCopy<CudaSpace, HostSpace>(arg_exec_space,
                                               RecordBase::m_alloc_ptr, &header,
                                               sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<Kokkos::CudaUVMSpace, void>::SharedAllocationRecord(
    const Kokkos::CudaUVMSpace &arg_space, const std::string &arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::CudaUVMSpace, void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  this->base_t::_fill_host_accessible_header_info(*base_t::m_alloc_ptr,
                                                  arg_label);
}

SharedAllocationRecord<Kokkos::CudaHostPinnedSpace, void>::
    SharedAllocationRecord(
        const Kokkos::CudaHostPinnedSpace &arg_space,
        const std::string &arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::CudaHostPinnedSpace,
                                  void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  this->base_t::_fill_host_accessible_header_info(*base_t::m_alloc_ptr,
                                                  arg_label);
}

// </editor-fold> end SharedAllocationRecord constructors }}}1
//==============================================================================

void cuda_prefetch_pointer(const Cuda &space, const void *ptr, size_t bytes,
                           bool to_device) {
  if ((ptr == nullptr) || (bytes == 0)) return;
  cudaPointerAttributes attr;
  KOKKOS_IMPL_CUDA_SAFE_CALL((
      space.impl_internal_space_instance()->cuda_pointer_get_attributes_wrapper(
          &attr, ptr)));
  // I measured this and it turns out prefetching towards the host slows
  // DualView syncs down. Probably because the latency is not too bad in the
  // first place for the pull down. If we want to change that provde
  // cudaCpuDeviceId as the device if to_device is false
  bool is_managed = attr.type == cudaMemoryTypeManaged;
  if (to_device && is_managed &&
      space.cuda_device_prop().concurrentManagedAccess) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (space.impl_internal_space_instance()->cuda_mem_prefetch_async_wrapper(
            ptr, bytes, space.cuda_device())));
  }
}

}  // namespace Impl
}  // namespace Kokkos

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

namespace Kokkos {
namespace Impl {

// To avoid additional compilation cost for something that's (mostly?) not
// performance sensitive, we explicity instantiate these CRTP base classes here,
// where we have access to the associated *_timpl.hpp header files.
template class SharedAllocationRecordCommon<Kokkos::CudaSpace>;
template class HostInaccessibleSharedAllocationRecordCommon<Kokkos::CudaSpace>;
template class SharedAllocationRecordCommon<Kokkos::CudaUVMSpace>;
template class SharedAllocationRecordCommon<Kokkos::CudaHostPinnedSpace>;

}  // end namespace Impl
}  // end namespace Kokkos

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================

#else
void KOKKOS_CORE_SRC_CUDA_CUDASPACE_PREVENT_LINK_ERROR() {}
#endif  // KOKKOS_ENABLE_CUDA
