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
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_CUDASPACE_HPP
#define KOKKOS_CUDASPACE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Core_fwd.hpp>

#include <iosfwd>
#include <typeinfo>
#include <string>
#include <memory>

#include <Kokkos_HostSpace.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <impl/Kokkos_Profiling_Interface.hpp>

#include <Cuda/Kokkos_Cuda_abort.hpp>

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
extern "C" bool kokkos_impl_cuda_pin_uvm_to_host();
extern "C" void kokkos_impl_cuda_set_pin_uvm_to_host(bool);
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <typename T>
struct is_cuda_type_space : public std::false_type {};

}  // namespace Impl

/** \brief  Cuda on-device memory management */

class CudaSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space    = CudaSpace;
  using execution_space = Kokkos::Cuda;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using size_type = unsigned int;

  /*--------------------------------*/

  CudaSpace();

 private:
  CudaSpace(int device_id, cudaStream_t stream);

 public:
  CudaSpace(CudaSpace&& rhs)      = default;
  CudaSpace(const CudaSpace& rhs) = default;
  CudaSpace& operator=(CudaSpace&& rhs) = default;
  CudaSpace& operator=(const CudaSpace& rhs) = default;
  ~CudaSpace()                               = default;

  /**\brief  Allocate untracked memory in the cuda space */
  void* allocate(const Cuda& exec_space, const size_t arg_alloc_size) const;
  void* allocate(const Cuda& exec_space, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the cuda space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static CudaSpace impl_create(int device_id, cudaStream_t stream) {
    return CudaSpace(device_id, stream);
  }

 private:
  void* impl_allocate(const Cuda& exec_space, const char* arg_label,
                      const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
  void* impl_allocate(const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
  void impl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                       const size_t arg_alloc_size,
                       const size_t arg_logical_size = 0,
                       const Kokkos::Tools::SpaceHandle =
                           Kokkos::Tools::make_space_handle(name())) const;

 public:
  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; }

 private:
  int m_device;
  cudaStream_t m_stream;

  static constexpr const char* m_name = "Cuda";
};

template <>
struct Impl::is_cuda_type_space<CudaSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda memory that is accessible to Host execution space
 *          through Cuda's unified virtual memory (UVM) runtime.
 */
class CudaUVMSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space    = CudaUVMSpace;
  using execution_space = Cuda;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = unsigned int;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  /** \brief  If UVM capability is available */
  KOKKOS_DEPRECATED static bool available();
#endif

  /*--------------------------------*/

  /*--------------------------------*/

  CudaUVMSpace();

 private:
  CudaUVMSpace(int device_id, cudaStream_t stream);

 public:
  CudaUVMSpace(CudaUVMSpace&& rhs)      = default;
  CudaUVMSpace(const CudaUVMSpace& rhs) = default;
  CudaUVMSpace& operator=(CudaUVMSpace&& rhs) = default;
  CudaUVMSpace& operator=(const CudaUVMSpace& rhs) = default;
  ~CudaUVMSpace()                                  = default;

  /**\brief  Allocate untracked memory in the cuda space */
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace&, const size_t arg_alloc_size) const {
    return allocate(arg_alloc_size);
  }
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace&, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return allocate(arg_label, arg_alloc_size, arg_logical_size);
  }
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the cuda space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

 private:
  void* impl_allocate(const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
  void impl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                       const size_t arg_alloc_size,
                       const size_t arg_logical_size = 0,
                       const Kokkos::Tools::SpaceHandle =
                           Kokkos::Tools::make_space_handle(name())) const;

 public:
  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; }

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
  static bool cuda_pin_uvm_to_host();
  static void cuda_set_pin_uvm_to_host(bool val);
#endif
  /*--------------------------------*/

  static CudaUVMSpace impl_create(int device_id, cudaStream_t stream) {
    return CudaUVMSpace(device_id, stream);
  }

 private:
  int m_device;
  cudaStream_t m_stream;

#ifdef KOKKOS_IMPL_DEBUG_CUDA_PIN_UVM_TO_HOST
  static bool kokkos_impl_cuda_pin_uvm_to_host_v;
#endif
  static constexpr const char* m_name = "CudaUVM";
};

template <>
struct Impl::is_cuda_type_space<CudaUVMSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Host memory that is accessible to Cuda execution space
 *          through Cuda's host-pinned memory allocation.
 */
class CudaHostPinnedSpace {
 public:
  //! Tag this class as a kokkos memory space
  /** \brief  Memory is in HostSpace so use the HostSpace::execution_space */
  using execution_space = HostSpace::execution_space;
  using memory_space    = CudaHostPinnedSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = unsigned int;

  /*--------------------------------*/

  CudaHostPinnedSpace();

 private:
  CudaHostPinnedSpace(int device_id, cudaStream_t stream);

 public:
  CudaHostPinnedSpace(CudaHostPinnedSpace&& rhs)      = default;
  CudaHostPinnedSpace(const CudaHostPinnedSpace& rhs) = default;
  CudaHostPinnedSpace& operator=(CudaHostPinnedSpace&& rhs) = default;
  CudaHostPinnedSpace& operator=(const CudaHostPinnedSpace& rhs) = default;
  ~CudaHostPinnedSpace()                                         = default;

  /**\brief  Allocate untracked memory in the space */
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace&, const size_t arg_alloc_size) const {
    return allocate(arg_alloc_size);
  }
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace&, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return allocate(arg_label, arg_alloc_size, arg_logical_size);
  }
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static CudaHostPinnedSpace impl_create(int device_id, cudaStream_t stream) {
    return CudaHostPinnedSpace(device_id, stream);
  }

 private:
  void* impl_allocate(const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
  void impl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                       const size_t arg_alloc_size,
                       const size_t arg_logical_size = 0,
                       const Kokkos::Tools::SpaceHandle =
                           Kokkos::Tools::make_space_handle(name())) const;

 public:
  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; }

 private:
  int m_device;
  cudaStream_t m_stream;

  static constexpr const char* m_name = "CudaHostPinned";

  /*--------------------------------*/
};

template <>
struct Impl::is_cuda_type_space<CudaHostPinnedSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

cudaStream_t cuda_get_deep_copy_stream();

const std::unique_ptr<Kokkos::Cuda>& cuda_get_deep_copy_space(
    bool initialize = true);

static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaSpace,
                                              Kokkos::CudaSpace>::assignable);
static_assert(Kokkos::Impl::MemorySpaceAccess<
              Kokkos::CudaUVMSpace, Kokkos::CudaUVMSpace>::assignable);
static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::CudaHostPinnedSpace,
                                    Kokkos::CudaHostPinnedSpace>::assignable);

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::CudaSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::CudaUVMSpace> {
  // HostSpace::execution_space != CudaUVMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::CudaHostPinnedSpace> {
  // HostSpace::execution_space == CudaHostPinnedSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::CudaSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaSpace, Kokkos::CudaUVMSpace> {
  // CudaSpace::execution_space == CudaUVMSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaSpace, Kokkos::CudaHostPinnedSpace> {
  // CudaSpace::execution_space != CudaHostPinnedSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // CudaSpace::execution_space
  enum : bool { deepcopy = true };
};

//----------------------------------------
// CudaUVMSpace::execution_space == Cuda
// CudaUVMSpace accessible to both Cuda and Host

template <>
struct MemorySpaceAccess<Kokkos::CudaUVMSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };  // Cuda cannot access HostSpace
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaUVMSpace, Kokkos::CudaSpace> {
  // CudaUVMSpace::execution_space == CudaSpace::execution_space
  // Can access CudaUVMSpace from Host but cannot access CudaSpace from Host
  enum : bool { assignable = false };

  // CudaUVMSpace::execution_space can access CudaSpace
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaUVMSpace, Kokkos::CudaHostPinnedSpace> {
  // CudaUVMSpace::execution_space != CudaHostPinnedSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // CudaUVMSpace::execution_space
  enum : bool { deepcopy = true };
};

//----------------------------------------
// CudaHostPinnedSpace::execution_space == HostSpace::execution_space
// CudaHostPinnedSpace accessible to both Cuda and Host

template <>
struct MemorySpaceAccess<Kokkos::CudaHostPinnedSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };  // Cannot access from Cuda
  enum : bool { accessible = true };   // CudaHostPinnedSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaHostPinnedSpace, Kokkos::CudaSpace> {
  enum : bool { assignable = false };  // Cannot access from Host
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::CudaHostPinnedSpace, Kokkos::CudaUVMSpace> {
  enum : bool { assignable = false };  // different execution_space
  enum : bool { accessible = true };   // same accessibility
  enum : bool { deepcopy = true };
};

//----------------------------------------

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void DeepCopyCuda(void* dst, const void* src, size_t n);
void DeepCopyAsyncCuda(const Cuda& instance, void* dst, const void* src,
                       size_t n);
void DeepCopyAsyncCuda(void* dst, const void* src, size_t n);

template <class MemSpace>
struct DeepCopy<MemSpace, HostSpace, Cuda,
                std::enable_if_t<is_cuda_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyCuda(dst, src, n); }
  DeepCopy(const Cuda& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncCuda(instance, dst, src, n);
  }
};

template <class MemSpace>
struct DeepCopy<HostSpace, MemSpace, Cuda,
                std::enable_if_t<is_cuda_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyCuda(dst, src, n); }
  DeepCopy(const Cuda& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncCuda(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2>
struct DeepCopy<MemSpace1, MemSpace2, Cuda,
                std::enable_if_t<is_cuda_type_space<MemSpace1>::value &&
                                 is_cuda_type_space<MemSpace2>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyCuda(dst, src, n); }
  DeepCopy(const Cuda& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncCuda(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2, class ExecutionSpace>
struct DeepCopy<MemSpace1, MemSpace2, ExecutionSpace,
                std::enable_if_t<is_cuda_type_space<MemSpace1>::value &&
                                 is_cuda_type_space<MemSpace2>::value &&
                                 !std::is_same<ExecutionSpace, Cuda>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyCuda(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncCuda(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace1::name() + "Space, " +
        MemSpace2::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<MemSpace, HostSpace, ExecutionSpace,
                std::enable_if_t<is_cuda_type_space<MemSpace>::value &&
                                 !std::is_same<ExecutionSpace, Cuda>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyCuda(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncCuda(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace::name() +
        "Space, HostSpace, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<HostSpace, MemSpace, ExecutionSpace,
                std::enable_if_t<is_cuda_type_space<MemSpace>::value &&
                                 !std::is_same<ExecutionSpace, Cuda>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyCuda(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncCuda(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<HostSpace, ") + MemSpace::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
    Kokkos::CudaSpace);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::CudaUVMSpace);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::CudaHostPinnedSpace);

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #define KOKKOS_CUDASPACE_HPP */
