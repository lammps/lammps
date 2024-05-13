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

#ifndef KOKKOS_HIPSPACE_HPP
#define KOKKOS_HIPSPACE_HPP

#include <Kokkos_Core_fwd.hpp>

#include <iosfwd>
#include <typeinfo>
#include <string>
#include <cstddef>
#include <iosfwd>

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>  // HIP_SAFE_CALL

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>

#include <hip/hip_runtime_api.h>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <typename T>
struct is_hip_type_space : public std::false_type {};

}  // namespace Impl

/** \brief  HIP on-device memory management */

class HIPSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space    = HIPSpace;
  using execution_space = HIP;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using size_type = unsigned int;

  /*--------------------------------*/

  HIPSpace();
  HIPSpace(HIPSpace&& rhs)      = default;
  HIPSpace(const HIPSpace& rhs) = default;
  HIPSpace& operator=(HIPSpace&& rhs) = default;
  HIPSpace& operator=(const HIPSpace& rhs) = default;
  ~HIPSpace()                              = default;

  /**\brief  Allocate untracked memory in the hip space */
  // FIXME_HIP Use execution space instance
  void* allocate(const HIP&, const size_t arg_alloc_size) const {
    return allocate(arg_alloc_size);
  }
  // FIXME_HIP Use execution space instance
  void* allocate(const HIP&, const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return allocate(arg_label, arg_alloc_size, arg_logical_size);
  }
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the hip space */
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
  static constexpr const char* name() { return "HIP"; }

 private:
  int m_device;  ///< Which HIP device
};

template <>
struct Impl::is_hip_type_space<HIPSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
/** \brief  Host memory that is accessible to HIP execution space
 *          through HIP's host-pinned memory allocation.
 */
class HIPHostPinnedSpace {
 public:
  //! Tag this class as a kokkos memory space
  /** \brief  Memory is in HostSpace so use the HostSpace::execution_space */
  using execution_space = HostSpace::execution_space;
  using memory_space    = HIPHostPinnedSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = unsigned int;

  /*--------------------------------*/

  HIPHostPinnedSpace();
  HIPHostPinnedSpace(HIPHostPinnedSpace&& rhs)      = default;
  HIPHostPinnedSpace(const HIPHostPinnedSpace& rhs) = default;
  HIPHostPinnedSpace& operator=(HIPHostPinnedSpace&& rhs) = default;
  HIPHostPinnedSpace& operator=(const HIPHostPinnedSpace& rhs) = default;
  ~HIPHostPinnedSpace()                                        = default;

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
  static constexpr const char* name() { return "HIPHostPinned"; }

  /*--------------------------------*/
};

template <>
struct Impl::is_hip_type_space<HIPHostPinnedSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
/** \brief  Memory that is accessible to HIP execution space
 *          and host through HIP's memory page migration.
 */
class HIPManagedSpace {
 public:
  //! Tag this class as a kokkos memory space
  /** \brief  Memory is unified to both device and host via page migration
   *  and therefore able to be used by HostSpace::execution_space and
   *  DeviceSpace::execution_space.
   */
  //! tag this class as a kokkos memory space
  using memory_space    = HIPManagedSpace;
  using execution_space = HIP;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = unsigned int;

  /*--------------------------------*/

  HIPManagedSpace();
  HIPManagedSpace(HIPManagedSpace&& rhs)      = default;
  HIPManagedSpace(const HIPManagedSpace& rhs) = default;
  HIPManagedSpace& operator=(HIPManagedSpace&& rhs) = default;
  HIPManagedSpace& operator=(const HIPManagedSpace& rhs) = default;
  ~HIPManagedSpace()                                     = default;

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

  //  internal only method to determine whether page migration is supported
  bool impl_hip_driver_check_page_migration() const;

 private:
  int m_device;  ///< Which HIP device
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
  static constexpr const char* name() { return "HIPManaged"; }

  /*--------------------------------*/
};

template <>
struct Impl::is_hip_type_space<HIPManagedSpace> : public std::true_type {};

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

static_assert(Kokkos::Impl::MemorySpaceAccess<HIPSpace, HIPSpace>::assignable);

//----------------------------------------

template <>
struct MemorySpaceAccess<HostSpace, HIPSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HostSpace, HIPHostPinnedSpace> {
  // HostSpace::execution_space == HIPHostPinnedSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HostSpace, HIPManagedSpace> {
  // HostSpace::execution_space != HIPManagedSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<HIPSpace, HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPSpace, HIPHostPinnedSpace> {
  // HIPSpace::execution_space != HIPHostPinnedSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // HIPSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPSpace, HIPManagedSpace> {
  // HIPSpace::execution_space == HIPManagedSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

//----------------------------------------
// HIPHostPinnedSpace::execution_space == HostSpace::execution_space
// HIPHostPinnedSpace accessible to both HIP and Host

template <>
struct MemorySpaceAccess<HIPHostPinnedSpace, HostSpace> {
  enum : bool { assignable = false };  // Cannot access from HIP
  enum : bool { accessible = true };   // HIPHostPinnedSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPHostPinnedSpace, HIPSpace> {
  enum : bool { assignable = false };  // Cannot access from Host
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPHostPinnedSpace, HIPManagedSpace> {
  enum : bool { assignable = false };  // different exec_space
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

//----------------------------------------
// HIPManagedSpace::execution_space != HostSpace::execution_space
// HIPManagedSpace accessible to both HIP and Host

template <>
struct MemorySpaceAccess<HIPManagedSpace, HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };  // HIPHostPinnedSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPManagedSpace, HIPSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<HIPManagedSpace, HIPHostPinnedSpace> {
  enum : bool { assignable = false };  // different exec_space
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_HIPSPACE_HPP */
