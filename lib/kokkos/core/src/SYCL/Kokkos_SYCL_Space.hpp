// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_SYCLSPACE_HPP
#define KOKKOS_SYCLSPACE_HPP

#include <Kokkos_Core_fwd.hpp>

#ifdef KOKKOS_ENABLE_SYCL
#include <Kokkos_Concepts.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <SYCL/Kokkos_SYCL_Instance.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_Tools.hpp>

namespace Kokkos {

namespace Impl {
template <typename T>
struct is_sycl_type_space : public std::false_type {};
}  // namespace Impl

class SYCLDeviceUSMSpace {
 public:
  using execution_space = SYCL;
  using memory_space    = SYCLDeviceUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = Impl::SYCLInternal::size_type;

  SYCLDeviceUSMSpace();
  explicit SYCLDeviceUSMSpace(sycl::queue queue);

  void* allocate(const SYCL& exec_space,
                 const std::size_t arg_alloc_size) const;
  void* allocate(const SYCL& exec_space, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;
  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr const char* name() { return "SYCLDeviceUSM"; }

 private:
  sycl::queue m_queue;
};

class SYCLSharedUSMSpace {
 public:
  using execution_space = SYCL;
  using memory_space    = SYCLSharedUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = Impl::SYCLInternal::size_type;

  SYCLSharedUSMSpace();
  explicit SYCLSharedUSMSpace(sycl::queue queue);

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
  void* allocate(const SYCL& exec_space,
                 const std::size_t arg_alloc_size) const;
  void* allocate(const SYCL& exec_space, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;
  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr const char* name() { return "SYCLSharedUSM"; }

 private:
  sycl::queue m_queue;
};

class SYCLHostUSMSpace {
 public:
  using execution_space = HostSpace::execution_space;
  using memory_space    = SYCLHostUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = Impl::SYCLInternal::size_type;

  SYCLHostUSMSpace();
  explicit SYCLHostUSMSpace(sycl::queue queue);

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
  void* allocate(const SYCL& exec_space,
                 const std::size_t arg_alloc_size) const;
  void* allocate(const SYCL& exec_space, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;
  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr const char* name() { return "SYCLHostUSM"; }

 private:
  sycl::queue m_queue;
};

namespace Impl {

template <>
struct is_sycl_type_space<Kokkos::SYCLDeviceUSMSpace> : public std::true_type {
};

template <>
struct is_sycl_type_space<Kokkos::SYCLSharedUSMSpace> : public std::true_type {
};

template <>
struct is_sycl_type_space<Kokkos::SYCLHostUSMSpace> : public std::true_type {};

static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                    Kokkos::SYCLDeviceUSMSpace>::assignable);

static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                    Kokkos::SYCLSharedUSMSpace>::assignable);

static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                    Kokkos::SYCLDeviceUSMSpace>::assignable);

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::SYCLDeviceUSMSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::SYCLSharedUSMSpace> {
  // HostSpace::execution_space != SYCLSharedUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::SYCLHostUSMSpace> {
  // HostSpace::execution_space ==
  // SYCLHostUSMSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                         Kokkos::SYCLSharedUSMSpace> {
  // SYCLDeviceUSMSpace::execution_space == SYCLSharedUSMSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace, Kokkos::SYCLHostUSMSpace> {
  // SYCLDeviceUSMSpace::execution_space !=
  // SYCLHostUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // SYCLDeviceUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

//----------------------------------------
// SYCLSharedUSMSpace::execution_space == SYCL
// SYCLSharedUSMSpace accessible to both SYCL and Host

template <>
struct MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };  // SYCL cannot access HostSpace
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                         Kokkos::SYCLDeviceUSMSpace> {
  // SYCLSharedUSMSpace::execution_space == SYCLDeviceUSMSpace::execution_space
  // Can access SYCLSharedUSMSpace from Host but cannot access
  // SYCLDeviceUSMSpace from Host
  enum : bool { assignable = false };

  // SYCLSharedUSMSpace::execution_space can access SYCLDeviceUSMSpace
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace, Kokkos::SYCLHostUSMSpace> {
  // SYCLSharedUSMSpace::execution_space !=
  // SYCLHostUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // SYCLSharedUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLHostUSMSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };  // Cannot access from SYCL
  enum : bool { accessible = true };   // SYCLHostUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLHostUSMSpace, Kokkos::SYCLDeviceUSMSpace> {
  enum : bool { assignable = false };  // Cannot access from Host
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLHostUSMSpace, Kokkos::SYCLSharedUSMSpace> {
  enum : bool { assignable = false };  // different execution_space
  enum : bool { accessible = true };   // same accessibility
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                         Kokkos::ScratchMemorySpace<Kokkos::SYCL>> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl

}  // namespace Kokkos

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
    Kokkos::SYCLDeviceUSMSpace);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::SYCLSharedUSMSpace);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(Kokkos::SYCLHostUSMSpace);

#endif
#endif
