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
#ifndef KOKKOS_HBWSPACE_HPP
#define KOKKOS_HBWSPACE_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_HBWSPACE

#include <Kokkos_HostSpace.hpp>

namespace Kokkos {

namespace Experimental {

/// \class HBWSpace
/// \brief Memory management for host memory.
///
/// HBWSpace is a memory space that governs host memory.  "Host"
/// memory means the usual CPU-accessible memory.
class HBWSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space = HBWSpace;
  using size_type    = size_t;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
  using execution_space = Kokkos::DefaultHostExecutionSpace;

  //! This memory space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  /**\brief  Default memory space instance */
  HBWSpace();
  HBWSpace(const HBWSpace& rhs) = default;
  HBWSpace& operator=(const HBWSpace&) = default;
  ~HBWSpace()                          = default;

  /**\brief  Non-default memory space instance to choose allocation mechansim,
   * if available */

  enum AllocationMechanism {
    STD_MALLOC,
    POSIX_MEMALIGN,
    POSIX_MMAP,
    INTEL_MM_ALLOC
  };

  explicit HBWSpace(const AllocationMechanism&);

  /**\brief  Allocate untracked memory in the space */
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

 private:
  template <class, class, class, class>
  friend class LogicalMemorySpace;

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
  static constexpr const char* name() { return "HBW"; }

 private:
  AllocationMechanism m_alloc_mech;
  friend class Kokkos::Impl::SharedAllocationRecord<
      Kokkos::Experimental::HBWSpace, void>;
};

}  // namespace Experimental

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  friend Kokkos::Experimental::HBWSpace;

  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_ENABLE_DEBUG
  /**\brief  Root record for tracked allocations from this HBWSpace instance */
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HBWSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  SharedAllocationRecord(
      const Kokkos::Experimental::HBWSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  inline std::string get_label() const {
    return std::string(RecordBase::head()->m_label);
  }

  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::HBWSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size) {
    KOKKOS_IF_ON_HOST((return new SharedAllocationRecord(arg_space, arg_label,
                                                         arg_alloc_size);))
    KOKKOS_IF_ON_DEVICE(((void)arg_space; (void)arg_label; (void)arg_alloc_size;
                         return nullptr;))
  }

  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(const Kokkos::Experimental::HBWSpace& arg_space,
                                const std::string& arg_label,
                                const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::HBWSpace&,
                            bool detail = false);
};

}  // namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::HBWSpace,
                                    Kokkos::Experimental::HBWSpace>::assignable,
    "");

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::Experimental::HBWSpace> {
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HBWSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

}  // namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <>
struct DeepCopy<Kokkos::Experimental::HBWSpace, Kokkos::Experimental::HBWSpace,
                DefaultHostExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const DefaultHostExecutionSpace& exec, void* dst, const void* src,
           size_t n) {
    hostspace_parallel_deepcopy(exec, dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HBWSpace, Kokkos::Experimental::HBWSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<Kokkos::Experimental::HBWSpace, "
        "Kokkos::Experimental::HBWSpace,ExecutionSpace::DeepCopy: fence "
        "before copy");
    hostspace_parallel_deepcopy_async(dst, src, n);
  }
};

template <>
struct DeepCopy<HostSpace, Kokkos::Experimental::HBWSpace,
                DefaultHostExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const DefaultHostExecutionSpace& exec, void* dst, const void* src,
           size_t n) {
    hostspace_parallel_deepcopy(exec, dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::HBWSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::HBWSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    hostspace_parallel_deepcopy_async(copy_space, dst, src, n);
  }
};

template <>
struct DeepCopy<Kokkos::Experimental::HBWSpace, HostSpace,
                DefaultHostExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const DefaultHostExecutionSpace& exec, void* dst, const void* src,
           size_t n) {
    hostspace_parallel_deepcopy(exec, dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HBWSpace, HostSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    hostspace_parallel_deepcopy(dst, src, n);
  }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<Kokkos::Experimental::HBWSpace, HostSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    hostspace_parallel_deepcopy_async(dst, src, n);
  }
};

}  // namespace Impl

}  // namespace Kokkos

#endif
#endif  // #define KOKKOS_HBWSPACE_HPP
