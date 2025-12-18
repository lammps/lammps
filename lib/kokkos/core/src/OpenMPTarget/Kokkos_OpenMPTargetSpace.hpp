// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_OPENMPTARGETSPACE_HPP
#define KOKKOS_OPENMPTARGETSPACE_HPP

#include <cstring>
#include <string>
#include <iosfwd>
#include <typeinfo>

#include <Kokkos_Core_fwd.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_DeepCopy.hpp>

#ifdef KOKKOS_ENABLE_OPENMPTARGET

#include <OpenMPTarget/Kokkos_OpenMPTarget_Error.hpp>
#include <Kokkos_HostSpace.hpp>
#include <omp.h>

namespace Kokkos {
namespace Impl {

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::OpenMPTargetSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

//----------------------------------------
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {

/// \class OpenMPTargetSpace
/// \brief Memory management for host memory.
///
/// OpenMPTargetSpace is a memory space that governs host memory.  "Host"
/// memory means the usual CPU-accessible memory.
class OpenMPTargetSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space = OpenMPTargetSpace;
  using size_type    = unsigned;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
  using execution_space = Kokkos::Experimental::OpenMPTarget;

  //! This memory space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  /*--------------------------------*/

  /**\brief  Default memory space instance */
  OpenMPTargetSpace();
  OpenMPTargetSpace(OpenMPTargetSpace&& rhs)             = default;
  OpenMPTargetSpace(const OpenMPTargetSpace& rhs)        = default;
  OpenMPTargetSpace& operator=(OpenMPTargetSpace&&)      = default;
  OpenMPTargetSpace& operator=(const OpenMPTargetSpace&) = default;
  ~OpenMPTargetSpace()                                   = default;

  /**\brief  Allocate untracked memory in the space */
  // FIXME_OPENMPTARGET Use execution space instance
  void* allocate(const OpenMPTarget&, const size_t arg_alloc_size) const {
    return allocate(arg_alloc_size);
  }
  // FIXME_OPENMPTARGET Use execution space instance
  void* allocate(const OpenMPTarget&, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return allocate(arg_label, arg_alloc_size, arg_logical_size);
  }
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr const char* name() { return "OpenMPTargetSpace"; }

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
};
}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
    Kokkos::Experimental::OpenMPTargetSpace);

#endif
#endif /* #define KOKKOS_OPENMPTARGETSPACE_HPP */
