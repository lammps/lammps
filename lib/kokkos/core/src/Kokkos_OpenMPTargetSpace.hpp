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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
#endif
#ifndef KOKKOS_OPENMPTARGETSPACE_HPP
#define KOKKOS_OPENMPTARGETSPACE_HPP

#include <cstring>
#include <string>
#include <iosfwd>
#include <typeinfo>

#include <Kokkos_Core_fwd.hpp>

#ifdef KOKKOS_ENABLE_OPENMPTARGET

#include <OpenMPTarget/Kokkos_OpenMPTarget_Error.hpp>
#include <Kokkos_HostSpace.hpp>
#include <omp.h>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
// void init_lock_array_host_space();

/// \brief Acquire a lock for the address
///
/// This function tries to acquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully acquired the
/// function returns true. Otherwise it returns false.
// bool lock_address_host_space(void* ptr);

/// \brief Release lock for the address
///
/// This function releases the lock for the hash value derived
/// from the provided ptr. This function should only be called
/// after previously successfully acquiring a lock with
/// lock_address.
// void unlock_address_host_space(void* ptr);

}  // namespace Impl
}  // namespace Kokkos

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
  OpenMPTargetSpace(OpenMPTargetSpace&& rhs)      = default;
  OpenMPTargetSpace(const OpenMPTargetSpace& rhs) = default;
  OpenMPTargetSpace& operator=(OpenMPTargetSpace&&) = default;
  OpenMPTargetSpace& operator=(const OpenMPTargetSpace&) = default;
  ~OpenMPTargetSpace()                                   = default;

  /**\brief  Allocate untracked memory in the space */
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

  friend class Kokkos::Impl::SharedAllocationRecord<
      Kokkos::Experimental::OpenMPTargetSpace, void>;
};
}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>
    : public HostInaccessibleSharedAllocationRecordCommon<
          Kokkos::Experimental::OpenMPTargetSpace> {
 private:
  friend class HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenMPTargetSpace>;
  friend class SharedAllocationRecordCommon<
      Kokkos::Experimental::OpenMPTargetSpace>;
  friend Kokkos::Experimental::OpenMPTargetSpace;

  using base_t = HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenMPTargetSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  /**\brief  Root record for tracked allocations from this OpenMPTargetSpace
   * instance */
  static RecordBase s_root_record;

  const Kokkos::Experimental::OpenMPTargetSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  template <typename ExecutionSpace>
  SharedAllocationRecord(
      const ExecutionSpace& /*exec_space*/,
      const Kokkos::Experimental::OpenMPTargetSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate)
      : SharedAllocationRecord(arg_space, arg_label, arg_alloc_size,
                               arg_dealloc) {}

  SharedAllocationRecord(
      const Kokkos::Experimental::OpenMPTargetSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::OpenMPTargetSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc) {
    KOKKOS_IF_ON_HOST(
        (return new SharedAllocationRecord(arg_space, arg_label, arg_alloc);))
    KOKKOS_IF_ON_DEVICE(
        ((void)arg_space; (void)arg_label; (void)arg_alloc; return nullptr;))
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// TODO: implement all possible deep_copies
template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace,
                Kokkos::Experimental::OpenMPTargetSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // In the Release and RelWithDebInfo builds, the size of the memcpy should
    // be greater than zero to avoid error. omp_target_memcpy returns zero on
    // success.
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_default_device(),
                                       omp_get_default_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenMPTargetSpace, OpenMPTargetSpace>: fence "
        "before "
        "copy");
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_default_device(),
                                       omp_get_default_device()));
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace, HostSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_default_device(),
                                       omp_get_initial_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenMPTargetSpace, HostSpace>: fence before "
        "copy");
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_default_device(),
                                       omp_get_initial_device()));
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::OpenMPTargetSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_initial_device(),
                                       omp_get_default_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<HostSpace, OpenMPTargetSpace>: fence before "
        "copy");
    if (n > 0)
      OMPT_SAFE_CALL(omp_target_memcpy(dst, const_cast<void*>(src), n, 0, 0,
                                       omp_get_initial_device(),
                                       omp_get_default_device()));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif /* #define KOKKOS_OPENMPTARGETSPACE_HPP */
