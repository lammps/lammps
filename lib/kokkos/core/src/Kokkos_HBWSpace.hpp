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

#ifndef KOKKOS_HBWSPACE_HPP
#define KOKKOS_HBWSPACE_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_HBWSPACE

#include <Kokkos_HostSpace.hpp>

namespace Kokkos {

namespace Experimental {

namespace Impl {

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
void init_lock_array_hbw_space();

/// \brief Acquire a lock for the address
///
/// This function tries to acquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully acquired the
/// function returns true. Otherwise it returns false.
bool lock_address_hbw_space(void* ptr);

/// \brief Release lock for the address
///
/// This function releases the lock for the hash value derived
/// from the provided ptr. This function should only be called
/// after previously successfully acquiring a lock with
/// lock_address.
void unlock_address_hbw_space(void* ptr);

}  // namespace Impl

}  // namespace Experimental

}  // namespace Kokkos

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
  typedef HBWSpace memory_space;
  typedef size_t size_type;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
#if defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP)
  typedef Kokkos::OpenMP execution_space;
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS)
  typedef Kokkos::Threads execution_space;
#elif defined(KOKKOS_ENABLE_OPENMP)
  typedef Kokkos::OpenMP execution_space;
#elif defined(KOKKOS_ENABLE_THREADS)
  typedef Kokkos::Threads execution_space;
#elif defined(KOKKOS_ENABLE_SERIAL)
  typedef Kokkos::Serial execution_space;
#else
#error \
    "At least one of the following host execution spaces must be defined: Kokkos::OpenMP, Kokkos::Threads, or Kokkos::Serial.  You might be seeing this message if you disabled the Kokkos::Serial device explicitly using the Kokkos_ENABLE_Serial:BOOL=OFF CMake option, but did not enable any of the other host execution space devices."
#endif

  //! This memory space preferred device_type
  typedef Kokkos::Device<execution_space, memory_space> device_type;

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

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;

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

  typedef SharedAllocationRecord<void, void> RecordBase;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_DEBUG
  /**\brief  Root record for tracked allocations from this HBWSpace instance */
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HBWSpace m_space;

 protected:
  ~SharedAllocationRecord()
#if defined( \
    KOKKOS_IMPL_INTEL_WORKAROUND_NOEXCEPT_SPECIFICATION_VIRTUAL_FUNCTION)
      noexcept
#endif
      ;
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
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
#else
    return (SharedAllocationRecord*)0;
#endif
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
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HBWSpace, Kokkos::HostSpace> {
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy = true };
};

}  // namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HBWSpace, Kokkos::Experimental::HBWSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) { memcpy(dst, src, n); }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    memcpy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::HBWSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) { memcpy(dst, src, n); }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    memcpy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HBWSpace, HostSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) { memcpy(dst, src, n); }

  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    memcpy(dst, src, n);
  }
};

}  // namespace Impl

}  // namespace Kokkos

namespace Kokkos {

namespace Impl {

template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::HostSpace,
                                           Kokkos::Experimental::HBWSpace> {
  enum { value = true };
  inline static void verify(void) {}
  inline static void verify(const void*) {}
};

template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::Experimental::HBWSpace,
                                           Kokkos::HostSpace> {
  enum { value = true };
  inline static void verify(void) {}
  inline static void verify(const void*) {}
};

}  // namespace Impl

}  // namespace Kokkos

#endif
#endif  // #define KOKKOS_HBWSPACE_HPP
