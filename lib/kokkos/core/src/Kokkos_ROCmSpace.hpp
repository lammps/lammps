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

#ifndef KOKKOS_ROCMSPACE_HPP
#define KOKKOS_ROCMSPACE_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined(KOKKOS_ENABLE_ROCM)

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <Kokkos_HostSpace.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
/** \brief  ROCm on-device memory management */

class ROCmSpace {
 public:
  //! Tag this class as a kokkos memory space
  typedef ROCmSpace memory_space;
  typedef Kokkos::Experimental::ROCm execution_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  typedef unsigned int size_type;

  /*--------------------------------*/

  ROCmSpace();
  ROCmSpace(ROCmSpace&& rhs)      = default;
  ROCmSpace(const ROCmSpace& rhs) = default;
  ROCmSpace& operator=(ROCmSpace&& rhs) = default;
  ROCmSpace& operator=(const ROCmSpace& rhs) = default;
  ~ROCmSpace()                               = default;

  /**\brief  Allocate untracked memory in the rocm space */
  void* allocate(const size_t arg_alloc_size) const;

  /**\brief  Deallocate untracked memory in the rocm space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; };

  /*--------------------------------*/
  /** \brief  Error reporting for HostSpace attempt to access ROCmSpace */
  static void access_error();
  static void access_error(const void* const);

 private:
  int m_device;  ///< Which ROCm device

  static constexpr const char* m_name = "ROCm";
  friend class Kokkos::Impl::SharedAllocationRecord<
      Kokkos::Experimental::ROCmSpace, void>;
};

}  // namespace Experimental

namespace Impl {

void* rocm_device_allocate(int);
void* rocm_hostpinned_allocate(int);
void rocm_device_free(void*);

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
void init_lock_arrays_rocm_space();

/// \brief Retrieve the pointer to the lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* atomic_lock_array_rocm_space_ptr(bool deallocate = false);

/// \brief Retrieve the pointer to the scratch array for team and thread private
/// global memory.
///
/// Team and Thread private scratch allocations in
/// global memory are acquired via locks.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* scratch_lock_array_rocm_space_ptr(bool deallocate = false);

/// \brief Retrieve the pointer to the scratch array for unique identifiers.
///
/// Unique identifiers in the range 0-ROCm::concurrency
/// are provided via locks.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* threadid_lock_array_rocm_space_ptr(bool deallocate = false);
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
/** \brief  Host memory that is accessible to ROCm execution space
 *          through ROCm's host-pinned memory allocation.
 */
class ROCmHostPinnedSpace {
 public:
  //! Tag this class as a kokkos memory space
  /** \brief  Memory is in HostSpace so use the HostSpace::execution_space */
  typedef HostSpace::execution_space execution_space;
  typedef ROCmHostPinnedSpace memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;
  typedef unsigned int size_type;

  /*--------------------------------*/

  ROCmHostPinnedSpace();
  ROCmHostPinnedSpace(ROCmHostPinnedSpace&& rhs)      = default;
  ROCmHostPinnedSpace(const ROCmHostPinnedSpace& rhs) = default;
  ROCmHostPinnedSpace& operator=(ROCmHostPinnedSpace&& rhs) = default;
  ROCmHostPinnedSpace& operator=(const ROCmHostPinnedSpace& rhs) = default;
  ~ROCmHostPinnedSpace()                                         = default;

  /**\brief  Allocate untracked memory in the space */
  void* allocate(const size_t arg_alloc_size) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; };

 private:
  static constexpr const char* m_name = "ROCmHostPinned";

  /*--------------------------------*/
};
}  // namespace Experimental
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

static_assert(Kokkos::Impl::MemorySpaceAccess<
                  Kokkos::Experimental::ROCmSpace,
                  Kokkos::Experimental::ROCmSpace>::assignable,
              "");

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::Experimental::ROCmSpace> {
  enum { assignable = false };
  enum { accessible = false };
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::ROCmHostPinnedSpace> {
  // HostSpace::execution_space == ROCmHostPinnedSpace::execution_space
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::Experimental::ROCmSpace, Kokkos::HostSpace> {
  enum { assignable = false };
  enum { accessible = false };
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::ROCmSpace,
                         Kokkos::Experimental::ROCmHostPinnedSpace> {
  // ROCmSpace::execution_space != ROCmHostPinnedSpace::execution_space
  enum { assignable = false };
  enum { accessible = true };  // ROCmSpace::execution_space
  enum { deepcopy = true };
};

//----------------------------------------
// ROCmHostPinnedSpace::execution_space == HostSpace::execution_space
// ROCmHostPinnedSpace accessible to both ROCm and Host

template <>
struct MemorySpaceAccess<Kokkos::Experimental::ROCmHostPinnedSpace,
                         Kokkos::HostSpace> {
  enum { assignable = false };  // Cannot access from ROCm
  enum { accessible = true };   // ROCmHostPinnedSpace::execution_space
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::ROCmHostPinnedSpace,
                         Kokkos::Experimental::ROCmSpace> {
  enum { assignable = false };  // Cannot access from Host
  enum { accessible = false };
  enum { deepcopy = true };
};

};  // namespace Impl
//----------------------------------------

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

hc::completion_future DeepCopyAsyncROCm(void* dst, const void* src, size_t n);

template <>
struct DeepCopy<Kokkos::Experimental::ROCmSpace,
                Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace,
                Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace,
                Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmSpace,
                Kokkos::Experimental::ROCmSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::ROCmSpace,
                   Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCm>(
        dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    hc::completion_future fut = DeepCopyAsyncROCm(dst, src, n);
    fut.wait();
    //    DeepCopy (dst,src,n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopy(dst, src, n);
  }
};

template <>
struct DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
                Kokkos::Experimental::ROCmHostPinnedSpace,
                Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace, HostSpace,
                Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace,
                Kokkos::Experimental::ROCm> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::ROCm&, void* dst, const void* src,
           size_t);
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmSpace,
                Kokkos::Experimental::ROCmHostPinnedSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    hc::completion_future fut = DeepCopyAsyncROCm(dst, src, n);
    fut.wait();
    //    DeepCopyROCm (dst,src,n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
                Kokkos::Experimental::ROCmSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    hc::completion_future fut = DeepCopyAsyncROCm(dst, src, n);
    fut.wait();
    //    DeepCopyROCm (dst,src,n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
                Kokkos::Experimental::ROCmHostPinnedSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
                   Kokkos::Experimental::ROCmHostPinnedSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    //    hc::completion_future fut = DeepCopyAsyncROCm (dst,src,n);
    //    fut.wait();
    //    DeepCopyAsyncROCm (dst,src,n);
    DeepCopy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace, HostSpace,
                ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace, HostSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace,
                ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace,
                   Kokkos::Experimental::ROCm>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopy(dst, src, n);
  }
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** Running in ROCmSpace attempting to access HostSpace: error */
template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::Experimental::ROCmSpace,
                                           Kokkos::HostSpace> {
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify(void) {
    Kokkos::abort("ROCm code attempted to access HostSpace memory");
  }

  KOKKOS_INLINE_FUNCTION static void verify(const void*) {
    Kokkos::abort("ROCm code attempted to access HostSpace memory");
  }
};

/** Running in ROCmSpace accessing ROCmHostPinnedSpace: ok */
template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::Experimental::ROCmSpace,
    Kokkos::Experimental::ROCmHostPinnedSpace> {
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void*) {}
};

/** Running in ROCmSpace attempting to access an unknown space: error */
template <class OtherSpace>
struct VerifyExecutionCanAccessMemorySpace<
    typename std::enable_if<
        !is_same<Kokkos::Experimental::ROCmSpace, OtherSpace>::value,
        Kokkos::Experimental::ROCmSpace>::type,
    OtherSpace> {
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify(void) {
    Kokkos::abort("ROCm code attempted to access unknown Space memory");
  }

  KOKKOS_INLINE_FUNCTION static void verify(const void*) {
    Kokkos::abort("ROCm code attempted to access unknown Space memory");
  }
};

//----------------------------------------------------------------------------
/** Running in HostSpace attempting to access ROCmSpace */
template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::HostSpace,
                                           Kokkos::Experimental::ROCmSpace> {
  enum { value = false };
  inline static void verify(void) {
    Kokkos::Experimental::ROCmSpace::access_error();
  }
  inline static void verify(const void* p) {
    Kokkos::Experimental::ROCmSpace::access_error(p);
  }
};

/** Running in HostSpace accessing ROCmHostPinnedSpace is OK */
template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace> {
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void*) {}
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  typedef SharedAllocationRecord<void, void> RecordBase;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::ROCmSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord(
      const Kokkos::Experimental::ROCmSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  std::string get_label() const;

  static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::ROCmSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(
      const Kokkos::Experimental::ROCmSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::ROCmSpace&,
                            bool detail = false);
};

template <>
class SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  typedef SharedAllocationRecord<void, void> RecordBase;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::ROCmHostPinnedSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() : RecordBase(), m_space() {}

  SharedAllocationRecord(
      const Kokkos::Experimental::ROCmHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  std::string get_label() const;

  static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::ROCmHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);
  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(
      const Kokkos::Experimental::ROCmHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::ROCmHostPinnedSpace&,
                            bool detail = false);
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_ROCM ) */
#endif /* #define KOKKOS_ROCMSPACE_HPP */
