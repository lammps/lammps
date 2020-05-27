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

#ifndef KOKKOS_HIPSPACE_HPP
#define KOKKOS_HIPSPACE_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined(KOKKOS_ENABLE_HIP)

#include <iosfwd>
#include <typeinfo>
#include <string>
#include <cstddef>
#include <iosfwd>

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>

#include <impl/Kokkos_Profiling_Interface.hpp>

#include <hip/hip_runtime_api.h>
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
/** \brief  HIP on-device memory management */

class HIPSpace {
 public:
  //! Tag this class as a kokkos memory space
  using memory_space    = HIPSpace;
  using execution_space = Kokkos::Experimental::HIP;
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
  void* allocate(const size_t arg_alloc_size) const;

  /**\brief  Deallocate untracked memory in the hip space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return "HIP"; }

  /*--------------------------------*/
  /** \brief  Error reporting for HostSpace attempt to access HIPSpace */
  static void access_error();
  static void access_error(const void* const);

 private:
  int m_device;  ///< Which HIP device

  friend class Kokkos::Impl::SharedAllocationRecord<
      Kokkos::Experimental::HIPSpace, void>;
};

}  // namespace Experimental

namespace Impl {

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
void init_lock_arrays_hip_space();

/// \brief Retrieve the pointer to the lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* atomic_lock_array_hip_space_ptr(bool deallocate = false);

/// \brief Retrieve the pointer to the scratch array for team and thread private
/// global memory.
///
/// Team and Thread private scratch allocations in
/// global memory are acquired via locks.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* scratch_lock_array_hip_space_ptr(bool deallocate = false);

/// \brief Retrieve the pointer to the scratch array for unique identifiers.
///
/// Unique identifiers in the range 0-HIP::concurrency
/// are provided via locks.
/// This function retrieves the lock array pointer.
/// If the array is not yet allocated it will do so.
int* threadid_lock_array_hip_space_ptr(bool deallocate = false);
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
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
  void* allocate(const size_t arg_alloc_size) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return "HIPHostPinned"; }

  /*--------------------------------*/
};
}  // namespace Experimental
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

static_assert(
    Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                                    Kokkos::Experimental::HIPSpace>::assignable,
    "");

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::HostSpace, Kokkos::Experimental::HIPSpace> {
  enum { assignable = false };
  enum { accessible = false };
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::HIPHostPinnedSpace> {
  // HostSpace::execution_space == HIPHostPinnedSpace::execution_space
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace, Kokkos::HostSpace> {
  enum { assignable = false };
  enum { accessible = false };
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                         Kokkos::Experimental::HIPHostPinnedSpace> {
  // HIPSpace::execution_space != HIPHostPinnedSpace::execution_space
  enum { assignable = false };
  enum { accessible = true };  // HIPSpace::execution_space
  enum { deepcopy = true };
};

//----------------------------------------
// HIPHostPinnedSpace::execution_space == HostSpace::execution_space
// HIPHostPinnedSpace accessible to both HIP and Host

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPHostPinnedSpace,
                         Kokkos::HostSpace> {
  enum { assignable = false };  // Cannot access from HIP
  enum { accessible = true };   // HIPHostPinnedSpace::execution_space
  enum { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPHostPinnedSpace,
                         Kokkos::Experimental::HIPSpace> {
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

void DeepCopyAsyncHIP(void* dst, const void* src, size_t n);

template <>
struct DeepCopy<Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIPSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIPSpace,
                ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::HIPSpace,
                   Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIP>(
        dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <>
struct DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
                Kokkos::Experimental::HIPHostPinnedSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace, HostSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <>
struct DeepCopy<HostSpace, Kokkos::Experimental::HIPHostPinnedSpace,
                Kokkos::Experimental::HIP> {
  DeepCopy(void* dst, const void* src, size_t);
  DeepCopy(const Kokkos::Experimental::HIP&, void* dst, const void* src,
           size_t);
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPSpace,
                Kokkos::Experimental::HIPHostPinnedSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
                Kokkos::Experimental::HIPSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
                Kokkos::Experimental::HIPHostPinnedSpace, ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
                   Kokkos::Experimental::HIPHostPinnedSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace, HostSpace,
                ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace, HostSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::HIPHostPinnedSpace,
                ExecutionSpace> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    (void)DeepCopy<HostSpace, Kokkos::Experimental::HIPHostPinnedSpace,
                   Kokkos::Experimental::HIP>(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence();
    DeepCopyAsyncHIP(dst, src, n);
  }
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** Running in HIPSpace attempting to access HostSpace: error */
template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::Experimental::HIPSpace,
                                           Kokkos::HostSpace> {
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify(void) {
    Kokkos::abort("HIP code attempted to access HostSpace memory");
  }

  KOKKOS_INLINE_FUNCTION static void verify(const void*) {
    Kokkos::abort("HIP code attempted to access HostSpace memory");
  }
};

/** Running in HIPSpace accessing HIPHostPinnedSpace: ok */
template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIPHostPinnedSpace> {
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void*) {}
};

/** Running in HIPSpace attempting to access an unknown space: error */
template <class OtherSpace>
struct VerifyExecutionCanAccessMemorySpace<
    typename std::enable_if<
        !std::is_same<Kokkos::Experimental::HIPSpace, OtherSpace>::value,
        Kokkos::Experimental::HIPSpace>::type,
    OtherSpace> {
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify(void) {
    Kokkos::abort("HIP code attempted to access unknown Space memory");
  }

  KOKKOS_INLINE_FUNCTION static void verify(const void*) {
    Kokkos::abort("HIP code attempted to access unknown Space memory");
  }
};

//----------------------------------------------------------------------------
/** Running in HostSpace attempting to access HIPSpace */
template <>
struct VerifyExecutionCanAccessMemorySpace<Kokkos::HostSpace,
                                           Kokkos::Experimental::HIPSpace> {
  enum { value = false };
  inline static void verify(void) {
    Kokkos::Experimental::HIPSpace::access_error();
  }
  inline static void verify(const void* p) {
    Kokkos::Experimental::HIPSpace::access_error(p);
  }
};

/** Running in HostSpace accessing HIPHostPinnedSpace is OK */
template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::HostSpace, Kokkos::Experimental::HIPHostPinnedSpace> {
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
class SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  typedef SharedAllocationRecord<void, void> RecordBase;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HIPSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord(
      const Kokkos::Experimental::HIPSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  std::string get_label() const;

  static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::HIPSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(const Kokkos::Experimental::HIPSpace& arg_space,
                                const std::string& arg_label,
                                const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::HIPSpace&,
                            bool detail = false);
};

template <>
class SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  typedef SharedAllocationRecord<void, void> RecordBase;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HIPHostPinnedSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() : RecordBase(), m_space() {}

  SharedAllocationRecord(
      const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  std::string get_label() const;

  static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);
  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(
      const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::HIPHostPinnedSpace&,
                            bool detail = false);
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {
class HIPInternal;
}
/// \class HIP
/// \brief Kokkos device for multicore processors in the host memory space.
class HIP {
 public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as a kokkos execution space
  using execution_space = HIP;
  using memory_space    = HIPSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = HIPSpace::size_type;

  using scratch_memory_space = ScratchMemorySpace<HIP>;

  ~HIP() = default;
  HIP();
  //  explicit HIP( const int instance_id );

  HIP(HIP&&)      = default;
  HIP(const HIP&) = default;
  HIP& operator=(HIP&&) = default;
  HIP& operator=(const HIP&) = default;

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__HIP_ARCH__)
    return true;
#else
    return false;
#endif
  }

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void impl_static_fence();
  void fence() const;

  /// \brief Print configuration information to the given output stream.
  static void print_configuration(std::ostream&, const bool detail = false);

  /// \brief Free any resources being consumed by the device.
  static void impl_finalize();

  /** \brief  Initialize the device.
   *
   */
  struct SelectDevice {
    int hip_device_id;
    SelectDevice() : hip_device_id(0) {}
    explicit SelectDevice(int id) : hip_device_id(id) {}
  };

  int hip_device() const;

  static void impl_initialize(const SelectDevice = SelectDevice());

  static int impl_is_initialized();

  //  static size_type device_arch();

  //  static size_type detect_device_count();

  static int concurrency();
  static const char* name();

  inline Impl::HIPInternal* impl_internal_space_instance() const {
    return m_space_instance;
  }

  uint32_t impl_instance_id() const noexcept { return 0; }

 private:
  Impl::HIPInternal* m_space_instance;
};
}  // namespace Experimental
namespace Profiling {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::Experimental::HIP> {
  static constexpr DeviceType id = DeviceType::HIP;
};
}  // namespace Experimental
}  // namespace Profiling
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                         Kokkos::Experimental::HIP::scratch_memory_space> {
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy = false };
};

template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::Experimental::HIP::memory_space,
    Kokkos::Experimental::HIP::scratch_memory_space> {
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void*) {}
};

template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::HostSpace, Kokkos::Experimental::HIP::scratch_memory_space> {
  enum { value = false };
  inline static void verify(void) {
    Kokkos::Experimental::HIPSpace::access_error();
  }
  inline static void verify(const void* p) {
    Kokkos::Experimental::HIPSpace::access_error(p);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #if defined( KOKKOS_ENABLE_HIP ) */
#endif /* #define KOKKOS_HIPSPACE_HPP */
