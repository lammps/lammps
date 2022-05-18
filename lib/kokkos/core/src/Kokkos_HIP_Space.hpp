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
#include <HIP/Kokkos_HIP_Error.hpp>  // HIP_SAFE_CALL

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_ExecSpaceInitializer.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>

#include <hip/hip_runtime_api.h>
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <typename T>
struct is_hip_type_space : public std::false_type {};

}  // namespace Impl

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
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the hip space */
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
  static constexpr const char* name() { return "HIP"; }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  /*--------------------------------*/
  /** \brief  Error reporting for HostSpace attempt to access HIPSpace */
  KOKKOS_DEPRECATED static void access_error();
  KOKKOS_DEPRECATED static void access_error(const void* const);
#endif

 private:
  int m_device;  ///< Which HIP device

  friend class Kokkos::Impl::SharedAllocationRecord<
      Kokkos::Experimental::HIPSpace, void>;
};

}  // namespace Experimental

template <>
struct Impl::is_hip_type_space<Experimental::HIPSpace> : public std::true_type {
};

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
  static constexpr const char* name() { return "HIPHostPinned"; }

  /*--------------------------------*/
};
}  // namespace Experimental

template <>
struct Impl::is_hip_type_space<Experimental::HIPHostPinnedSpace>
    : public std::true_type {};

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
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::HIPHostPinnedSpace> {
  // HostSpace::execution_space == HIPHostPinnedSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

//----------------------------------------

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace, Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                         Kokkos::Experimental::HIPHostPinnedSpace> {
  // HIPSpace::execution_space != HIPHostPinnedSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };  // HIPSpace::execution_space
  enum : bool { deepcopy = true };
};

//----------------------------------------
// HIPHostPinnedSpace::execution_space == HostSpace::execution_space
// HIPHostPinnedSpace accessible to both HIP and Host

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPHostPinnedSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };  // Cannot access from HIP
  enum : bool { accessible = true };   // HIPHostPinnedSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPHostPinnedSpace,
                         Kokkos::Experimental::HIPSpace> {
  enum : bool { assignable = false };  // Cannot access from Host
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

};  // namespace Impl
//----------------------------------------

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void DeepCopyHIP(void* dst, const void* src, size_t n);
void DeepCopyAsyncHIP(const Kokkos::Experimental::HIP& instance, void* dst,
                      const void* src, size_t n);
void DeepCopyAsyncHIP(void* dst, const void* src, size_t n);

template <class MemSpace>
struct DeepCopy<MemSpace, HostSpace, Kokkos::Experimental::HIP,
                std::enable_if_t<is_hip_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace>
struct DeepCopy<HostSpace, MemSpace, Kokkos::Experimental::HIP,
                std::enable_if_t<is_hip_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2>
struct DeepCopy<MemSpace1, MemSpace2, Kokkos::Experimental::HIP,
                std::enable_if_t<is_hip_type_space<MemSpace1>::value &&
                                 is_hip_type_space<MemSpace2>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2, class ExecutionSpace>
struct DeepCopy<
    MemSpace1, MemSpace2, ExecutionSpace,
    std::enable_if_t<
        is_hip_type_space<MemSpace1>::value &&
        is_hip_type_space<MemSpace2>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
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
struct DeepCopy<
    MemSpace, HostSpace, ExecutionSpace,
    std::enable_if_t<
        is_hip_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
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
struct DeepCopy<
    HostSpace, MemSpace, ExecutionSpace,
    std::enable_if_t<
        is_hip_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
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

namespace Kokkos {
namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>
    : public HostInaccessibleSharedAllocationRecordCommon<
          Kokkos::Experimental::HIPSpace> {
 private:
  friend class SharedAllocationRecordCommon<Kokkos::Experimental::HIPSpace>;
  friend class HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::HIPSpace>;
  using base_t = HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::HIPSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

#ifdef KOKKOS_ENABLE_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HIPSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord(
      const Kokkos::Experimental::HIPSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &base_t::deallocate);
};

template <>
class SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>
    : public SharedAllocationRecordCommon<
          Kokkos::Experimental::HIPHostPinnedSpace> {
 private:
  friend class SharedAllocationRecordCommon<
      Kokkos::Experimental::HIPHostPinnedSpace>;
  using base_t =
      SharedAllocationRecordCommon<Kokkos::Experimental::HIPHostPinnedSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

#ifdef KOKKOS_ENABLE_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::HIPHostPinnedSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  SharedAllocationRecord(
      const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &base_t::deallocate);
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

  HIP();
  HIP(hipStream_t stream, bool manage_stream = false);

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__HIP_DEVICE_COMPILE__)
    return true;
#else
    return false;
#endif
  }

  /** \brief Wait until all dispatched functors complete.
   *
   * The parallel_for or parallel_reduce dispatch of a functor may return
   * asynchronously, before the functor completes. This method does not return
   * until all dispatched functors on this device have completed.
   */
  static void impl_static_fence();
  static void impl_static_fence(const std::string&);

  void fence() const;
  void fence(const std::string&) const;

  hipStream_t hip_stream() const;

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
  static hipDeviceProp_t const& hip_device_prop();

  static void impl_initialize(const SelectDevice = SelectDevice());

  static int impl_is_initialized();

  //  static size_type device_arch();

  static size_type detect_device_count();

  static int concurrency();
  static const char* name();

  inline Impl::HIPInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

  uint32_t impl_instance_id() const noexcept;

 private:
  Kokkos::Impl::HostSharedPtr<Impl::HIPInternal> m_space_instance;
};
}  // namespace Experimental
namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::Experimental::HIP> {
  static constexpr DeviceType id = DeviceType::HIP;
};
}  // namespace Experimental
}  // namespace Tools

namespace Impl {

class HIPSpaceInitializer : public Kokkos::Impl::ExecSpaceInitializerBase {
 public:
  HIPSpaceInitializer()  = default;
  ~HIPSpaceInitializer() = default;
  void initialize(const InitArguments& args) final;
  void finalize(const bool) final;
  void fence() final;
  void fence(const std::string&) final;
  void print_configuration(std::ostream& msg, const bool detail) final;
};

template <class DT, class... DP>
struct ZeroMemset<Kokkos::Experimental::HIP, DT, DP...> {
  ZeroMemset(const Kokkos::Experimental::HIP& exec_space,
             const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipMemsetAsync(
        dst.data(), 0,
        dst.size() * sizeof(typename View<DT, DP...>::value_type),
        exec_space.hip_stream()));
  }

  ZeroMemset(const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipMemset(dst.data(), 0,
                  dst.size() * sizeof(typename View<DT, DP...>::value_type)));
  }
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                         Kokkos::Experimental::HIP::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #if defined( KOKKOS_ENABLE_HIP ) */
#endif /* #define KOKKOS_HIPSPACE_HPP */
