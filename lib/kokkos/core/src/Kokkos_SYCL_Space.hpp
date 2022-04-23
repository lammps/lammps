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

namespace Experimental {

class SYCLDeviceUSMSpace {
 public:
  using execution_space = SYCL;
  using memory_space    = SYCLDeviceUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = Impl::SYCLInternal::size_type;

  SYCLDeviceUSMSpace();
  explicit SYCLDeviceUSMSpace(sycl::queue queue);

  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

 private:
  template <class, class, class, class>
  friend class LogicalMemorySpace;

 public:
  static constexpr const char* name() { return "SYCLDeviceUSM"; };

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

  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

 private:
  template <class, class, class, class>
  friend class LogicalMemorySpace;

 public:
  static constexpr const char* name() { return "SYCLSharedUSM"; };

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

  void* allocate(const std::size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  void deallocate(void* const arg_alloc_ptr,
                  const std::size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

 private:
  template <class, class, class, class>
  friend class LogicalMemorySpace;

 public:
  static constexpr const char* name() { return "SYCLHostUSM"; };

 private:
  sycl::queue m_queue;
};

}  // namespace Experimental

namespace Impl {

template <>
struct is_sycl_type_space<Kokkos::Experimental::SYCLDeviceUSMSpace>
    : public std::true_type {};

template <>
struct is_sycl_type_space<Kokkos::Experimental::SYCLSharedUSMSpace>
    : public std::true_type {};

template <>
struct is_sycl_type_space<Kokkos::Experimental::SYCLHostUSMSpace>
    : public std::true_type {};

static_assert(Kokkos::Impl::MemorySpaceAccess<
                  Kokkos::Experimental::SYCLDeviceUSMSpace,
                  Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
              "");

static_assert(Kokkos::Impl::MemorySpaceAccess<
                  Kokkos::Experimental::SYCLSharedUSMSpace,
                  Kokkos::Experimental::SYCLSharedUSMSpace>::assignable,
              "");

static_assert(Kokkos::Impl::MemorySpaceAccess<
                  Kokkos::Experimental::SYCLDeviceUSMSpace,
                  Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
              "");

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::SYCLDeviceUSMSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::SYCLSharedUSMSpace> {
  // HostSpace::execution_space != SYCLSharedUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::SYCLHostUSMSpace> {
  // HostSpace::execution_space ==
  // Experimental::SYCLHostUSMSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::Experimental::SYCLSharedUSMSpace> {
  // SYCLDeviceUSMSpace::execution_space == SYCLSharedUSMSpace::execution_space
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::Experimental::SYCLHostUSMSpace> {
  // Experimental::SYCLDeviceUSMSpace::execution_space !=
  // Experimental::SYCLHostUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool {
    accessible = true
  };  // Experimental::SYCLDeviceUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

//----------------------------------------
// SYCLSharedUSMSpace::execution_space == SYCL
// SYCLSharedUSMSpace accessible to both SYCL and Host

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLSharedUSMSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };  // SYCL cannot access HostSpace
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLSharedUSMSpace,
                         Kokkos::Experimental::SYCLDeviceUSMSpace> {
  // SYCLSharedUSMSpace::execution_space == SYCLDeviceUSMSpace::execution_space
  // Can access SYCLSharedUSMSpace from Host but cannot access
  // SYCLDeviceUSMSpace from Host
  enum : bool { assignable = false };

  // SYCLSharedUSMSpace::execution_space can access SYCLDeviceUSMSpace
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLSharedUSMSpace,
                         Kokkos::Experimental::SYCLHostUSMSpace> {
  // Experimental::SYCLSharedUSMSpace::execution_space !=
  // Experimental::SYCLHostUSMSpace::execution_space
  enum : bool { assignable = false };
  enum : bool {
    accessible = true
  };  // Experimental::SYCLSharedUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };  // Cannot access from SYCL
  enum : bool {
    accessible = true
  };  // Experimental::SYCLHostUSMSpace::execution_space
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                         Kokkos::Experimental::SYCLDeviceUSMSpace> {
  enum : bool { assignable = false };  // Cannot access from Host
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                         Kokkos::Experimental::SYCLSharedUSMSpace> {
  enum : bool { assignable = false };  // different execution_space
  enum : bool { accessible = true };   // same accessibility
  enum : bool { deepcopy = true };
};

template <>
struct MemorySpaceAccess<
    Kokkos::Experimental::SYCLDeviceUSMSpace,
    Kokkos::ScratchMemorySpace<Kokkos::Experimental::SYCL>> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl

namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>
    : public HostInaccessibleSharedAllocationRecordCommon<
          Kokkos::Experimental::SYCLDeviceUSMSpace> {
 private:
  friend class SharedAllocationRecordCommon<
      Kokkos::Experimental::SYCLDeviceUSMSpace>;
  friend class HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::SYCLDeviceUSMSpace>;
  using base_t = HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::SYCLDeviceUSMSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord(SharedAllocationRecord&&)      = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(SharedAllocationRecord&&) = delete;

#ifdef KOKKOS_ENABLE_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::SYCLDeviceUSMSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord(
      const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &base_t::deallocate);
};

template <>
class SharedAllocationRecord<Kokkos::Experimental::SYCLSharedUSMSpace, void>
    : public SharedAllocationRecordCommon<
          Kokkos::Experimental::SYCLSharedUSMSpace> {
 private:
  friend class SharedAllocationRecordCommon<
      Kokkos::Experimental::SYCLSharedUSMSpace>;
  using base_t =
      SharedAllocationRecordCommon<Kokkos::Experimental::SYCLSharedUSMSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord(SharedAllocationRecord&&)      = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(SharedAllocationRecord&&) = delete;

  static RecordBase s_root_record;

  const Kokkos::Experimental::SYCLSharedUSMSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord() = default;

  SharedAllocationRecord(
      const Kokkos::Experimental::SYCLSharedUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &base_t::deallocate);
};

template <>
class SharedAllocationRecord<Kokkos::Experimental::SYCLHostUSMSpace, void>
    : public SharedAllocationRecordCommon<
          Kokkos::Experimental::SYCLHostUSMSpace> {
 private:
  friend class SharedAllocationRecordCommon<
      Kokkos::Experimental::SYCLHostUSMSpace>;
  using base_t =
      SharedAllocationRecordCommon<Kokkos::Experimental::SYCLHostUSMSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord(SharedAllocationRecord&&)      = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(SharedAllocationRecord&&) = delete;

  static RecordBase s_root_record;

  const Kokkos::Experimental::SYCLHostUSMSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord() = default;

  SharedAllocationRecord(
      const Kokkos::Experimental::SYCLHostUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &base_t::deallocate);
};

}  // namespace Impl

}  // namespace Kokkos

#endif
#endif
