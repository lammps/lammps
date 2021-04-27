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
#include <SYCL/Kokkos_SYCL_Instance.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_Tools.hpp>

namespace Kokkos {
namespace Experimental {

class SYCLDeviceUSMSpace {
 public:
  using execution_space = SYCL;
  using memory_space    = SYCLDeviceUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = Impl::SYCLInternal::size_type;

  SYCLDeviceUSMSpace();

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
  static constexpr const char* name() { return "SYCLDeviceUSM"; };

 private:
  int m_device;
};
}  // namespace Experimental

namespace Impl {
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
struct MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

}  // namespace Impl

namespace Impl {

template <>
class SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>
    : public SharedAllocationRecord<void, void> {
 private:
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase*);

#ifdef KOKKOS_ENABLE_DEBUG
  static RecordBase s_root_record;
#endif

  const Kokkos::Experimental::SYCLDeviceUSMSpace m_space;

 protected:
  ~SharedAllocationRecord();

  SharedAllocationRecord(
      const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  std::string get_label() const;

  static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(
      const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size);

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size);

  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr);

  static SharedAllocationRecord* get_record(void* arg_alloc_ptr);

  static void print_records(std::ostream&,
                            const Kokkos::Experimental::SYCLDeviceUSMSpace&,
                            bool detail = false);
};

}  // namespace Impl

}  // namespace Kokkos

#endif
#endif
