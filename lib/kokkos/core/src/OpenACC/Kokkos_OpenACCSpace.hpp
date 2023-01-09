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

#ifndef KOKKOS_OPENACC_SPACE_HPP
#define KOKKOS_OPENACC_SPACE_HPP

#include <Kokkos_Concepts.hpp>

#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include <openacc.h>
#include <iosfwd>

namespace Kokkos::Experimental {

class OpenACC;

class OpenACCSpace {
 public:
  using memory_space    = OpenACCSpace;
  using execution_space = OpenACC;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using size_type = size_t;

  OpenACCSpace() = default;

  /**\brief  Allocate untracked memory in the space */
  void* allocate(const Kokkos::Experimental::OpenACC& exec_space,
                 const size_t arg_alloc_size) const;
  void* allocate(const Kokkos::Experimental::OpenACC& exec_space,
                 const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr char const* name() { return "OpenACCSpace"; }

 private:
  void* impl_allocate(const Kokkos::Experimental::OpenACC& exec_space,
                      const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
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

}  // namespace Kokkos::Experimental

/*--------------------------------------------------------------------------*/

template <>
struct Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::Experimental::OpenACCSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                                       Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool { accessible = false };
  enum : bool { deepcopy = true };
};

template <>
struct Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                                       Kokkos::Experimental::OpenACCSpace> {
  enum : bool { assignable = true };
  enum : bool { accessible = true };
  enum : bool { deepcopy = true };
};
/*--------------------------------------------------------------------------*/

template <>
class Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace,
                                           void>
    : public HostInaccessibleSharedAllocationRecordCommon<
          Kokkos::Experimental::OpenACCSpace> {
 private:
  friend class HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenACCSpace>;
  friend class SharedAllocationRecordCommon<Kokkos::Experimental::OpenACCSpace>;
  friend Kokkos::Experimental::OpenACCSpace;

  using base_t = HostInaccessibleSharedAllocationRecordCommon<
      Kokkos::Experimental::OpenACCSpace>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  /**\brief  Root record for tracked allocations from this OpenACCSpace
   * instance */
  static RecordBase s_root_record;

  const Kokkos::Experimental::OpenACCSpace m_space;

 protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  template <typename ExecutionSpace>
  SharedAllocationRecord(
      const ExecutionSpace& /*exec_space*/,
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate)
      : SharedAllocationRecord(arg_space, arg_label, arg_alloc_size,
                               arg_dealloc) {}

  SharedAllocationRecord(
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate);

 public:
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const Kokkos::Experimental::OpenACCSpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size) {
    if (acc_on_device(acc_device_host)) {
      return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
    } else {
      return nullptr;
    }
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// FIXME_OPENACC: Need to update the DeepCopy implementations below to support
// multiple execution space instances.
// The current OpenACC backend implementation assumes that there is only one
// device execution space instance, and all the device operations (e.g., memory
// transfers, kernel launches, etc.) are implemented to be synchronous, which
// does not violate the Kokkos execution semantics with the single execution
// space instance.
template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::Experimental::OpenACCSpace,
                              ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // The behavior of acc_memcpy_device when bytes argument is zero is
    // clarified only in the latest OpenACC specification (V3.2), and thus the
    // value checking is added as a safeguard. (The current NVHPC (V22.5)
    // supports OpenACC V2.7.)
    if (n > 0) acc_memcpy_device(dst, const_cast<void*>(src), n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    if (n > 0) acc_memcpy_device(dst, const_cast<void*>(src), n);
  }
};

template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::HostSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) acc_memcpy_to_device(dst, const_cast<void*>(src), n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    if (n > 0) acc_memcpy_to_device(dst, const_cast<void*>(src), n);
  }
};

template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<
    Kokkos::HostSpace, Kokkos::Experimental::OpenACCSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) acc_memcpy_from_device(dst, const_cast<void*>(src), n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence();
    if (n > 0) acc_memcpy_from_device(dst, const_cast<void*>(src), n);
  }
};

#endif
