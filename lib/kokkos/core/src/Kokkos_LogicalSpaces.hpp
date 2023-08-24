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
#ifndef KOKKOS_LOGICALSPACES_HPP
#define KOKKOS_LOGICALSPACES_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_Profiling.hpp>
#include <cstring>
namespace Kokkos {
namespace Experimental {
struct DefaultMemorySpaceNamer {
  static constexpr const char* get_name() {
    return "DefaultLogicalMemorySpaceName";
  }
};

struct LogicalSpaceSharesAccess {
  struct shared_access {};
  struct no_shared_access {};
};

/// \class LogicalMemorySpace
/// \brief
///
/// LogicalMemorySpace is a space that is identical to another space,
/// but differentiable by name and template argument
template <class BaseSpace, class DefaultBaseExecutionSpace = void,
          class Namer                = DefaultMemorySpaceNamer,
          class SharesAccessWithBase = LogicalSpaceSharesAccess::shared_access>
class LogicalMemorySpace {
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  // [DZP] For some reason I don't yet know, using LogicalMemorySpaces
  // inside an OpenMPTarget build causes errors in the
  // SharedAllocationRecords of other types. This is my way of erroring
  // a build if we instantiate a LogicalMemSpace in an OMPTarget build
  static_assert(!std::is_same<BaseSpace, BaseSpace>::value,
                "Can't use LogicalMemorySpaces in an OpenMPTarget build, we're "
                "debugging memory issues");
#endif
 public:
  //! Tag this class as a kokkos memory space
  using memory_space = LogicalMemorySpace<BaseSpace, DefaultBaseExecutionSpace,
                                          Namer, SharesAccessWithBase>;
  using size_type    = typename BaseSpace::size_type;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).

  using execution_space =
      std::conditional_t<std::is_void<DefaultBaseExecutionSpace>::value,
                         typename BaseSpace::execution_space,
                         DefaultBaseExecutionSpace>;

  using device_type = Kokkos::Device<execution_space, memory_space>;

  LogicalMemorySpace() = default;

  template <typename... Args>
  LogicalMemorySpace(Args&&... args) : underlying_space((Args &&) args...) {}

  /**\brief  Allocate untracked memory in the space */
  void* allocate(const size_t arg_alloc_size) const {
    return allocate("[unlabeled]", arg_alloc_size);
  }
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
  }

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr,
                  const size_t arg_alloc_size) const {
    deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
  }
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const {
    impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
  }

  /**\brief Return Name of the MemorySpace */
  constexpr static const char* name() { return Namer::get_name(); }

 private:
  BaseSpace underlying_space;
  template <class, class, class, class>
  friend class LogicalMemorySpace;
  friend class Kokkos::Impl::SharedAllocationRecord<memory_space, void>;

  void* impl_allocate(const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      Kokkos::Tools::SpaceHandle arg_handle =
                          Kokkos::Tools::make_space_handle(name())) const {
    return underlying_space.impl_allocate(arg_label, arg_alloc_size,
                                          arg_logical_size, arg_handle);
  }
  void impl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                       const size_t arg_alloc_size,
                       const size_t arg_logical_size = 0,
                       const Kokkos::Tools::SpaceHandle arg_handle =
                           Kokkos::Tools::make_space_handle(name())) const {
    underlying_space.impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size,
                                     arg_logical_size, arg_handle);
  }
};
}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <typename BaseSpace, typename DefaultBaseExecutionSpace, class Namer,
          typename OtherSpace>
struct MemorySpaceAccess<
    Kokkos::Experimental::LogicalMemorySpace<
        BaseSpace, DefaultBaseExecutionSpace, Namer,
        Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>,
    OtherSpace> {
  enum { assignable = MemorySpaceAccess<BaseSpace, OtherSpace>::assignable };
  enum { accessible = MemorySpaceAccess<BaseSpace, OtherSpace>::accessible };
  enum { deepcopy = MemorySpaceAccess<BaseSpace, OtherSpace>::deepcopy };
};

template <typename BaseSpace, typename DefaultBaseExecutionSpace, class Namer,
          typename OtherSpace>
struct MemorySpaceAccess<
    OtherSpace,
    Kokkos::Experimental::LogicalMemorySpace<
        BaseSpace, DefaultBaseExecutionSpace, Namer,
        Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>> {
  enum { assignable = MemorySpaceAccess<OtherSpace, BaseSpace>::assignable };
  enum { accessible = MemorySpaceAccess<OtherSpace, BaseSpace>::accessible };
  enum { deepcopy = MemorySpaceAccess<OtherSpace, BaseSpace>::deepcopy };
};

template <typename BaseSpace, typename DefaultBaseExecutionSpace, class Namer>
struct MemorySpaceAccess<
    Kokkos::Experimental::LogicalMemorySpace<
        BaseSpace, DefaultBaseExecutionSpace, Namer,
        Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>,
    Kokkos::Experimental::LogicalMemorySpace<
        BaseSpace, DefaultBaseExecutionSpace, Namer,
        Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>> {
  enum { assignable = true };
  enum { accessible = true };
  enum { deepcopy = true };
};

}  // namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {
template <class BaseSpace, class DefaultBaseExecutionSpace, class Namer,
          class SharesAccessSemanticsWithBase>
class SharedAllocationRecord<Kokkos::Experimental::LogicalMemorySpace<
                                 BaseSpace, DefaultBaseExecutionSpace, Namer,
                                 SharesAccessSemanticsWithBase>,
                             void> : public SharedAllocationRecord<void, void> {
 private:
  using SpaceType =
      Kokkos::Experimental::LogicalMemorySpace<BaseSpace,
                                               DefaultBaseExecutionSpace, Namer,
                                               SharesAccessSemanticsWithBase>;
  using RecordBase = SharedAllocationRecord<void, void>;

  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  static void deallocate(RecordBase* arg_rec) {
    delete static_cast<SharedAllocationRecord*>(arg_rec);
  }

#ifdef KOKKOS_ENABLE_DEBUG
  /**\brief  Root record for tracked allocations from this
   * LogicalMemorySpace instance */
  static RecordBase s_root_record;
#endif

  const SpaceType m_space;

 protected:
  ~SharedAllocationRecord() {
    m_space.deallocate(RecordBase::m_alloc_ptr->m_label,
                       SharedAllocationRecord<void, void>::m_alloc_ptr,
                       SharedAllocationRecord<void, void>::m_alloc_size,
                       (SharedAllocationRecord<void, void>::m_alloc_size -
                        sizeof(SharedAllocationHeader)));
  }
  SharedAllocationRecord() = default;

  template <typename ExecutionSpace>
  SharedAllocationRecord(
      const ExecutionSpace& /*exec_space*/, const SpaceType& arg_space,
      const std::string& arg_label, const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate)
      : SharedAllocationRecord(arg_space, arg_label, arg_alloc_size,
                               arg_dealloc) {}

  SharedAllocationRecord(
      const SpaceType& arg_space, const std::string& arg_label,
      const size_t arg_alloc_size,
      const RecordBase::function_type arg_dealloc = &deallocate)
      : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
            &SharedAllocationRecord<SpaceType, void>::s_root_record,
#endif
            Impl::checked_allocation_with_header(arg_space, arg_label,
                                                 arg_alloc_size),
            sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
            arg_label),
        m_space(arg_space) {
    // Fill in the Header information
    RecordBase::m_alloc_ptr->m_record =
        static_cast<SharedAllocationRecord<void, void>*>(this);

    strncpy(RecordBase::m_alloc_ptr->m_label, arg_label.c_str(),
            SharedAllocationHeader::maximum_label_length - 1);
    // Set last element zero, in case c_str is too long
    RecordBase::m_alloc_ptr
        ->m_label[SharedAllocationHeader::maximum_label_length - 1] = '\0';
  }

 public:
  inline std::string get_label() const {
    return std::string(RecordBase::head()->m_label);
  }
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const SpaceType& arg_space, const std::string& arg_label,
      const size_t arg_alloc_size) {
    KOKKOS_IF_ON_HOST((return new SharedAllocationRecord(arg_space, arg_label,
                                                         arg_alloc_size);))
    KOKKOS_IF_ON_DEVICE(((void)arg_space; (void)arg_label; (void)arg_alloc_size;
                         return nullptr;))
  }

  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(const SpaceType& arg_space,
                                const std::string& arg_label,
                                const size_t arg_alloc_size) {
    if (!arg_alloc_size) return (void*)nullptr;

    SharedAllocationRecord* const r =
        allocate(arg_space, arg_label, arg_alloc_size);

    RecordBase::increment(r);

    return r->data();
  }

  /**\brief  Reallocate tracked memory in the space */
  static void* reallocate_tracked(void* const arg_alloc_ptr,
                                  const size_t arg_alloc_size) {
    SharedAllocationRecord* const r_old = get_record(arg_alloc_ptr);
    SharedAllocationRecord* const r_new =
        allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

    Kokkos::Impl::DeepCopy<SpaceType, SpaceType>(
        r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));
    Kokkos::fence(
        "SharedAllocationRecord<Kokkos::Experimental::LogicalMemorySpace, "
        "void>::reallocate_tracked: fence after copying data");

    RecordBase::increment(r_new);
    RecordBase::decrement(r_old);

    return r_new->data();
  }
  /**\brief  Deallocate tracked memory in the space */
  static void deallocate_tracked(void* const arg_alloc_ptr) {
    if (arg_alloc_ptr != nullptr) {
      SharedAllocationRecord* const r = get_record(arg_alloc_ptr);

      RecordBase::decrement(r);
    }
  }

  static SharedAllocationRecord* get_record(void* alloc_ptr) {
    using Header     = SharedAllocationHeader;
    using RecordHost = SharedAllocationRecord<SpaceType, void>;

    SharedAllocationHeader const* const head =
        alloc_ptr ? Header::get_header(alloc_ptr)
                  : (SharedAllocationHeader*)nullptr;
    RecordHost* const record =
        head ? static_cast<RecordHost*>(head->m_record) : (RecordHost*)nullptr;

    if (!alloc_ptr || record->m_alloc_ptr != head) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::SharedAllocationRecord< LogicalMemorySpace<> , "
          "void >::get_record ERROR"));
    }

    return record;
  }
#ifdef KOKKOS_ENABLE_DEBUG
  static void print_records(std::ostream& s, const SpaceType&,
                            bool detail = false) {
    SharedAllocationRecord<void, void>::print_host_accessible_records(
        s, "HostSpace", &s_root_record, detail);
  }
#else
  static void print_records(std::ostream&, const SpaceType&,
                            bool detail = false) {
    (void)detail;
    throw_runtime_exception(
        "SharedAllocationRecord<HostSpace>::print_records only works "
        "with KOKKOS_ENABLE_DEBUG enabled");
  }
#endif
};
#ifdef KOKKOS_ENABLE_DEBUG
/**\brief  Root record for tracked allocations from this LogicalSpace
 * instance */
template <class BaseSpace, class DefaultBaseExecutionSpace, class Namer,
          class SharesAccessSemanticsWithBase>
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::Experimental::LogicalMemorySpace<
                               BaseSpace, DefaultBaseExecutionSpace, Namer,
                               SharesAccessSemanticsWithBase>,
                           void>::s_root_record;
#endif

}  // namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <class Namer, class BaseSpace, class DefaultBaseExecutionSpace,
          class SharesAccess, class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::LogicalMemorySpace<
                    BaseSpace, DefaultBaseExecutionSpace, Namer, SharesAccess>,
                Kokkos::Experimental::LogicalMemorySpace<
                    BaseSpace, DefaultBaseExecutionSpace, Namer, SharesAccess>,
                ExecutionSpace> {
  DeepCopy(void* dst, void* src, size_t n) {
    DeepCopy<BaseSpace, BaseSpace, ExecutionSpace>(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, void* src, size_t n) {
    DeepCopy<BaseSpace, BaseSpace, ExecutionSpace>(exec, dst, src, n);
  }
};

template <class Namer, class BaseSpace, class DefaultBaseExecutionSpace,
          class SharesAccess, class ExecutionSpace, class SourceSpace>
struct DeepCopy<SourceSpace,
                Kokkos::Experimental::LogicalMemorySpace<
                    BaseSpace, DefaultBaseExecutionSpace, Namer, SharesAccess>,
                ExecutionSpace> {
  DeepCopy(void* dst, void* src, size_t n) {
    DeepCopy<SourceSpace, BaseSpace, ExecutionSpace>(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, void* src, size_t n) {
    DeepCopy<SourceSpace, BaseSpace, ExecutionSpace>(exec, dst, src, n);
  }
};

template <class Namer, class BaseSpace, class DefaultBaseExecutionSpace,
          class SharesAccess, class ExecutionSpace, class DestinationSpace>
struct DeepCopy<Kokkos::Experimental::LogicalMemorySpace<
                    BaseSpace, DefaultBaseExecutionSpace, Namer, SharesAccess>,
                DestinationSpace, ExecutionSpace> {
  DeepCopy(void* dst, void* src, size_t n) {
    DeepCopy<BaseSpace, DestinationSpace, ExecutionSpace>(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, void* src, size_t n) {
    DeepCopy<BaseSpace, DestinationSpace, ExecutionSpace>(exec, dst, src, n);
  }
};
}  // namespace Impl

}  // namespace Kokkos
#endif  // KOKKOS_LOGICALSPACES_HPP
