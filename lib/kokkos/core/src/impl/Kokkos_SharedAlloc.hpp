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

#ifndef KOKKOS_SHARED_ALLOC_HPP
#define KOKKOS_SHARED_ALLOC_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Error.hpp>  // Impl::throw_runtime_exception

#include <cstdint>
#include <string>

namespace Kokkos {
namespace Impl {

template <class MemorySpace = void, class DestroyFunctor = void>
class SharedAllocationRecord;

template <class MemorySpace>
class SharedAllocationRecordCommon;

class SharedAllocationHeader {
 private:
  using Record = SharedAllocationRecord<void, void>;

#if defined(KOKKOS_ARCH_VEGA) || defined(KOKKOS_ARCH_NAVI)
  static constexpr unsigned maximum_label_length =
      (1u << 8 /* 256 */) - sizeof(Record*);
#else
  static constexpr unsigned maximum_label_length =
      (1u << 7 /* 128 */) - sizeof(Record*);
#endif

  template <class, class>
  friend class SharedAllocationRecord;
  template <class>
  friend class SharedAllocationRecordCommon;
  template <class>
  friend class HostInaccessibleSharedAllocationRecordCommon;

  Record* m_record;
  char m_label[maximum_label_length];

 public:
  /* Given user memory get pointer to the header */
  KOKKOS_INLINE_FUNCTION static const SharedAllocationHeader* get_header(
      void const* alloc_ptr) {
    return reinterpret_cast<SharedAllocationHeader const*>(
        static_cast<char const*>(alloc_ptr) - sizeof(SharedAllocationHeader));
  }

  KOKKOS_INLINE_FUNCTION
  const char* label() const { return m_label; }
};

template <>
class SharedAllocationRecord<void, void> {
 protected:
#if defined(KOKKOS_ARCH_VEGA) || defined(KOKKOS_ARCH_NAVI)
  static_assert(sizeof(SharedAllocationHeader) == (1u << 8 /* 256 */),
                "sizeof(SharedAllocationHeader) != 256");
#else
  static_assert(sizeof(SharedAllocationHeader) == (1u << 7 /* 128 */),
                "sizeof(SharedAllocationHeader) != 128");
#endif

  template <class, class>
  friend class SharedAllocationRecord;
  template <class>
  friend class SharedAllocationRecordCommon;
  template <class>
  friend class HostInaccessibleSharedAllocationRecordCommon;

  using function_type = void (*)(SharedAllocationRecord<void, void>*);

  SharedAllocationHeader* const m_alloc_ptr;
  size_t const m_alloc_size;
  function_type const m_dealloc;
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord* const m_root;
  SharedAllocationRecord* m_prev;
  SharedAllocationRecord* m_next;
#endif
  int m_count;
  std::string m_label;

  SharedAllocationRecord(SharedAllocationRecord&&)      = delete;
  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(SharedAllocationRecord&&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  /**\brief  Construct and insert into 'arg_root' tracking set.
   *         use_count is zero.
   */
  SharedAllocationRecord(
#ifdef KOKKOS_ENABLE_DEBUG
      SharedAllocationRecord* arg_root,
#endif
      SharedAllocationHeader* arg_alloc_ptr, size_t arg_alloc_size,
      function_type arg_dealloc, const std::string& label);
 private:
  static inline thread_local int t_tracking_enabled = 1;

 public:
  virtual std::string get_label() const { return std::string("Unmanaged"); }

#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
  static KOKKOS_FUNCTION int tracking_enabled() {
    KOKKOS_IF_ON_HOST(return t_tracking_enabled;)
    KOKKOS_IF_ON_DEVICE(return 0;)
  }
#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma pop
#endif

  /**\brief A host process thread claims and disables the
   *        shared allocation tracking flag.
   */
  static void tracking_disable() { t_tracking_enabled = 0; }

  /**\brief A host process thread releases and enables the
   *        shared allocation tracking flag.
   */
  static void tracking_enable() { t_tracking_enabled = 1; }

  virtual ~SharedAllocationRecord() = default;

  SharedAllocationRecord()
      : m_alloc_ptr(nullptr),
        m_alloc_size(0),
        m_dealloc(nullptr)
#ifdef KOKKOS_ENABLE_DEBUG
        ,
        m_root(this),
        m_prev(this),
        m_next(this)
#endif
        ,
        m_count(0) {
  }

  static constexpr unsigned maximum_label_length =
      SharedAllocationHeader::maximum_label_length;

  KOKKOS_INLINE_FUNCTION
  const SharedAllocationHeader* head() const { return m_alloc_ptr; }

  /* User's memory begins at the end of the header */
  KOKKOS_INLINE_FUNCTION
  void* data() const { return static_cast<void*>(m_alloc_ptr + 1); }

  /* User's memory begins at the end of the header */
  size_t size() const { return m_alloc_size - sizeof(SharedAllocationHeader); }

  /* Cannot be 'constexpr' because 'm_count' is volatile */
  int use_count() const { return *static_cast<const volatile int*>(&m_count); }

  /* Increment use count */
  static void increment(SharedAllocationRecord*);

  /* Decrement use count. If 1->0 then remove from the tracking list and invoke
   * m_dealloc */
  static SharedAllocationRecord* decrement(SharedAllocationRecord*);

  /* Given a root record and data pointer find the record */
  static SharedAllocationRecord* find(SharedAllocationRecord* const,
                                      void* const);

  /*  Sanity check for the whole set of records to which the input record
   * belongs. Locks the set's insert/erase operations until the sanity check is
   * complete.
   */
  static bool is_sane(SharedAllocationRecord*);

  /*  Print host-accessible records */
  static void print_host_accessible_records(
      std::ostream&, const char* const space_name,
      const SharedAllocationRecord* const root, const bool detail);
};

template <class MemorySpace>
class SharedAllocationRecordCommon : public SharedAllocationRecord<void, void> {
 private:
  using derived_t     = SharedAllocationRecord<MemorySpace, void>;
  using record_base_t = SharedAllocationRecord<void, void>;
  derived_t& self() { return *static_cast<derived_t*>(this); }
  derived_t const& self() const { return *static_cast<derived_t const*>(this); }

 protected:
  using record_base_t::record_base_t;

  void _fill_host_accessible_header_info(SharedAllocationHeader& arg_header,
                                         std::string const& arg_label);

  static void deallocate(record_base_t* arg_rec);

 public:
  static auto allocate(MemorySpace const& arg_space,
                       std::string const& arg_label, size_t arg_alloc_size)
      -> derived_t*;
  /**\brief  Allocate tracked memory in the space */
  static void* allocate_tracked(MemorySpace const& arg_space,
                                std::string const& arg_alloc_label,
                                size_t arg_alloc_size);
  /**\brief  Reallocate tracked memory in the space */
  static void deallocate_tracked(void* arg_alloc_ptr);
  /**\brief  Deallocate tracked memory in the space */
  static void* reallocate_tracked(void* arg_alloc_ptr, size_t arg_alloc_size);
  static auto get_record(void* alloc_ptr) -> derived_t*;
  std::string get_label() const;
  static void print_records(std::ostream& s, MemorySpace const&,
                            bool detail = false);
};

template <class MemorySpace>
class HostInaccessibleSharedAllocationRecordCommon
    : public SharedAllocationRecordCommon<MemorySpace> {
 private:
  using base_t        = SharedAllocationRecordCommon<MemorySpace>;
  using derived_t     = SharedAllocationRecord<MemorySpace, void>;
  using record_base_t = SharedAllocationRecord<void, void>;

 protected:
  using base_t::base_t;

 public:
  static void print_records(std::ostream& s, MemorySpace const&,
                            bool detail = false);
  static auto get_record(void* alloc_ptr) -> derived_t*;
  std::string get_label() const;
};

namespace {

/* Taking the address of this function so make sure it is unique */
template <class MemorySpace, class DestroyFunctor>
void deallocate(SharedAllocationRecord<void, void>* record_ptr) {
  using base_type = SharedAllocationRecord<MemorySpace, void>;
  using this_type = SharedAllocationRecord<MemorySpace, DestroyFunctor>;

  this_type* const ptr =
      static_cast<this_type*>(static_cast<base_type*>(record_ptr));

  ptr->m_destroy.destroy_shared_allocation();

  delete ptr;
}

}  // namespace

/*
 *  Memory space specialization of SharedAllocationRecord< Space , void >
 * requires :
 *
 *  SharedAllocationRecord< Space , void > : public SharedAllocationRecord< void
 * , void >
 *  {
 *    // delete allocated user memory via static_cast to this type.
 *    static void deallocate( const SharedAllocationRecord<void,void> * );
 *    Space m_space ;
 *  }
 */
template <class MemorySpace, class DestroyFunctor>
class SharedAllocationRecord
    : public SharedAllocationRecord<MemorySpace, void> {
 private:
  template <typename ExecutionSpace>
  SharedAllocationRecord(const ExecutionSpace& execution_space,
                         const MemorySpace& arg_space,
                         const std::string& arg_label, const size_t arg_alloc)
      /*  Allocate user memory as [ SharedAllocationHeader , user_memory ] */
      : SharedAllocationRecord<MemorySpace, void>(
            execution_space, arg_space, arg_label, arg_alloc,
            &Kokkos::Impl::deallocate<MemorySpace, DestroyFunctor>),
        m_destroy() {}

  SharedAllocationRecord(const MemorySpace& arg_space,
                         const std::string& arg_label, const size_t arg_alloc)
      /*  Allocate user memory as [ SharedAllocationHeader , user_memory ] */
      : SharedAllocationRecord<MemorySpace, void>(
            arg_space, arg_label, arg_alloc,
            &Kokkos::Impl::deallocate<MemorySpace, DestroyFunctor>),
        m_destroy() {}

  SharedAllocationRecord()                              = delete;
  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

 public:
  DestroyFunctor m_destroy;

  // Allocate with a zero use count.  Incrementing the use count from zero to
  // one inserts the record into the tracking list.  Decrementing the count from
  // one to zero removes from the tracking list and deallocates.
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const MemorySpace& arg_space, const std::string& arg_label,
      const size_t arg_alloc) {
    KOKKOS_IF_ON_HOST(
        (return new SharedAllocationRecord(arg_space, arg_label, arg_alloc);))
    KOKKOS_IF_ON_DEVICE(
        ((void)arg_space; (void)arg_label; (void)arg_alloc; return nullptr;))
  }

  template <typename ExecutionSpace>
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const ExecutionSpace& exec_space, const MemorySpace& arg_space,
      const std::string& arg_label, const size_t arg_alloc) {
    KOKKOS_IF_ON_HOST(
        (return new SharedAllocationRecord(exec_space, arg_space, arg_label,
                                           arg_alloc);))
    KOKKOS_IF_ON_DEVICE(((void)exec_space; (void)arg_space; (void)arg_label;
                         (void)arg_alloc; return nullptr;))
  }
};

template <class MemorySpace>
class SharedAllocationRecord<MemorySpace, void>
    : public SharedAllocationRecord<void, void> {};

union SharedAllocationTracker {
 private:
  using Record = SharedAllocationRecord<void, void>;

  enum : uintptr_t { DO_NOT_DEREF_FLAG = 0x01ul };

  // The allocation record resides in Host memory space
  uintptr_t m_record_bits;
  Record* m_record;

 public:
  // Use macros instead of inline functions to reduce
  // pressure on compiler optimization by reducing
  // number of symbols and inline functions.

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT          \
  KOKKOS_IF_ON_HOST((if (!(m_record_bits & DO_NOT_DEREF_FLAG)) { \
    Record::increment(m_record);                                 \
  }))

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT          \
  KOKKOS_IF_ON_HOST((if (!(m_record_bits & DO_NOT_DEREF_FLAG)) { \
    Record::decrement(m_record);                                 \
  }))

#define KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(rhs,               \
                                                        override_tracking) \
  (((!override_tracking) || (rhs.m_record_bits & DO_NOT_DEREF_FLAG) ||     \
    (!Record::tracking_enabled()))                                         \
       ? rhs.m_record_bits | DO_NOT_DEREF_FLAG                             \
       : rhs.m_record_bits)

  /** \brief  Assign a specialized record */
  inline void assign_allocated_record_to_uninitialized(Record* arg_record) {
    if (arg_record) {
      Record::increment(m_record = arg_record);
    } else {
      m_record_bits = DO_NOT_DEREF_FLAG;
    }
  }

  template <class MemorySpace>
  constexpr SharedAllocationRecord<MemorySpace, void>* get_record() const
      noexcept {
    return (m_record_bits & DO_NOT_DEREF_FLAG)
               ? nullptr
               : static_cast<SharedAllocationRecord<MemorySpace, void>*>(
                     m_record);
  }

  template <class MemorySpace>
  std::string get_label() const {
    return (m_record_bits == DO_NOT_DEREF_FLAG)
               ? std::string()
               : reinterpret_cast<SharedAllocationRecord<MemorySpace, void>*>(
                     m_record_bits & ~DO_NOT_DEREF_FLAG)
                     ->get_label();
  }

  KOKKOS_INLINE_FUNCTION
  int use_count() const {
    KOKKOS_IF_ON_HOST((Record* const tmp = reinterpret_cast<Record*>(
                           m_record_bits & ~DO_NOT_DEREF_FLAG);
                       return (tmp ? tmp->use_count() : 0);))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  KOKKOS_INLINE_FUNCTION bool has_record() const {
    return (m_record_bits & (~DO_NOT_DEREF_FLAG)) != 0;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void clear() {
    // If this is tracking then must decrement
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    // Reset to default constructed value.
    m_record_bits = DO_NOT_DEREF_FLAG;
  }

  // Copy:
  KOKKOS_FORCEINLINE_FUNCTION
  ~SharedAllocationTracker(){KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT}

  KOKKOS_FORCEINLINE_FUNCTION constexpr SharedAllocationTracker()
      : m_record_bits(DO_NOT_DEREF_FLAG) {}

  // Move:

  KOKKOS_FORCEINLINE_FUNCTION
  SharedAllocationTracker(SharedAllocationTracker&& rhs)
      : m_record_bits(rhs.m_record_bits) {
    rhs.m_record_bits = DO_NOT_DEREF_FLAG;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  SharedAllocationTracker& operator=(SharedAllocationTracker&& rhs) {
    auto swap_tmp     = m_record_bits;
    m_record_bits     = rhs.m_record_bits;
    rhs.m_record_bits = swap_tmp;
    return *this;
  }

  // Copy:

  KOKKOS_FORCEINLINE_FUNCTION
  SharedAllocationTracker(const SharedAllocationTracker& rhs)
      : m_record_bits(KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(
            rhs, true)){KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT}

        /** \brief  Copy construction may disable tracking. */
        KOKKOS_FORCEINLINE_FUNCTION SharedAllocationTracker(
            const SharedAllocationTracker& rhs, const bool enable_tracking)
      : m_record_bits(KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(
            rhs,
            enable_tracking)){KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT}

        KOKKOS_FORCEINLINE_FUNCTION SharedAllocationTracker
        &
        operator=(const SharedAllocationTracker& rhs) {
    // If this is tracking then must decrement
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    m_record_bits = KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(rhs, true);
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
    return *this;
  }

  /*  The following functions (assign_direct and assign_force_disable)
   *  are the result of deconstructing the
   *  KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS macro.  This
   *  allows the caller to do the check for tracking enabled and managed
   *  apart from the assignment of the record because the tracking
   *  enabled / managed question may be important for other tasks as well
   */

  /** \brief  Copy assignment without the carry bits logic
   *         This assumes that externally defined tracking is explicitly enabled
   */
  KOKKOS_FORCEINLINE_FUNCTION
  void assign_direct(const SharedAllocationTracker& rhs) {
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    m_record_bits = rhs.m_record_bits;
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
  }

  /** \brief  Copy assignment without the increment
   *         we cannot assume that current record is unmanaged
   *         but with externally defined tracking explicitly disabled
   *         we can go straight to the do not deref flag     */
  KOKKOS_FORCEINLINE_FUNCTION
  void assign_force_disable(const SharedAllocationTracker& rhs) {
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    m_record_bits = rhs.m_record_bits | DO_NOT_DEREF_FLAG;
  }

  // report if record is tracking or not
  KOKKOS_FORCEINLINE_FUNCTION
  bool tracking_enabled() { return (!(m_record_bits & DO_NOT_DEREF_FLAG)); }

  /** \brief  Copy assignment may disable tracking */
  KOKKOS_FORCEINLINE_FUNCTION
  void assign(const SharedAllocationTracker& rhs, const bool enable_tracking) {
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    m_record_bits =
        KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(rhs, enable_tracking);
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
  }

#undef KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
#undef KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
};

} /* namespace Impl */
} /* namespace Kokkos */
#endif
