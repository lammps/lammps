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

#ifndef KOKKOS_SHARED_ALLOC_HPP
#define KOKKOS_SHARED_ALLOC_HPP

#include <cstdint>
#include <string>

namespace Kokkos {
namespace Impl {

template <class MemorySpace = void, class DestroyFunctor = void>
class SharedAllocationRecord;

class SharedAllocationHeader {
 private:
  typedef SharedAllocationRecord<void, void> Record;

  static constexpr unsigned maximum_label_length =
      (1u << 7 /* 128 */) - sizeof(Record*);

  template <class, class>
  friend class SharedAllocationRecord;

  Record* m_record;
  char m_label[maximum_label_length];

 public:
  /* Given user memory get pointer to the header */
  KOKKOS_INLINE_FUNCTION static const SharedAllocationHeader* get_header(
      void* alloc_ptr) {
    return reinterpret_cast<SharedAllocationHeader*>(
        reinterpret_cast<char*>(alloc_ptr) - sizeof(SharedAllocationHeader));
  }

  KOKKOS_INLINE_FUNCTION
  const char* label() const { return m_label; }
};

template <>
class SharedAllocationRecord<void, void> {
 protected:
  static_assert(sizeof(SharedAllocationHeader) == (1u << 7 /* 128 */),
                "sizeof(SharedAllocationHeader) != 128");

  template <class, class>
  friend class SharedAllocationRecord;

  using function_type = void (*)(SharedAllocationRecord<void, void>*);

  SharedAllocationHeader* const m_alloc_ptr;
  size_t const m_alloc_size;
  function_type const m_dealloc;
#ifdef KOKKOS_DEBUG
  SharedAllocationRecord* const m_root;
  SharedAllocationRecord* m_prev;
  SharedAllocationRecord* m_next;
#endif
  int m_count;

  SharedAllocationRecord(SharedAllocationRecord&&)      = delete;
  SharedAllocationRecord(const SharedAllocationRecord&) = delete;
  SharedAllocationRecord& operator=(SharedAllocationRecord&&) = delete;
  SharedAllocationRecord& operator=(const SharedAllocationRecord&) = delete;

  /**\brief  Construct and insert into 'arg_root' tracking set.
   *         use_count is zero.
   */
  SharedAllocationRecord(
#ifdef KOKKOS_DEBUG
      SharedAllocationRecord* arg_root,
#endif
      SharedAllocationHeader* arg_alloc_ptr, size_t arg_alloc_size,
      function_type arg_dealloc);
 private:
  static __thread int t_tracking_enabled;

 public:
  virtual std::string get_label() const { return std::string("Unmanaged"); }

  static int tracking_enabled() { return t_tracking_enabled; }

  /**\brief A host process thread claims and disables the
   *        shared allocation tracking flag.
   */
  static void tracking_disable() { t_tracking_enabled = 0; }

  /**\brief A host process thread releases and enables the
   *        shared allocation tracking flag.
   */
  static void tracking_enable() { t_tracking_enabled = 1; }

  virtual ~SharedAllocationRecord() {}

  SharedAllocationRecord()
      : m_alloc_ptr(nullptr),
        m_alloc_size(0),
        m_dealloc(nullptr)
#ifdef KOKKOS_DEBUG
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
  void* data() const { return reinterpret_cast<void*>(m_alloc_ptr + 1); }

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

namespace {

/* Taking the address of this function so make sure it is unique */
template <class MemorySpace, class DestroyFunctor>
void deallocate(SharedAllocationRecord<void, void>* record_ptr) {
  typedef SharedAllocationRecord<MemorySpace, void> base_type;
  typedef SharedAllocationRecord<MemorySpace, DestroyFunctor> this_type;

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
  // one to zero removes from the trakcing list and deallocates.
  KOKKOS_INLINE_FUNCTION static SharedAllocationRecord* allocate(
      const MemorySpace& arg_space, const std::string& arg_label,
      const size_t arg_alloc) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    return new SharedAllocationRecord(arg_space, arg_label, arg_alloc);
#else
    return (SharedAllocationRecord*)0;
#endif
  }
};

template <class MemorySpace>
class SharedAllocationRecord<MemorySpace, void>
    : public SharedAllocationRecord<void, void> {};

union SharedAllocationTracker {
 private:
  typedef SharedAllocationRecord<void, void> Record;

  enum : uintptr_t { DO_NOT_DEREF_FLAG = 0x01ul };

  // The allocation record resides in Host memory space
  uintptr_t m_record_bits;
  Record* m_record;

 public:
  // Use macros instead of inline functions to reduce
  // pressure on compiler optimization by reducing
  // number of symbols and inline functons.

#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_ENABLED Record::tracking_enabled()

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT \
  if (!(m_record_bits & DO_NOT_DEREF_FLAG)) Record::increment(m_record);

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT \
  if (!(m_record_bits & DO_NOT_DEREF_FLAG)) Record::decrement(m_record);

#else

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_ENABLED 0

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT /* */

#define KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT /* */

#endif

#define KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(rhs,               \
                                                        override_tracking) \
  (((!override_tracking) || (rhs.m_record_bits & DO_NOT_DEREF_FLAG) ||     \
    (!KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_ENABLED))                      \
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
               ? (SharedAllocationRecord<MemorySpace, void>*)0
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
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    Record* const tmp =
        reinterpret_cast<Record*>(m_record_bits & ~DO_NOT_DEREF_FLAG);
    return (tmp ? tmp->use_count() : 0);
#else
    return 0;
#endif
  }

  KOKKOS_INLINE_FUNCTION
  bool has_record() const {
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

  /** \brief  Copy assignment may disable tracking */
  KOKKOS_FORCEINLINE_FUNCTION
  void assign(const SharedAllocationTracker& rhs, const bool enable_tracking) {
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
    m_record_bits =
        KOKKOS_IMPL_SHARED_ALLOCATION_CARRY_RECORD_BITS(rhs, enable_tracking);
    KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
  }

#undef KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_ENABLED
#undef KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_INCREMENT
#undef KOKKOS_IMPL_SHARED_ALLOCATION_TRACKER_DECREMENT
};

} /* namespace Impl */
} /* namespace Kokkos */

#endif
