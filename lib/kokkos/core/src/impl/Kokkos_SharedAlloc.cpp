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

#include <Kokkos_Core.hpp>

namespace Kokkos {
namespace Impl {

KOKKOS_THREAD_LOCAL int SharedAllocationRecord<void, void>::t_tracking_enabled =
    1;

#ifdef KOKKOS_DEBUG
bool SharedAllocationRecord<void, void>::is_sane(
    SharedAllocationRecord<void, void>* arg_record) {
  SharedAllocationRecord* const root = arg_record ? arg_record->m_root : 0;

  bool ok = root != 0 && root->use_count() == 0;

  if (ok) {
    SharedAllocationRecord* root_next             = 0;
    static constexpr SharedAllocationRecord* zero = nullptr;
    // Lock the list:
    while ((root_next = Kokkos::atomic_exchange(&root->m_next, zero)) ==
           nullptr)
      ;

    for (SharedAllocationRecord* rec = root_next; ok && rec != root;
         rec                         = rec->m_next) {
      const bool ok_non_null =
          rec && rec->m_prev && (rec == root || rec->m_next);
      const bool ok_root = ok_non_null && rec->m_root == root;
      const bool ok_prev_next =
          ok_non_null &&
          (rec->m_prev != root ? rec->m_prev->m_next == rec : root_next == rec);
      const bool ok_next_prev = ok_non_null && rec->m_next->m_prev == rec;
      const bool ok_count     = ok_non_null && 0 <= rec->use_count();

      ok = ok_root && ok_prev_next && ok_next_prev && ok_count;

      if (!ok) {
        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string =
              "Kokkos::Impl::SharedAllocationRecord failed is_sane: "
              "rec(0x%.12lx){ m_count(%d) m_root(0x%.12lx) m_next(0x%.12lx) "
              "m_prev(0x%.12lx) m_next->m_prev(0x%.12lx) "
              "m_prev->m_next(0x%.12lx) }\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string =
              "Kokkos::Impl::SharedAllocationRecord failed is_sane: "
              "rec(0x%.12llx){ m_count(%d) m_root(0x%.12llx) m_next(0x%.12llx) "
              "m_prev(0x%.12llx) m_next->m_prev(0x%.12llx) "
              "m_prev->m_next(0x%.12llx) }\n";
        }

        fprintf(stderr, format_string, reinterpret_cast<uintptr_t>(rec),
                rec->use_count(), reinterpret_cast<uintptr_t>(rec->m_root),
                reinterpret_cast<uintptr_t>(rec->m_next),
                reinterpret_cast<uintptr_t>(rec->m_prev),
                reinterpret_cast<uintptr_t>(
                    rec->m_next != nullptr ? rec->m_next->m_prev : nullptr),
                reinterpret_cast<uintptr_t>(rec->m_prev != rec->m_root
                                                ? rec->m_prev->m_next
                                                : root_next));
      }
    }

    if (nullptr != Kokkos::atomic_exchange(&root->m_next, root_next)) {
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::Impl::SharedAllocationRecord failed is_sane unlocking");
    }
  }
  return ok;
}

#else

bool SharedAllocationRecord<void, void>::is_sane(
    SharedAllocationRecord<void, void>*) {
  Kokkos::Impl::throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord::is_sane only works with "
      "KOKKOS_DEBUG enabled");
  return false;
}
#endif  //#ifdef KOKKOS_DEBUG

#ifdef KOKKOS_DEBUG
SharedAllocationRecord<void, void>* SharedAllocationRecord<void, void>::find(
    SharedAllocationRecord<void, void>* const arg_root,
    void* const arg_data_ptr) {
  SharedAllocationRecord* root_next             = 0;
  static constexpr SharedAllocationRecord* zero = nullptr;

  // Lock the list:
  while ((root_next = Kokkos::atomic_exchange(&arg_root->m_next, zero)) ==
         nullptr)
    ;

  // Iterate searching for the record with this data pointer

  SharedAllocationRecord* r = root_next;

  while ((r != arg_root) && (r->data() != arg_data_ptr)) {
    r = r->m_next;
  }

  if (r == arg_root) {
    r = 0;
  }

  if (nullptr != Kokkos::atomic_exchange(&arg_root->m_next, root_next)) {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::SharedAllocationRecord failed locking/unlocking");
  }
  return r;
}
#else
SharedAllocationRecord<void, void>* SharedAllocationRecord<void, void>::find(
    SharedAllocationRecord<void, void>* const, void* const) {
  Kokkos::Impl::throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord::find only works with KOKKOS_DEBUG "
      "enabled");
  return nullptr;
}
#endif

/**\brief  Construct and insert into 'arg_root' tracking set.
 *         use_count is zero.
 */
SharedAllocationRecord<void, void>::SharedAllocationRecord(
#ifdef KOKKOS_DEBUG
    SharedAllocationRecord<void, void>* arg_root,
#endif
    SharedAllocationHeader* arg_alloc_ptr, size_t arg_alloc_size,
    SharedAllocationRecord<void, void>::function_type arg_dealloc)
    : m_alloc_ptr(arg_alloc_ptr),
      m_alloc_size(arg_alloc_size),
      m_dealloc(arg_dealloc)
#ifdef KOKKOS_DEBUG
      ,
      m_root(arg_root),
      m_prev(0),
      m_next(0)
#endif
      ,
      m_count(0) {
  if (nullptr != arg_alloc_ptr) {
#ifdef KOKKOS_DEBUG
    // Insert into the root double-linked list for tracking
    //
    // before:  arg_root->m_next == next ; next->m_prev == arg_root
    // after:   arg_root->m_next == this ; this->m_prev == arg_root ;
    //              this->m_next == next ; next->m_prev == this

    m_prev                                        = m_root;
    static constexpr SharedAllocationRecord* zero = nullptr;

    // Read root->m_next and lock by setting to nullptr
    while ((m_next = Kokkos::atomic_exchange(&m_root->m_next, zero)) == nullptr)
      ;

    m_next->m_prev = this;

    // memory fence before completing insertion into linked list
    Kokkos::memory_fence();

    if (nullptr != Kokkos::atomic_exchange(&m_root->m_next, this)) {
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::Impl::SharedAllocationRecord failed locking/unlocking");
    }
#endif

  } else {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::SharedAllocationRecord given nullptr allocation");
  }
}

void SharedAllocationRecord<void, void>::increment(
    SharedAllocationRecord<void, void>* arg_record) {
  const int old_count = Kokkos::atomic_fetch_add(&arg_record->m_count, 1);

  if (old_count < 0) {  // Error
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::SharedAllocationRecord failed increment");
  }
}

SharedAllocationRecord<void, void>* SharedAllocationRecord<
    void, void>::decrement(SharedAllocationRecord<void, void>* arg_record) {
  const int old_count = Kokkos::atomic_fetch_sub(&arg_record->m_count, 1);

  if (old_count == 1) {
    if (!Kokkos::is_initialized()) {
      std::stringstream ss;
      ss << "Kokkos allocation \"";
      ss << arg_record->get_label();
      ss << "\" is being deallocated after Kokkos::finalize was called\n";
      auto s = ss.str();
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      std::cerr << s;
      std::cerr << "This behavior is incorrect Kokkos usage, and will crash in "
                   "future releases\n";
#else
      Kokkos::Impl::throw_runtime_exception(s);
#endif
    }

#ifdef KOKKOS_DEBUG
    // before:  arg_record->m_prev->m_next == arg_record  &&
    //          arg_record->m_next->m_prev == arg_record
    //
    // after:   arg_record->m_prev->m_next == arg_record->m_next  &&
    //          arg_record->m_next->m_prev == arg_record->m_prev

    SharedAllocationRecord* root_next             = 0;
    static constexpr SharedAllocationRecord* zero = nullptr;

    // Lock the list:
    while ((root_next = Kokkos::atomic_exchange(&arg_record->m_root->m_next,
                                                zero)) == nullptr)
      ;

    arg_record->m_next->m_prev = arg_record->m_prev;

    if (root_next != arg_record) {
      arg_record->m_prev->m_next = arg_record->m_next;
    } else {
      // before:  arg_record->m_root == arg_record->m_prev
      // after:   arg_record->m_root == arg_record->m_next
      root_next = arg_record->m_next;
    }

    Kokkos::memory_fence();

    // Unlock the list:
    if (nullptr !=
        Kokkos::atomic_exchange(&arg_record->m_root->m_next, root_next)) {
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::Impl::SharedAllocationRecord failed decrement unlocking");
    }

    arg_record->m_next = 0;
    arg_record->m_prev = 0;
#endif

    function_type d = arg_record->m_dealloc;
    (*d)(arg_record);
    arg_record = nullptr;
  } else if (old_count < 1) {  // Error
    fprintf(stderr,
            "Kokkos::Impl::SharedAllocationRecord '%s' failed decrement count "
            "= %d\n",
            arg_record->m_alloc_ptr->m_label, old_count);
    fflush(stderr);
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::SharedAllocationRecord failed decrement count");
  }

  return arg_record;
}

#ifdef KOKKOS_DEBUG
void SharedAllocationRecord<void, void>::print_host_accessible_records(
    std::ostream& s, const char* const space_name,
    const SharedAllocationRecord* const root, const bool detail) {
  const SharedAllocationRecord<void, void>* r = root;

  char buffer[256];

  if (detail) {
    do {
      // Formatting dependent on sizeof(uintptr_t)
      const char* format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) {
        format_string =
            "%s addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx + "
            "%.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
        format_string =
            "%s addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ "
            "0x%.12llx + %.8ld ] count(%d) dealloc(0x%.12llx) %s\n";
      }

      snprintf(buffer, 256, format_string, space_name,
               reinterpret_cast<uintptr_t>(r),
               reinterpret_cast<uintptr_t>(r->m_prev),
               reinterpret_cast<uintptr_t>(r->m_next),
               reinterpret_cast<uintptr_t>(r->m_alloc_ptr), r->m_alloc_size,
               r->use_count(), reinterpret_cast<uintptr_t>(r->m_dealloc),
               r->m_alloc_ptr->m_label);
      s << buffer;
      r = r->m_next;
    } while (r != root);
  } else {
    do {
      if (r->m_alloc_ptr) {
        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string = "%s [ 0x%.12lx + %ld ] %s\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string = "%s [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf(buffer, 256, format_string, space_name,
                 reinterpret_cast<uintptr_t>(r->data()), r->size(),
                 r->m_alloc_ptr->m_label);
      } else {
        snprintf(buffer, 256, "%s [ 0 + 0 ]\n", space_name);
      }
      s << buffer;
      r = r->m_next;
    } while (r != root);
  }
}
#else
void SharedAllocationRecord<void, void>::print_host_accessible_records(
    std::ostream&, const char* const, const SharedAllocationRecord* const,
    const bool) {
  Kokkos::Impl::throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord::print_host_accessible_records"
      " only works with KOKKOS_DEBUG enabled");
}
#endif

} /* namespace Impl */
} /* namespace Kokkos */
