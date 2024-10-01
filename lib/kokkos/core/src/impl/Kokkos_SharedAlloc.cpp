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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
bool SharedAllocationRecord<void, void>::is_sane(
    SharedAllocationRecord<void, void>* arg_record) {
  SharedAllocationRecord* const root =
      arg_record ? arg_record->m_root : nullptr;

  bool ok = root != nullptr && root->use_count() == 0;

  if (ok) {
    SharedAllocationRecord* root_next             = nullptr;
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
      "KOKKOS_ENABLE_DEBUG enabled");
  return false;
}
#endif  //#ifdef KOKKOS_ENABLE_DEBUG

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>* SharedAllocationRecord<void, void>::find(
    SharedAllocationRecord<void, void>* const arg_root,
    void* const arg_data_ptr) {
  SharedAllocationRecord* root_next             = nullptr;
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
    r = nullptr;
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
      "Kokkos::Impl::SharedAllocationRecord::find only works with "
      "KOKKOS_ENABLE_DEBUG "
      "enabled");
  return nullptr;
}
#endif

/**\brief  Construct and insert into 'arg_root' tracking set.
 *         use_count is zero.
 */
SharedAllocationRecord<void, void>::SharedAllocationRecord(
#ifdef KOKKOS_ENABLE_DEBUG
    SharedAllocationRecord<void, void>* arg_root,
#endif
    SharedAllocationHeader* arg_alloc_ptr, size_t arg_alloc_size,
    SharedAllocationRecord<void, void>::function_type arg_dealloc,
    const std::string& label)
    : m_alloc_ptr(arg_alloc_ptr),
      m_alloc_size(arg_alloc_size),
      m_dealloc(arg_dealloc)
#ifdef KOKKOS_ENABLE_DEBUG
      ,
      m_root(arg_root),
      m_prev(nullptr),
      m_next(nullptr)
#endif
      ,
      m_count(0),
      m_label(label) {
  if (nullptr != arg_alloc_ptr) {
#ifdef KOKKOS_ENABLE_DEBUG
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
    if (is_finalized()) {
      std::stringstream ss;
      ss << "Kokkos allocation \"";
      ss << arg_record->get_label();
      ss << "\" is being deallocated after Kokkos::finalize was called\n";
      auto s = ss.str();
      Kokkos::Impl::throw_runtime_exception(s);
    }

#ifdef KOKKOS_ENABLE_DEBUG
    // before:  arg_record->m_prev->m_next == arg_record  &&
    //          arg_record->m_next->m_prev == arg_record
    //
    // after:   arg_record->m_prev->m_next == arg_record->m_next  &&
    //          arg_record->m_next->m_prev == arg_record->m_prev

    SharedAllocationRecord* root_next             = nullptr;
    static constexpr SharedAllocationRecord* zero = nullptr;

    // Lock the list:
    while ((root_next = Kokkos::atomic_exchange(&arg_record->m_root->m_next,
                                                zero)) == nullptr)
      ;
    // We need a memory_fence() here so that the following update
    // is properly sequenced
    Kokkos::memory_fence();

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

    arg_record->m_next = nullptr;
    arg_record->m_prev = nullptr;
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

#ifdef KOKKOS_ENABLE_DEBUG
void SharedAllocationRecord<void, void>::print_host_accessible_records(
    std::ostream& s, const char* const space_name,
    const SharedAllocationRecord* const root, const bool detail) {
  // Print every node except the root, which does not represent an actual
  // allocation.
  const SharedAllocationRecord<void, void>* r = root->m_next;

  std::ios_base::fmtflags saved_flags = s.flags();
#define KOKKOS_PAD_HEX(ptr)                              \
  "0x" << std::hex << std::setw(12) << std::setfill('0') \
       << reinterpret_cast<uintptr_t>(ptr)
  if (detail) {
    while (r != root) {
      s << space_name << " addr( " << KOKKOS_PAD_HEX(r) << " ) list ( "
        << KOKKOS_PAD_HEX(r->m_prev) << ' ' << KOKKOS_PAD_HEX(r->m_next)
        << " ) extent[ " << KOKKOS_PAD_HEX(r->m_alloc_ptr) << " + " << std::dec
        << std::setw(8) << r->m_alloc_size << " ] count(" << r->use_count()
        << ") dealloc(" << KOKKOS_PAD_HEX(r->m_dealloc) << ") "
        << r->m_alloc_ptr->m_label << '\n';

      r = r->m_next;
    }
  } else {
    while (r != root) {
      s << space_name << " [ " << KOKKOS_PAD_HEX(r->data()) << " + " << std::dec
        << r->size() << " ] " << r->m_alloc_ptr->m_label << '\n';
      r = r->m_next;
    }
  }
#undef KOKKOS_PAD_HEX
  s.flags(saved_flags);
}
#else
void SharedAllocationRecord<void, void>::print_host_accessible_records(
    std::ostream&, const char* const, const SharedAllocationRecord* const,
    const bool) {
  Kokkos::Impl::throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord::print_host_accessible_records"
      " only works with KOKKOS_ENABLE_DEBUG enabled");
}
#endif

void fill_host_accessible_header_info(
    SharedAllocationRecord<void, void>* arg_record,
    SharedAllocationHeader& arg_header, std::string const& arg_label) {
  // Fill in the Header information, directly accessible on the host

  arg_header.m_record = arg_record;

  strncpy(arg_header.m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  arg_header.m_label[SharedAllocationHeader::maximum_label_length - 1] = '\0';
}

} /* namespace Impl */
} /* namespace Kokkos */
