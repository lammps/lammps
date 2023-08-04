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

/** @file Kokkos_MemorySpace.cpp
 *
 *  Operations common to memory space instances, or at least default
 *  implementations thereof.
 */

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <impl/Kokkos_MemorySpace.hpp>

#include <iostream>
#include <string>
#include <sstream>

namespace Kokkos {
namespace Impl {

void safe_throw_allocation_with_header_failure(
    std::string const& space_name, std::string const& label,
    Kokkos::Experimental::RawMemoryAllocationFailure const& failure) {
  auto generate_failure_message = [&](std::ostream& o) {
    o << "Kokkos failed to allocate memory for label \"" << label
      << "\".  Allocation using MemorySpace named \"" << space_name
      << "\" failed with the following error:  ";
    failure.print_error_message(o);
    if (failure.failure_mode() ==
        Kokkos::Experimental::RawMemoryAllocationFailure::FailureMode::
            AllocationNotAligned) {
      // TODO: delete the misaligned memory?
      o << "Warning: Allocation failed due to misalignment; memory may "
           "be leaked.\n";
    }
    o.flush();
  };
  try {
    std::ostringstream sstr;
    generate_failure_message(sstr);
    Kokkos::Impl::throw_runtime_exception(sstr.str());
  } catch (std::bad_alloc const&) {
    // Probably failed to allocate the string because we're so close to out
    // of memory. Try printing to std::cerr instead
    try {
      generate_failure_message(std::cerr);
    } catch (std::bad_alloc const&) {
      // oh well, we tried...
    }
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos encountered an allocation failure, then another allocation "
        "failure while trying to create the error message.");
  }
}

}  // end namespace Impl
}  // end namespace Kokkos
