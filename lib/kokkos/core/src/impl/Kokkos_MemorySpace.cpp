/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
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

/** @file Kokkos_MemorySpace.cpp
 *
 *  Operations common to memory space instances, or at least default
 *  implementations thereof.
 */

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
