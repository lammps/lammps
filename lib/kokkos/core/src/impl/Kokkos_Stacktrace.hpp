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

#ifndef KOKKOS_STACKTRACE_HPP
#define KOKKOS_STACKTRACE_HPP

#include <functional>
#include <ostream>
#include <string>

namespace Kokkos {
namespace Impl {

/// \brief Return the demangled version of the input symbol, or the
///   original input if demangling is not possible.
std::string demangle(const std::string& name);

/// \brief Save the current stacktrace.
///
/// You may only save one stacktrace at a time.  If you call this
/// twice, the second call will overwrite the result of the first
/// call.
void save_stacktrace();

/// \brief Print the raw form of the currently saved stacktrace, if
///   any, to the given output stream.
void print_saved_stacktrace(std::ostream& out);

/// \brief Print the currently saved, demangled stacktrace, if any, to
///   the given output stream.
///
/// Demangling is best effort only.
void print_demangled_saved_stacktrace(std::ostream& out);

/// \brief Set the std::terminate handler so that it prints the
///   currently saved stack trace, then calls user_post.
///
/// This is useful if you want to call, say, MPI_Abort instead of
/// std::abort.  The MPI Standard frowns upon calling MPI functions
/// without including their header file, and Kokkos does not depend on
/// MPI, so there's no way for Kokkos to depend on MPI_Abort in a
/// portable way.
void set_kokkos_terminate_handler(std::function<void()> user_post = nullptr);

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_STACKTRACE_HPP
