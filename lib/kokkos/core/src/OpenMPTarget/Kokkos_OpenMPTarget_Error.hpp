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

#ifndef KOKKOS_OPENMPTARGET_ERROR_HPP
#define KOKKOS_OPENMPTARGET_ERROR_HPP

#include <impl/Kokkos_Error.hpp>
#include <sstream>

namespace Kokkos {
namespace Impl {

inline void ompt_internal_safe_call(int e, const char* name,
                                    const char* file = nullptr,
                                    const int line   = 0) {
  if (e != 0) {
    std::ostringstream out;
    out << name << " return value of " << e << " indicates failure";
    if (file) {
      out << " " << file << ":" << line;
    }
    throw_runtime_exception(out.str());
  }
}

#define KOKKOS_IMPL_OMPT_SAFE_CALL(call) \
  Kokkos::Impl::ompt_internal_safe_call(call, #call, __FILE__, __LINE__)

}  // namespace Impl
}  // namespace Kokkos

#endif
