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

#ifndef KOKKOS_GIT_VERSION_INFO_H
#define KOKKOS_GIT_VERSION_INFO_H

#include <string>

namespace Kokkos {
namespace Impl {

extern std::string GIT_BRANCH;
extern std::string GIT_COMMIT_HASH;
extern std::string GIT_CLEAN_STATUS;
extern std::string GIT_COMMIT_DESCRIPTION;
extern std::string GIT_COMMIT_DATE;

}  // namespace Impl
}  // namespace Kokkos

#endif
