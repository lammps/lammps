// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
