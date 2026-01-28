// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSP_SCOPED_REGION_HPP
#define KOKKOSP_SCOPED_REGION_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_SCOPEDREGION
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Profiling.hpp>

#include <string>

namespace Kokkos::Profiling {

class [[nodiscard]] ScopedRegion {
 public:
  ScopedRegion(ScopedRegion const &)            = delete;
  ScopedRegion &operator=(ScopedRegion const &) = delete;

#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  [[nodiscard]]
#endif
  explicit ScopedRegion(std::string const &name) {
    Kokkos::Profiling::pushRegion(name);
  }
  ~ScopedRegion() { Kokkos::Profiling::popRegion(); }
};

}  // namespace Kokkos::Profiling

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_SCOPEDREGION
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_SCOPEDREGION
#endif
#endif
