// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSP_PROFILE_SECTION_HPP
#define KOKKOSP_PROFILE_SECTION_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_PROFILESECTION
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Profiling.hpp>

#include <string>

namespace Kokkos::Profiling {

class [[nodiscard]] ProfilingSection {
  uint32_t sectionID;

 public:
  ProfilingSection(ProfilingSection const&)            = delete;
  ProfilingSection& operator=(ProfilingSection const&) = delete;

#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  [[nodiscard]]
#endif
  explicit ProfilingSection(const std::string& sectionName) {
    Kokkos::Profiling::createProfileSection(sectionName, &sectionID);
  }

  void start() { Kokkos::Profiling::startSection(sectionID); }

  void stop() { Kokkos::Profiling::stopSection(sectionID); }

  ~ProfilingSection() { Kokkos::Profiling::destroyProfileSection(sectionID); }
};

}  // namespace Kokkos::Profiling

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_PROFILESECTION
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_PROFILESECTION
#endif
#endif
