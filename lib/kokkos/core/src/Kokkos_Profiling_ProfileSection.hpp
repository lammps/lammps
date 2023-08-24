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

#ifndef KOKKOSP_PROFILE_SECTION_HPP
#define KOKKOSP_PROFILE_SECTION_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_PROFILING_PROFILESECTION
#endif

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_Profiling.hpp>

#include <string>

namespace Kokkos {
namespace Profiling {

class ProfilingSection {
 public:
  ProfilingSection(ProfilingSection const&) = delete;
  ProfilingSection& operator=(ProfilingSection const&) = delete;

  ProfilingSection(const std::string& sectionName) {
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::createProfileSection(sectionName, &secID);
    }
  }

  void start() {
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::startSection(secID);
    }
  }

  void stop() {
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::stopSection(secID);
    }
  }

  ~ProfilingSection() {
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::destroyProfileSection(secID);
    }
  }

 protected:
  uint32_t secID;
};

}  // namespace Profiling
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PROFILING_PROFILESECTION
#endif
#endif
