// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// This file calls most of the basic Kokkos primitives. When combined with a
// testing library this tests that our shared-library loading based profiling
// mechanisms work

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <cstdint>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
  Kokkos::initialize(argc, argv);
  {
    // This test only uses host kernel launch mechanisms. This is to allow for
    // the test to run on platforms where CUDA lambda launch isn't supported.
    // This is safe because this test only seeks to test that the dlsym-based
    // tool loading mechanisms work, all of which happens completely
    // independently of the enabled backends
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    using memory_space    = typename execution_space::memory_space;
    Kokkos::View<int*, memory_space> src_view("source", 10);
    Kokkos::View<int*, memory_space> dst_view("destination", 10);
    Kokkos::deep_copy(dst_view, src_view);
    Kokkos::parallel_for("parallel_for",
                         Kokkos::RangePolicy<execution_space>(0, 1),
                         [=](int i) { (void)i; });
    int result;
    Kokkos::parallel_reduce(
        "parallel_reduce", Kokkos::RangePolicy<execution_space>(0, 1),
        [=](int i, int& hold_result) { hold_result += i; }, result);
    Kokkos::parallel_scan("parallel_scan",
                          Kokkos::RangePolicy<execution_space>(0, 1),
                          [=](const int i, int& hold_result, const bool final) {
                            if (final) {
                              hold_result += i;
                            }
                          });
    Kokkos::Profiling::pushRegion("push_region");
    Kokkos::Profiling::popRegion();
    uint32_t sectionId;
    Kokkos::Profiling::createProfileSection("created_section", &sectionId);
    Kokkos::Profiling::startSection(sectionId);
    Kokkos::Profiling::stopSection(sectionId);
    Kokkos::Profiling::destroyProfileSection(sectionId);
    Kokkos::Profiling::markEvent("profiling_event");
    Kokkos::Tools::declareMetadata("dogs", "good");
  }
  Kokkos::finalize();
}
