// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TIMER_HPP
#define KOKKOS_TIMER_HPP

#include <Kokkos_Macros.hpp>

#include <chrono>

namespace Kokkos {

class Timer {
  using Clock = std::chrono::high_resolution_clock;
  Clock::time_point start_;

 public:
  Timer(const Timer&)            = delete;
  Timer& operator=(const Timer&) = delete;

  Timer() { reset(); }

  void reset() { start_ = Clock::now(); }

  double seconds() const {
    using namespace std::chrono;
    auto const now = Clock::now();
    return duration_cast<duration<double>>(now - start_).count();
  }
};

}  // namespace Kokkos

#endif
