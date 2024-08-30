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

#ifndef KOKKOS_TIMER_HPP
#define KOKKOS_TIMER_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_TIMER
#endif

#include <Kokkos_Macros.hpp>
// gcc 10.3.0 with CUDA doesn't support std::chrono,
// see https://github.com/kokkos/kokkos/issues/4334
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU == 1030) && \
    defined(KOKKOS_COMPILER_NVCC)
#include <sys/time.h>
#else
#include <chrono>
#endif

namespace Kokkos {

/** \brief  Time since construction */

#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU == 1030) && \
    defined(KOKKOS_COMPILER_NVCC)
class Timer {
 private:
  struct timeval m_old;

 public:
  inline void reset() { gettimeofday(&m_old, nullptr); }

  inline ~Timer() = default;

  inline Timer() { reset(); }

  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;

  inline double seconds() const {
    struct timeval m_new;

    gettimeofday(&m_new, nullptr);

    return ((double)(m_new.tv_sec - m_old.tv_sec)) +
           ((double)(m_new.tv_usec - m_old.tv_usec) * 1.0e-6);
  }
};
#else
class Timer {
 private:
  std::chrono::high_resolution_clock::time_point m_old;

 public:
  inline void reset() { m_old = std::chrono::high_resolution_clock::now(); }

  inline ~Timer() = default;

  inline Timer() { reset(); }

  Timer(const Timer&);
  Timer& operator=(const Timer&);

  inline double seconds() const {
    std::chrono::high_resolution_clock::time_point m_new =
        std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double> >(m_new -
                                                                      m_old)
        .count();
  }
};
#endif

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_TIMER
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_TIMER
#endif
#endif /* #ifndef KOKKOS_TIMER_HPP */
