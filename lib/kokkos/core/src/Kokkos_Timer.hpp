/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
