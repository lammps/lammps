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

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_COMPILER_INTEL) && (KOKKOS_COMPILER_INTEL < 1800)

namespace {

// error: expression must have a constant value
//   std::enable_if_t<!has_deprecated_cuda_impl_initialize_v<T>>
constexpr bool
test_compiler_upgrade_needed_for_detection_idiom_and_variable_template() {
  return true;
}
static_assert(
    test_compiler_upgrade_needed_for_detection_idiom_and_variable_template(),
    "Intel C++ compiler is awesome");

}  // namespace

#else

// The purpose of this compile-only test is twofold:
// 1. mimic Legion's use of Kokkos implementation details for initializing the
//    exectution environment
// 2. demonstrate how to leverage SFINAE to support Kokkos version through the
//    ExecutionSpace::impl_initialize breaking change before release 3.7
namespace {
#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, "")  // FIXME C++17

#ifdef KOKKOS_ENABLE_CUDA
template <class T>
using deprecated_cuda_impl_initialize_t =
    decltype(T::impl_initialize(typename T::SelectDevice(0), 1));

template <class T>
constexpr bool has_deprecated_cuda_impl_initialize_v =
    Kokkos::is_detected<deprecated_cuda_impl_initialize_t, T>::value;

template <class T>
std::enable_if_t<has_deprecated_cuda_impl_initialize_v<T> >
legion_initialize_kokkos_cuda() {
  int cuda_device_id = 0;
  int num_instances  = 1;
  T::impl_initialize(typename T::SelectDevice(cuda_device_id), num_instances);
}

template <class T>
std::enable_if_t<!has_deprecated_cuda_impl_initialize_v<T> >
legion_initialize_kokkos_cuda() {
  int cuda_device_id = 0;
  auto const settings =
      Kokkos::InitializationSettings().set_device_id(cuda_device_id);
  T::impl_initialize(settings);
}

STATIC_ASSERT(std::is_void<
              decltype(legion_initialize_kokkos_cuda<Kokkos::Cuda>())>::value);
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template <class T>
using deprecated_openmp_impl_initialize_t = decltype(T::impl_initialize(0));

template <class T>
constexpr bool has_deprecated_openmp_impl_initialize_v =
    Kokkos::is_detected<deprecated_openmp_impl_initialize_t, T>::value;

template <class T>
std::enable_if_t<has_deprecated_openmp_impl_initialize_v<T> >
legion_initialize_kokkos_openmp() {
  int thread_count = -1;
  T::impl_initialize(thread_count);
}

template <class T>
std::enable_if_t<!has_deprecated_openmp_impl_initialize_v<T> >
legion_initialize_kokkos_openmp() {
  int thread_count = -1;
  auto const settings =
      Kokkos::InitializationSettings().set_num_threads(thread_count);
  T::impl_initialize(settings);
}

STATIC_ASSERT(std::is_void<decltype(
                  legion_initialize_kokkos_openmp<Kokkos::OpenMP>())>::value);

#endif

#ifdef KOKKOS_ENABLE_SERIAL
template <class T>
using deprecated_serial_impl_initialize_t = decltype(T::impl_initialize());

template <class T>
constexpr bool has_deprecated_serial_impl_initialize_v =
    Kokkos::is_detected<deprecated_serial_impl_initialize_t, T>::value;

template <class T>
std::enable_if_t<has_deprecated_serial_impl_initialize_v<T> >
legion_initialize_kokkos_serial() {
  T::impl_initialize();
}

template <class T>
std::enable_if_t<!has_deprecated_serial_impl_initialize_v<T> >
legion_initialize_kokkos_serial() {
  Kokkos::InitializationSettings settings;
  T::impl_initialize(settings);
}

STATIC_ASSERT(std::is_void<decltype(
                  legion_initialize_kokkos_serial<Kokkos::Serial>())>::value);
#endif

}  // namespace

#endif
