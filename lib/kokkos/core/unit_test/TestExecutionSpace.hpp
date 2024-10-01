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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

template <class ExecutionSpace>
struct CheckClassWithExecutionSpaceAsDataMemberIsCopyable {
  Kokkos::DefaultExecutionSpace device;
  Kokkos::DefaultHostExecutionSpace host;

  KOKKOS_FUNCTION void operator()(int i, int& e) const { e += i; }

  CheckClassWithExecutionSpaceAsDataMemberIsCopyable() {
    int errors;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this,
                            errors);
    EXPECT_EQ(errors, 0);
  }
};

// FIXME_OPENMPTARGET nvlink error: Undefined reference to
// '_ZSt25__throw_bad_function_callv' in
// '/tmp/TestOpenMPTarget_ExecutionSpace-434d81.cubin'
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, execution_space_as_class_data_member) {
  CheckClassWithExecutionSpaceAsDataMemberIsCopyable<TEST_EXECSPACE>();
}
#endif

constexpr bool test_execspace_explicit_construction() {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifdef KOKKOS_ENABLE_SERIAL
  static_assert(std::is_convertible_v<Kokkos::NewInstance, Kokkos::Serial>);
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  static_assert(std::is_convertible_v<int, Kokkos::OpenMP>);
#endif
#ifdef KOKKOS_ENABLE_CUDA
  static_assert(std::is_convertible_v<cudaStream_t, Kokkos::Cuda>);
#endif
#ifdef KOKKOS_ENABLE_HIP
  static_assert(std::is_convertible_v<hipStream_t, Kokkos::HIP>);
#endif
#ifdef KOKKOS_ENABLE_HPX
  static_assert(std::is_convertible_v<Kokkos::Experimental::HPX::instance_mode,
                                      Kokkos::Experimental::HPX>);
  static_assert(
      std::is_convertible_v<hpx::execution::experimental::unique_any_sender<>&&,
                            Kokkos::Experimental::HPX>);
#endif
#else
#ifdef KOKKOS_ENABLE_SERIAL
  static_assert(!std::is_convertible_v<Kokkos::NewInstance, Kokkos::Serial>);
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  static_assert(!std::is_convertible_v<int, Kokkos::OpenMP>);
#endif
#ifdef KOKKOS_ENABLE_CUDA
  static_assert(!std::is_convertible_v<cudaStream_t, Kokkos::Cuda>);
#endif
#ifdef KOKKOS_ENABLE_HIP
  static_assert(!std::is_convertible_v<hipStream_t, Kokkos::HIP>);
#endif
#ifdef KOKKOS_ENABLE_HPX
  static_assert(!std::is_convertible_v<Kokkos::Experimental::HPX::instance_mode,
                                       Kokkos::Experimental::HPX>);
  static_assert(!std::is_convertible_v<
                hpx::execution::experimental::unique_any_sender<>&&,
                Kokkos::Experimental::HPX>);
#endif
#endif

#ifdef KOKKOS_ENABLE_OPENACC
  static_assert(!std::is_convertible_v<int, Kokkos::Experimental::OpenACC>);
#endif
#ifdef KOKKOS_ENABLE_SYCL
  static_assert(
      !std::is_convertible_v<sycl::queue, Kokkos::Experimental::SYCL>);
#endif

  return true;
}

static_assert(test_execspace_explicit_construction());

}  // namespace
