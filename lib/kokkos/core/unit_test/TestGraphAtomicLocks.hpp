// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

namespace {

struct InitTag {};
struct WorkTag {};

template <typename ViewType>
struct TestFunctor {
  ViewType data;

  static_assert(ViewType::rank() == 0);
  static_assert(
      std::is_same_v<typename ViewType::value_type, Kokkos::complex<double>>);

  template <typename T>
  KOKKOS_FUNCTION void operator()(const InitTag, const T) const {
    data() = Kokkos::complex<double>{0., 0.};
  }

  template <typename T>
  KOKKOS_FUNCTION void operator()(const WorkTag, const T) const {
    Kokkos::atomic_add(&data(), Kokkos::complex<double>{1., 1.});
  }
};

// This test serves to ensure that lock-based atomic operations work in
// a graph on device. In particular, this test serves to ensure that the
// lock arrays needed for such operations be on device.
TEST(TEST_CATEGORY, graph_lock_based_atomic_op) {
  TEST_EXECSPACE ex{};

  // Don't initialize here to avoid that the initialization triggers a kernel
  // launch that ensures that the lock arrays are on device. We want to make
  // sure they are on device even without a preceding kernel launch.
  Kokkos::View<Kokkos::complex<double>, TEST_EXECSPACE> result(
      Kokkos::view_alloc(Kokkos::WithoutInitializing));

  auto graph = Kokkos::Experimental::create_graph(ex, [&](const auto& root) {
    root.then_parallel_for(
            Kokkos::RangePolicy<TEST_EXECSPACE, InitTag>(0, 1),
            TestFunctor<Kokkos::View<Kokkos::complex<double>, TEST_EXECSPACE>>{
                result})
        .then_parallel_for(
            Kokkos::RangePolicy<TEST_EXECSPACE, WorkTag>(0, 100),
            TestFunctor<Kokkos::View<Kokkos::complex<double>, TEST_EXECSPACE>>{
                result});
  });

  graph.submit(ex);

  // Check.
  Kokkos::complex<double> result_h;
  Kokkos::deep_copy(ex, result_h, result);
  ex.fence();
  ASSERT_FLOAT_EQ(result_h.real(), 100.);
  ASSERT_FLOAT_EQ(result_h.imag(), 100.);
}

}  // end namespace
