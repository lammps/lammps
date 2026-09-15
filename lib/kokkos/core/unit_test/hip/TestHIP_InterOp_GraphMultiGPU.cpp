// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "Kokkos_Core.hpp"
#include "Kokkos_Graph.hpp"
#include "impl/Kokkos_DeviceManagement.hpp"
#include "TestHIP_Category.hpp"

#include <algorithm>
#include <vector>

namespace {

class TEST_CATEGORY_FIXTURE(multi_gpu) : public ::testing::Test {
 public:
  using graph_t = Kokkos::Experimental::Graph<Kokkos::HIP>;

 public:
  void SetUp() override {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&num_devices));

    if (num_devices < 2)
      GTEST_SKIP()
          << "Skipping multi-gpu testing since current machine only has "
          << num_devices << " GPU.";

    execs.reserve(num_devices);

    for (int idevice = 0; idevice < num_devices; ++idevice) {
      KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(idevice));
      hipStream_t stream = nullptr;
      KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&stream));
      execs.push_back(Kokkos::HIP(stream));
    }
  }

  void TearDown() override {
    while (!execs.empty()) {
      const auto stream = execs.back().hip_stream();
      execs.pop_back();
      KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(stream));
    }
  }

 protected:
  int num_devices;
  std::vector<Kokkos::HIP> execs;
};

// FIXME Copy/pasted from TestGraph.hpp.
template <typename ViewType, size_t Count>
struct SizedFunctor {
 public:
  static constexpr size_t count = Count;

  ViewType data;

  SizedFunctor(ViewType data_) : data(std::move(data_)) {}

  KOKKOS_FUNCTION void operator()() const { ++data(); }

 private:
  std::byte unused[count] = {};
};

// Create a graph with a then node on the first device that is not the default
// device. It should exercise the global launch mechanism.
TEST_F(TEST_CATEGORY_FIXTURE(multi_gpu),
       then_force_global_launch_on_non_default_device) {
  using view_t = Kokkos::View<int, Kokkos::HIPSpace>;
  using functor_t =
      SizedFunctor<view_t, Kokkos::Impl::HIPTraits::ConstantMemoryUsage + 1>;

  const auto exec_it = std::ranges::find_if(
      this->execs, [default_device_id =
                        Kokkos::HIP().hip_device()](const Kokkos::HIP& exec) {
        return exec.hip_device() != default_device_id;
      });
  ASSERT_NE(exec_it, std::cend(this->execs));

  const view_t data(Kokkos::view_alloc(
      "data on device " + std::to_string(exec_it->hip_device()), *exec_it));

  const auto device_handle = Kokkos::Experimental::get_device_handle(*exec_it);

  const auto graph = Kokkos::Experimental::create_graph(
      device_handle, [&](const auto& root) { root.then(functor_t(data)); });

  graph.submit(*exec_it);

  int value = 0;

  Kokkos::deep_copy(*exec_it, value, data);
  exec_it->fence();

  ASSERT_EQ(value, 1);
}

}  // namespace
