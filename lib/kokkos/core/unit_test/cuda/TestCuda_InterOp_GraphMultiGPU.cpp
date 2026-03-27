// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "Kokkos_Core.hpp"
#include "Kokkos_Graph.hpp"
#include "TestCuda_Category.hpp"

namespace {

class TEST_CATEGORY_FIXTURE(multi_gpu) : public ::testing::Test {
 public:
  using graph_t          = Kokkos::Experimental::Graph<Kokkos::Cuda>;
  using graph_node_ref_t = Kokkos::Experimental::GraphNodeRef<Kokkos::Cuda>;

 public:
  void SetUp() override {
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&num_devices));

    if (num_devices < 2)
      GTEST_SKIP()
          << "Skipping multi-gpu testing since current machine only has "
          << num_devices << " GPU.";

    execs.reserve(num_devices);

    for (int idevice = 0; idevice < num_devices; ++idevice) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaSetDevice(idevice));
      cudaStream_t stream = nullptr;
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamCreate(&stream));
      execs.push_back(Kokkos::Cuda(stream));
    }
  }

  void TearDown() override {
    while (!execs.empty()) {
      const auto stream = execs.back().cuda_stream();
      execs.pop_back();
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamDestroy(stream));
    }
  }

 protected:
  int num_devices;
  std::vector<Kokkos::Cuda> execs;
};

template <typename ViewTypeA, typename ViewTypeB>
struct ThenFunctor {
  ViewTypeA data_A;
  ViewTypeB data_B;

  KOKKOS_FUNCTION void operator()() const {
    ++data_A();
    ++data_B();
  }
};

TEST_F(TEST_CATEGORY_FIXTURE(multi_gpu), then) {
  using view_t        = Kokkos::View<unsigned int, Kokkos::CudaSpace>;
  using shared_view_t = Kokkos::View<unsigned int, Kokkos::SharedSpace>;

  const shared_view_t shared(Kokkos::view_alloc("shared data"));

  std::vector<view_t> views(num_devices);

  // Create the view on each device, but not yet the nodes.
  for (int idevice = 0; idevice < num_devices; ++idevice) {
    views.at(idevice) = view_t(Kokkos::view_alloc(
        "data on device " + std::to_string(idevice), execs.at(idevice)));
  }

  graph_t graph{};

  graph_node_ref_t pred = graph.root_node();

  // Now, create the nodes. Doing it in a loop separated from the view
  // construction ensures that the graph itself correctly calls cudaSetDevice
  // when adding a node.
  for (int idevice = 0; idevice < num_devices; ++idevice) {
    pred = pred.then(
        "then node device " + std::to_string(idevice), execs.at(idevice),
        ThenFunctor<view_t, shared_view_t>{views.at(idevice), shared});
  }

  // Fence execs on all devices but the device 0 to make sure views are properly
  // allocated and initialized.
  for (int idevice = 1; idevice < num_devices; ++idevice) {
    execs.at(idevice).fence();
  }

  graph.submit(execs.at(0));

  execs.at(0).fence();

  ASSERT_EQ(shared(), num_devices);

  for (int idevice = 0; idevice < num_devices; ++idevice) {
    ASSERT_EQ(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                                  views.at(idevice))(),
              1);
  }
}

}  // namespace
