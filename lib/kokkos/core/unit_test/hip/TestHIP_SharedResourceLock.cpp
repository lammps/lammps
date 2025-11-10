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

#include <random>

#include <Kokkos_Core.hpp>
#include <TestHIP_Category.hpp>

template <typename ViewType>
struct AddTo {
  static_assert(ViewType::rank() == 0);

  ViewType data;
  typename ViewType::value_type value;

  template <typename T>
  KOKKOS_FUNCTION void operator()(const T) const {
    data() += value;
  }

  std::byte unused[Kokkos::Impl::HIPTraits::ConstantMemoryUseThreshold] = {};
};

template <typename ViewType>
struct ThreadWorkOnConstantMemory {
  using functor_t = AddTo<ViewType>;

  static_assert(sizeof(functor_t) >
                Kokkos::Impl::HIPTraits::ConstantMemoryUseThreshold);
  static_assert(sizeof(functor_t) <
                Kokkos::Impl::HIPTraits::ConstantMemoryUsage);

  std::vector<std::chrono::milliseconds> sleep_for;
  std::vector<Kokkos::HIP> execs;
  ViewType data;
  typename ViewType::value_type value;

  void operator()() && {
    ASSERT_EQ(execs.size(), sleep_for.size());
    size_t irep = 0;
    for (auto&& exec : execs) {
      std::this_thread::sleep_for(sleep_for.at(irep++));
      Kokkos::parallel_for(Kokkos::RangePolicy(exec, 0, 1),
                           functor_t{data, value});
    }
  }
};

struct TEST_CATEGORY_FIXTURE(constant_memory) : public ::testing::Test {
  static constexpr size_t nthreads_per_device = 2 << 5;
  static constexpr size_t nreps_per_thread    = 2 << 3;

  using view_t = Kokkos::View<size_t[nthreads_per_device], Kokkos::SharedSpace>;
  using view_r0_t = Kokkos::View<size_t, Kokkos::SharedSpace>;

  using work_t = ThreadWorkOnConstantMemory<view_r0_t>;

  static int get_device_count() {
    int num_devices = 0;
    KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&num_devices));
    return num_devices;
  }

  int num_devices{get_device_count()};
  std::random_device rnd_dev{};
  std::mt19937 twister{rnd_dev()};
  std::uniform_int_distribution<size_t> dist_value{1, nthreads_per_device};

  // Sleep as long as the launch latency [µs].
  std::uniform_int_distribution<size_t> dist_sleep_for{1, 5};

  std::vector<std::array<std::thread, nthreads_per_device>> threads_per_device{
      static_cast<size_t>(num_devices)};

  view_t accu_per_device{Kokkos::view_alloc("per-device counter")};
  std::vector<size_t> expd_per_device = std::vector<size_t>(num_devices);

  std::vector<hipStream_t> streams_parent{static_cast<size_t>(num_devices),
                                          nullptr};

  auto get_value(const int device) {
    const auto value = dist_value(twister);
    expd_per_device.at(device) += nreps_per_thread * value;
    return value;
  }

  auto get_sleep_for() {
    std::vector<std::chrono::milliseconds> sleep_for(nreps_per_thread);
    std::generate(sleep_for.begin(), sleep_for.end(), [&]() {
      return std::chrono::milliseconds{dist_sleep_for(twister)};
    });
    return sleep_for;
  }

  // Wait for all threads to finish.
  void finalize() {
    for (int device = 0; device < num_devices; ++device) {
      for (auto& t : threads_per_device.at(device)) {
        t.join();
      }

      KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(device));
      KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize());

      EXPECT_EQ(accu_per_device(device), expd_per_device.at(device)) << device;

      KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(streams_parent.at(device)));
    }
  }
};

TEST_F(TEST_CATEGORY_FIXTURE(constant_memory), many_streams_per_device) {
  ASSERT_GT(this->get_device_count(), 0);

  for (int device = 0; device < this->num_devices; ++device) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(device));
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&this->streams_parent[device]));

    const auto execs = Kokkos::Experimental::partition_space(
        Kokkos::HIP(this->streams_parent.at(device),
                    Kokkos::Impl::ManageStream::no),
        std::vector<size_t>(nreps_per_thread, 1));

    for (auto& t : this->threads_per_device[device]) {
      t = std::thread(work_t{this->get_sleep_for(), execs,
                             Kokkos::subview(this->accu_per_device, device),
                             this->get_value(device)});
    }
  }

  this->finalize();
}

TEST_F(TEST_CATEGORY_FIXTURE(constant_memory), one_stream_per_device) {
  ASSERT_GT(this->get_device_count(), 0);

  for (int device = 0; device < this->num_devices; ++device) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(device));
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&this->streams_parent[device]));

    const std::vector<Kokkos::HIP> execs{
        nreps_per_thread, Kokkos::HIP(this->streams_parent.at(device),
                                      Kokkos::Impl::ManageStream::no)};

    for (auto& t : this->threads_per_device[device]) {
      t = std::thread(work_t{this->get_sleep_for(), execs,
                             Kokkos::subview(this->accu_per_device, device),
                             this->get_value(device)});
    }
  }

  this->finalize();
}
