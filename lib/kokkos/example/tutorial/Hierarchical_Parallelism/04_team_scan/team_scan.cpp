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

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Timer.hpp>
#include <cstdio>
#include <cstdlib>

using Device = Kokkos::DefaultExecutionSpace;
using Host   = Kokkos::HostSpace::execution_space;

using team_policy = Kokkos::TeamPolicy<Device>;
using team_member = team_policy::member_type;

static const int TEAM_SIZE = 16;

struct find_2_tuples {
  int chunk_size;
  Kokkos::View<const int*> data;
  Kokkos::View<int**> histogram;

  find_2_tuples(int chunk_size_, Kokkos::DualView<int*> data_,
                Kokkos::DualView<int**> histogram_)
      : chunk_size(chunk_size_),
        data(data_.d_view),
        histogram(histogram_.d_view) {
    data_.sync<Device>();
    histogram_.sync<Device>();
    histogram_.modify<Device>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& dev) const {
    Kokkos::View<int**, Kokkos::MemoryUnmanaged> l_histogram(
        dev.team_shmem(), TEAM_SIZE, TEAM_SIZE);
    Kokkos::View<int*, Kokkos::MemoryUnmanaged> l_data(dev.team_shmem(),
                                                       chunk_size + 1);

    const int i = dev.league_rank() * chunk_size;
    for (int j = dev.team_rank(); j < chunk_size + 1; j += dev.team_size())
      l_data(j) = data(i + j);

    for (int k = dev.team_rank(); k < TEAM_SIZE; k += dev.team_size())
      for (int l = 0; l < TEAM_SIZE; l++) l_histogram(k, l) = 0;
    dev.team_barrier();

    for (int j = 0; j < chunk_size; j++) {
      for (int k = dev.team_rank(); k < TEAM_SIZE; k += dev.team_size())
        for (int l = 0; l < TEAM_SIZE; l++) {
          if ((l_data(j) == k) && (l_data(j + 1) == l)) l_histogram(k, l)++;
        }
    }

    for (int k = dev.team_rank(); k < TEAM_SIZE; k += dev.team_size())
      for (int l = 0; l < TEAM_SIZE; l++) {
        Kokkos::atomic_fetch_add(&histogram(k, l), l_histogram(k, l));
      }
    dev.team_barrier();
  }
  size_t team_shmem_size(int team_size) const {
    return Kokkos::View<int**, Kokkos::MemoryUnmanaged>::shmem_size(TEAM_SIZE,
                                                                    TEAM_SIZE) +
           Kokkos::View<int*, Kokkos::MemoryUnmanaged>::shmem_size(chunk_size +
                                                                   1);
  }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg, args);

  {
    int chunk_size = 1024;
    int nchunks    = 100000;  // 1024*1024;
    Kokkos::DualView<int*> data("data", nchunks * chunk_size + 1);

    srand(1231093);

    for (int i = 0; i < (int)data.extent(0); i++) {
      data.h_view(i) = rand() % TEAM_SIZE;
    }
    data.modify<Host>();
    data.sync<Device>();

    Kokkos::DualView<int**> histogram("histogram", TEAM_SIZE, TEAM_SIZE);

    Kokkos::Timer timer;
    // threads/team is automatically limited to maximum supported by the device.
    int const concurrency = Device::execution_space().concurrency();
    int team_size         = TEAM_SIZE;
    if (team_size > concurrency) team_size = concurrency;
    Kokkos::parallel_for(team_policy(nchunks, team_size),
                         find_2_tuples(chunk_size, data, histogram));
    Kokkos::fence();
    double time = timer.seconds();

    histogram.sync<Host>();

    printf("Time: %f \n\n", time);
    int sum = 0;
    for (int k = 0; k < TEAM_SIZE; k++) {
      for (int l = 0; l < TEAM_SIZE; l++) {
        printf("%i ", histogram.h_view(k, l));
        sum += histogram.h_view(k, l);
      }
      printf("\n");
    }
    printf("Result: %i %i\n", sum, chunk_size * nchunks);
  }
  Kokkos::finalize();
}
