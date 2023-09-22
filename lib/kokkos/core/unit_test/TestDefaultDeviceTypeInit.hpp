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

#ifdef KOKKOS_ENABLE_OPENMP
#include <omp.h>
#endif
#include <set>
#if !defined(KOKKOS_ENABLE_CUDA) || defined(__CUDACC__)

namespace Test {

namespace Impl {

std::set<char*> delete_these;
void cleanup_memory() {
  for (auto x : delete_these) {
    delete[] x;
  }
}

char** init_kokkos_args(bool do_threads, bool do_numa, bool do_device,
                        bool do_other, bool do_tune, int& nargs,
                        Kokkos::InitArguments& init_args) {
  nargs = (do_threads ? 1 : 0) + (do_numa ? 1 : 0) + (do_device ? 1 : 0) +
          (do_other ? 4 : 0) + (do_tune ? 1 : 0);

  char** args_kokkos      = new char*[nargs];
  const int max_args_size = 45;
  for (int i = 0; i < nargs; i++) {
    args_kokkos[i] = new char[max_args_size];
    delete_these.insert(args_kokkos[i]);
  }

  int threads_idx = do_other ? 1 : 0;
  int numa_idx    = (do_other ? 3 : 0) + (do_threads ? 1 : 0);
  int device_idx =
      (do_other ? 3 : 0) + (do_threads ? 1 : 0) + (do_numa ? 1 : 0);
  int tune_idx = (do_other ? 4 : 0) + (do_threads ? 1 : 0) + (do_numa ? 1 : 0) +
                 (do_device ? 1 : 0);

  if (do_threads) {
    int nthreads = 3;

#ifdef KOKKOS_ENABLE_OPENMP
    if (omp_get_max_threads() < nthreads) {
      nthreads = omp_get_max_threads();
    }
#elif defined(KOKKOS_ENABLE_HPX)
    const int concurrency = std::thread::hardware_concurrency();
    if (concurrency < nthreads) {
      nthreads = concurrency;
    }
#endif

    if (Kokkos::hwloc::available()) {
      if (Kokkos::hwloc::get_available_threads_per_core() <
          static_cast<unsigned>(nthreads))
        nthreads = Kokkos::hwloc::get_available_threads_per_core() *
                   Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::Serial, Kokkos::DefaultExecutionSpace>::value ||
        std::is_same<Kokkos::Serial,
                     Kokkos::DefaultHostExecutionSpace>::value) {
      nthreads = 1;
    }
#endif

    init_args.num_threads = nthreads;
    snprintf(args_kokkos[threads_idx], max_args_size, "--threads=%i", nthreads);
  }

  if (do_numa) {
    int numa = 1;
    if (Kokkos::hwloc::available()) {
      numa = Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::Serial, Kokkos::DefaultExecutionSpace>::value ||
        std::is_same<Kokkos::Serial,
                     Kokkos::DefaultHostExecutionSpace>::value) {
      numa = 1;
    }
#endif

    init_args.num_numa = numa;
    snprintf(args_kokkos[numa_idx], max_args_size, "--numa=%i", numa);
  }

  if (do_device) {
    init_args.device_id = 0;
    snprintf(args_kokkos[device_idx], max_args_size, "--device-id=%i", 0);
  }

  if (do_other) {
    snprintf(args_kokkos[0], max_args_size, "--dummyarg=1");
    snprintf(args_kokkos[threads_idx + (do_threads ? 1 : 0)], max_args_size,
             "--dummy2arg");
    snprintf(args_kokkos[threads_idx + (do_threads ? 1 : 0) + 1], max_args_size,
             "dummy3arg");
    snprintf(args_kokkos[device_idx + (do_device ? 1 : 0)], max_args_size,
             "dummy4arg=1");
  }

  if (do_tune) {
    init_args.tune_internals = true;
    snprintf(args_kokkos[tune_idx], max_args_size, "--kokkos-tune-internals");
  }

  return args_kokkos;
}

Kokkos::InitArguments init_initstruct(bool do_threads, bool do_numa,
                                      bool do_device, bool do_tune) {
  Kokkos::InitArguments args;

  if (do_threads) {
    int nthreads = 3;

#ifdef KOKKOS_ENABLE_OPENMP
    if (omp_get_max_threads() < nthreads) {
      nthreads = omp_get_max_threads();
    }
#elif defined(KOKKOS_ENABLE_HPX)
    const int concurrency = std::thread::hardware_concurrency();
    if (concurrency < nthreads) {
      nthreads = concurrency;
    }
#endif

    if (Kokkos::hwloc::available()) {
      if (Kokkos::hwloc::get_available_threads_per_core() <
          static_cast<unsigned>(nthreads)) {
        nthreads = Kokkos::hwloc::get_available_threads_per_core() *
                   Kokkos::hwloc::get_available_numa_count();
      }
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::Serial, Kokkos::DefaultExecutionSpace>::value ||
        std::is_same<Kokkos::Serial,
                     Kokkos::DefaultHostExecutionSpace>::value) {
      nthreads = 1;
    }
#endif

    args.num_threads = nthreads;
  }

  if (do_numa) {
    int numa = 1;
    if (Kokkos::hwloc::available()) {
      numa = Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::Serial, Kokkos::DefaultExecutionSpace>::value ||
        std::is_same<Kokkos::Serial,
                     Kokkos::DefaultHostExecutionSpace>::value) {
      numa = 1;
    }
#endif

    args.num_numa = numa;
  }

  if (do_device) {
    args.device_id = 0;
  }

  if (do_tune) {
    args.tune_internals = true;
  }

  return args;
}

void check_correct_initialization(const Kokkos::InitArguments& argstruct) {
  ASSERT_EQ(Kokkos::DefaultExecutionSpace::impl_is_initialized(), 1);
  ASSERT_EQ(Kokkos::HostSpace::execution_space::impl_is_initialized(), 1);

  // Figure out the number of threads the HostSpace ExecutionSpace should have
  // initialized to.
  int expected_nthreads = argstruct.num_threads;

#ifdef KOKKOS_ENABLE_OPENMP
  if (std::is_same<Kokkos::HostSpace::execution_space, Kokkos::OpenMP>::value) {
    // use openmp default num threads
    if (expected_nthreads < 0 ||
        (expected_nthreads == 0 && !Kokkos::hwloc::available())) {
      expected_nthreads = omp_get_max_threads();
    }
    // use hwloc if available
    else if (expected_nthreads == 0 && Kokkos::hwloc::available()) {
      expected_nthreads = Kokkos::hwloc::get_available_numa_count() *
                          Kokkos::hwloc::get_available_cores_per_numa() *
                          Kokkos::hwloc::get_available_threads_per_core();
    }
  }
#endif

  if (expected_nthreads < 1) {
    if (Kokkos::hwloc::available()) {
      expected_nthreads = Kokkos::hwloc::get_available_numa_count() *
                          Kokkos::hwloc::get_available_cores_per_numa() *
                          Kokkos::hwloc::get_available_threads_per_core();
    } else {
      expected_nthreads = 1;
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Serial>::value ||
        std::is_same<Kokkos::DefaultHostExecutionSpace,
                     Kokkos::Serial>::value) {
      expected_nthreads = 1;
    }
#endif

#ifdef KOKKOS_ENABLE_HPX
    // HPX uses all cores on machine by default. Skip this test.
    if (std::is_same<Kokkos::DefaultExecutionSpace,
                     Kokkos::Experimental::HPX>::value ||
        std::is_same<Kokkos::DefaultHostExecutionSpace,
                     Kokkos::Experimental::HPX>::value) {
      return;
    }
#endif
  }

  int expected_numa = argstruct.num_numa;

  if (expected_numa < 1) {
    if (Kokkos::hwloc::available()) {
      expected_numa = Kokkos::hwloc::get_available_numa_count();
    } else {
      expected_numa = 1;
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Serial>::value ||
        std::is_same<Kokkos::DefaultHostExecutionSpace, Kokkos::Serial>::value)
      expected_numa = 1;
#endif
  }

  ASSERT_EQ(Kokkos::HostSpace::execution_space().impl_thread_pool_size(),
            expected_nthreads);

#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Cuda>::value) {
    int device;
    cudaGetDevice(&device);

    int expected_device = argstruct.device_id;
    if (argstruct.device_id < 0) {
      expected_device = Kokkos::Cuda().cuda_device();
    }

    ASSERT_EQ(expected_device, device);
  }
#endif
  ASSERT_EQ(argstruct.tune_internals, Kokkos::tune_internals());
}

// TODO: Add check whether correct number of threads are actually started.
void test_no_arguments() {
  Kokkos::initialize();
  check_correct_initialization(Kokkos::InitArguments());
  Kokkos::finalize();
}

void test_commandline_args(int nargs, char** args,
                           const Kokkos::InitArguments& argstruct) {
  Kokkos::initialize(nargs, args);
  check_correct_initialization(argstruct);
  Kokkos::finalize();
}

void test_initstruct_args(const Kokkos::InitArguments& args) {
  Kokkos::initialize(args);
  check_correct_initialization(args);
  Kokkos::finalize();
}

}  // namespace Impl

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_01
TEST(defaultdevicetypeinit, no_args) { Impl::test_no_arguments(); }
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_02
TEST(defaultdevicetypeinit, commandline_args_empty) {
  Kokkos::InitArguments argstruct;
  int nargs   = 0;
  char** args = Impl::init_kokkos_args(false, false, false, false, false, nargs,
                                       argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_03
TEST(defaultdevicetypeinit, commandline_args_other) {
  Kokkos::InitArguments argstruct;
  int nargs   = 0;
  char** args = Impl::init_kokkos_args(false, false, false, true, false, nargs,
                                       argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_04
TEST(defaultdevicetypeinit, commandline_args_nthreads) {
  Kokkos::InitArguments argstruct;
  int nargs   = 0;
  char** args = Impl::init_kokkos_args(true, false, false, false, false, nargs,
                                       argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_05
TEST(defaultdevicetypeinit, commandline_args_nthreads_numa) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(true, true, false, false, false, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();

  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_06
TEST(defaultdevicetypeinit, commandline_args_nthreads_numa_device) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(true, true, true, false, false, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();

  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_07
TEST(defaultdevicetypeinit, commandline_args_nthreads_device) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(true, false, true, false, false, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_08
TEST(defaultdevicetypeinit, commandline_args_numa_device) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(false, true, true, false, false, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_09
TEST(defaultdevicetypeinit, commandline_args_device) {
  Kokkos::InitArguments argstruct;
  int nargs   = 0;
  char** args = Impl::init_kokkos_args(false, false, true, false, false, nargs,
                                       argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);

  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_10
TEST(defaultdevicetypeinit, commandline_args_nthreads_numa_device_other) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(true, true, true, true, false, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);
  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_11
TEST(defaultdevicetypeinit, commandline_args_nthreads_numa_device_other_tune) {
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args =
      Impl::init_kokkos_args(true, true, true, true, true, nargs, argstruct);
  Impl::test_commandline_args(nargs, args, argstruct);
  Impl::cleanup_memory();
  delete[] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_12
TEST(defaultdevicetypeinit, initstruct_default) {
  Kokkos::InitArguments args;
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_13
TEST(defaultdevicetypeinit, initstruct_nthreads) {
  Kokkos::InitArguments args = Impl::init_initstruct(true, false, false, false);
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_14
TEST(defaultdevicetypeinit, initstruct_nthreads_numa) {
  Kokkos::InitArguments args = Impl::init_initstruct(true, true, false, false);
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_15
TEST(defaultdevicetypeinit, initstruct_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(false, false, true, false);
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_16
TEST(defaultdevicetypeinit, initstruct_nthreads_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(true, false, true, false);
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_17
TEST(defaultdevicetypeinit, initstruct_nthreads_numa_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(true, true, true, false);
  Impl::test_initstruct_args(args);
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_18
TEST(defaultdevicetypeinit, initstruct_nthreads_numa_device_tune) {
  Kokkos::InitArguments args = Impl::init_initstruct(true, true, true, true);
  Impl::test_initstruct_args(args);
}
#endif

}  // namespace Test

#endif
