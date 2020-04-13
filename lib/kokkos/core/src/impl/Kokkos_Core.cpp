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
#include <impl/Kokkos_Error.hpp>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stack>
#include <cerrno>
#include <unistd.h>

//----------------------------------------------------------------------------

namespace {
bool g_is_initialized = false;
bool g_show_warnings  = true;
std::stack<std::function<void()> > finalize_hooks;
}  // namespace

namespace Kokkos {
namespace Impl {
namespace {

bool is_unsigned_int(const char* str) {
  const size_t len = strlen(str);
  for (size_t i = 0; i < len; ++i) {
    if (!isdigit(str[i])) {
      return false;
    }
  }
  return true;
}
void initialize_internal(const InitArguments& args) {
// This is an experimental setting
// For KNL in Flat mode this variable should be set, so that
// memkind allocates high bandwidth memory correctly.
#ifdef KOKKOS_ENABLE_HBWSPACE
  setenv("MEMKIND_HBW_NODES", "1", 0);
#endif

  if (args.disable_warnings) {
    g_show_warnings = false;
  }

  // Protect declarations, to prevent "unused variable" warnings.
#if defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS) || \
    defined(KOKKOS_ENABLE_OPENMPTARGET) || defined(KOKKOS_ENABLE_HPX)
  const int num_threads = args.num_threads;
#endif
#if defined(KOKKOS_ENABLE_THREADS) || defined(KOKKOS_ENABLE_OPENMPTARGET)
  const int use_numa = args.num_numa;
#endif
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_ROCM)
  int use_gpu           = args.device_id;
  const int ndevices    = args.ndevices;
  const int skip_device = args.skip_device;
  // if the exact device is not set, but ndevices was given, assign round-robin
  // using on-node MPI rank
  if (use_gpu < 0 && ndevices >= 0) {
    auto local_rank_str = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");  // OpenMPI
    if (!local_rank_str)
      local_rank_str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK");  // MVAPICH2
    if (!local_rank_str)
      local_rank_str = std::getenv("SLURM_LOCALID");  // SLURM
    if (local_rank_str) {
      auto local_rank = std::atoi(local_rank_str);
      use_gpu         = local_rank % ndevices;
    } else {
      // user only gave us ndevices, but the MPI environment variable wasn't
      // set. start with GPU 0 at this point
      use_gpu = 0;
    }
    // shift assignments over by one so no one is assigned to "skip_device"
    if (use_gpu >= skip_device) ++use_gpu;
  }
#endif  // defined( KOKKOS_ENABLE_CUDA )

#if defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same<Kokkos::OpenMP, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::OpenMP, Kokkos::HostSpace::execution_space>::value) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    Kokkos::OpenMP::initialize(num_threads);
#else
    Kokkos::OpenMP::impl_initialize(num_threads);
#endif
  } else {
    // std::cout << "Kokkos::initialize() fyi: OpenMP enabled but not
    // initialized" << std::endl ;
  }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (std::is_same<Kokkos::Threads, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Threads,
                   Kokkos::HostSpace::execution_space>::value) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    if (num_threads > 0) {
      if (use_numa > 0) {
        Kokkos::Threads::initialize(num_threads, use_numa);
      } else {
        Kokkos::Threads::initialize(num_threads);
      }
    } else {
      Kokkos::Threads::initialize();
    }
#else
    if (num_threads > 0) {
      if (use_numa > 0) {
        Kokkos::Threads::impl_initialize(num_threads, use_numa);
      } else {
        Kokkos::Threads::impl_initialize(num_threads);
      }
    } else {
      Kokkos::Threads::impl_initialize();
    }
#endif
    // std::cout << "Kokkos::initialize() fyi: Pthread enabled and initialized"
    // << std::endl ;
  } else {
    // std::cout << "Kokkos::initialize() fyi: Pthread enabled but not
    // initialized" << std::endl ;
  }
#endif

#if defined(KOKKOS_ENABLE_HPX)
  if (std::is_same<Kokkos::Experimental::HPX,
                   Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Experimental::HPX,
                   Kokkos::HostSpace::execution_space>::value) {
    if (num_threads > 0) {
      Kokkos::Experimental::HPX::impl_initialize(num_threads);
    } else {
      Kokkos::Experimental::HPX::impl_initialize();
    }
    // std::cout << "Kokkos::initialize() fyi: HPX enabled and initialized" <<
    // std::endl ;
  } else {
    // std::cout << "Kokkos::initialize() fyi: HPX enabled but not initialized"
    // << std::endl ;
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  // Prevent "unused variable" warning for 'args' input struct.  If
  // Serial::initialize() ever needs to take arguments from the input
  // struct, you may remove this line of code.
  (void)args;

  // Always initialize Serial if it is configure time enabled
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  Kokkos::Serial::initialize();
#else
  Kokkos::Serial::impl_initialize();
#endif
#endif

#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  if (Impl::is_same<Kokkos::Experimental::OpenMPTarget,
                    Kokkos::DefaultExecutionSpace>::value) {
    if (num_threads > 0) {
      if (use_numa > 0) {
        Kokkos::Experimental::OpenMPTarget::initialize(num_threads, use_numa);
      } else {
        Kokkos::Experimental::OpenMPTarget::initialize(num_threads);
      }
    } else {
      Kokkos::Experimental::OpenMPTarget::initialize();
    }
    // std::cout << "Kokkos::initialize() fyi: OpenMP enabled and initialized"
    // << std::endl ;
  } else {
    // std::cout << "Kokkos::initialize() fyi: OpenMP enabled but not
    // initialized" << std::endl ;
  }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<Kokkos::Cuda, Kokkos::DefaultExecutionSpace>::value ||
      0 < use_gpu) {
    if (use_gpu > -1) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      Kokkos::Cuda::initialize(Kokkos::Cuda::SelectDevice(use_gpu));
#else
      Kokkos::Cuda::impl_initialize(Kokkos::Cuda::SelectDevice(use_gpu));
#endif
    } else {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      Kokkos::Cuda::initialize();
#else
      Kokkos::Cuda::impl_initialize();
#endif
    }
    // std::cout << "Kokkos::initialize() fyi: Cuda enabled and initialized" <<
    // std::endl ;
  }
#endif

#if defined(KOKKOS_ENABLE_ROCM)
  if (std::is_same<Kokkos::Experimental::ROCm,
                   Kokkos::DefaultExecutionSpace>::value ||
      0 < use_gpu) {
    if (use_gpu > -1) {
      Kokkos::Experimental::ROCm::initialize(
          Kokkos::Experimental::ROCm::SelectDevice(use_gpu));
    } else {
      Kokkos::Experimental::ROCm::initialize();
    }
    std::cout << "Kokkos::initialize() fyi: ROCm enabled and initialized"
              << std::endl;
  }
#endif

#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::initialize();
#else
  if (getenv("KOKKOS_PROFILE_LIBRARY") != nullptr) {
    std::cerr << "Kokkos::initialize() warning: Requested Kokkos Profiling, "
                 "but Kokkos was built without Profiling support"
              << std::endl;
  }
#endif
  g_is_initialized = true;
}

void finalize_internal(const bool all_spaces = false) {
  typename decltype(finalize_hooks)::size_type numSuccessfulCalls = 0;
  while (!finalize_hooks.empty()) {
    auto f = finalize_hooks.top();
    try {
      f();
    } catch (...) {
      std::cerr << "Kokkos::finalize: A finalize hook (set via "
                   "Kokkos::push_finalize_hook) threw an exception that it did "
                   "not catch."
                   "  Per std::atexit rules, this results in std::terminate.  "
                   "This is "
                   "finalize hook number "
                << numSuccessfulCalls
                << " (1-based indexing) "
                   "out of "
                << finalize_hooks.size()
                << " to call.  Remember that "
                   "Kokkos::finalize calls finalize hooks in reverse order "
                   "from how they "
                   "were pushed."
                << std::endl;
      std::terminate();
    }
    finalize_hooks.pop();
    ++numSuccessfulCalls;
  }

#if defined(KOKKOS_ENABLE_PROFILING)
  Kokkos::Profiling::finalize();
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<Kokkos::Cuda, Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    if (Kokkos::Cuda::is_initialized()) Kokkos::Cuda::finalize();
#else
    if (Kokkos::Cuda::impl_is_initialized()) Kokkos::Cuda::impl_finalize();
#endif
  }
#else
  (void)all_spaces;
#endif

#if defined(KOKKOS_ENABLE_ROCM)
  if (std::is_same<Kokkos::Experimental::ROCm,
                   Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
    if (Kokkos::Experimental::ROCm::is_initialized())
      Kokkos::Experimental::ROCm::finalize();
  }
#endif

#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  if (std::is_same<Kokkos::Experimental::OpenMPTarget,
                   Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
    if (Kokkos::Experimental::OpenMPTarget::is_initialized())
      Kokkos::Experimental::OpenMPTarget::finalize();
  }
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same<Kokkos::OpenMP, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::OpenMP, Kokkos::HostSpace::execution_space>::value ||
      all_spaces) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    if (Kokkos::OpenMP::is_initialized()) Kokkos::OpenMP::finalize();
#else
    if (Kokkos::OpenMP::impl_is_initialized()) Kokkos::OpenMP::impl_finalize();
#endif
  }
#endif

#if defined(KOKKOS_ENABLE_HPX)
  if (std::is_same<Kokkos::Experimental::HPX,
                   Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Experimental::HPX,
                   Kokkos::HostSpace::execution_space>::value ||
      all_spaces) {
    if (Kokkos::Experimental::HPX::impl_is_initialized())
      Kokkos::Experimental::HPX::impl_finalize();
  }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (std::is_same<Kokkos::Threads, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Threads,
                   Kokkos::HostSpace::execution_space>::value ||
      all_spaces) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    if (Kokkos::Threads::is_initialized()) Kokkos::Threads::finalize();
#else
    if (Kokkos::Threads::impl_is_initialized())
      Kokkos::Threads::impl_finalize();
#endif
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  if (Kokkos::Serial::is_initialized()) Kokkos::Serial::finalize();
#else
  if (Kokkos::Serial::impl_is_initialized()) Kokkos::Serial::impl_finalize();
#endif
#endif

  g_is_initialized = false;
  g_show_warnings  = true;
}

void fence_internal() {
#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<Kokkos::Cuda, Kokkos::DefaultExecutionSpace>::value) {
    Kokkos::Cuda::impl_static_fence();
  }
#endif

#if defined(KOKKOS_ENABLE_ROCM)
  if (std::is_same<Kokkos::Experimental::ROCm,
                   Kokkos::DefaultExecutionSpace>::value) {
    Kokkos::Experimental::ROCm().fence();
  }
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same<Kokkos::OpenMP, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::OpenMP, Kokkos::HostSpace::execution_space>::value) {
    Kokkos::OpenMP::impl_static_fence();
  }
#endif

#if defined(KOKKOS_ENABLE_HPX)
  Kokkos::Experimental::HPX::impl_static_fence();
#endif

#if defined(KOKKOS_ENABLE_THREADS)
  if (std::is_same<Kokkos::Threads, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Threads,
                   Kokkos::HostSpace::execution_space>::value) {
    Kokkos::Threads::impl_static_fence();
  }
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
  if (std::is_same<Kokkos::Serial, Kokkos::DefaultExecutionSpace>::value ||
      std::is_same<Kokkos::Serial, Kokkos::HostSpace::execution_space>::value) {
    Kokkos::Serial::impl_static_fence();
  }
#endif
}

bool check_arg(char const* arg, char const* expected) {
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  if (arg_len < exp_len) return false;
  if (std::strncmp(arg, expected, exp_len) != 0) return false;
  if (arg_len == exp_len) return true;
  /* if expected is "--threads", ignore "--threads-for-application"
     by checking this character          ---------^
     to see if it continues to make a longer name */
  if (std::isalnum(arg[exp_len]) || arg[exp_len] == '-' ||
      arg[exp_len] == '_') {
    return false;
  }
  return true;
}

bool check_int_arg(char const* arg, char const* expected, int* value) {
  if (!check_arg(arg, expected)) return false;
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  bool okay           = true;
  if (arg_len == exp_len || arg[exp_len] != '=') okay = false;
  char const* number = arg + exp_len + 1;
  if (!Impl::is_unsigned_int(number) || strlen(number) == 0) okay = false;
  *value = std::atoi(number);
  if (!okay) {
    std::ostringstream ss;
    ss << "Error: expecting an '=INT' after command line argument '" << expected
       << "'";
    ss << ". Raised by Kokkos::initialize(int narg, char* argc[]).";
    Impl::throw_runtime_exception(ss.str());
  }
  return true;
}

}  // namespace

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

void initialize(int& narg, char* arg[]) {
  int num_threads       = -1;
  int numa              = -1;
  int device            = -1;
  int ndevices          = -1;
  int skip_device       = 9999;
  bool disable_warnings = false;

  int kokkos_threads_found  = 0;
  int kokkos_numa_found     = 0;
  int kokkos_device_found   = 0;
  int kokkos_ndevices_found = 0;

  int iarg = 0;

  while (iarg < narg) {
    if (Impl::check_int_arg(arg[iarg], "--kokkos-threads", &num_threads)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_threads_found = 1;
      narg--;
    } else if (!kokkos_threads_found &&
               Impl::check_int_arg(arg[iarg], "--threads", &num_threads)) {
      iarg++;
    } else if (Impl::check_int_arg(arg[iarg], "--kokkos-numa", &numa)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_numa_found = 1;
      narg--;
    } else if (!kokkos_numa_found &&
               Impl::check_int_arg(arg[iarg], "--numa", &numa)) {
      iarg++;
    } else if (Impl::check_int_arg(arg[iarg], "--kokkos-device", &device)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_device_found = 1;
      narg--;
    } else if (!kokkos_device_found &&
               Impl::check_int_arg(arg[iarg], "--device", &device)) {
      iarg++;
    } else if (Impl::check_arg(arg[iarg], "--kokkos-ndevices") ||
               Impl::check_arg(arg[iarg], "--ndevices")) {
      // Find the number of device (expecting --device=XX)
      if (!((strncmp(arg[iarg], "--kokkos-ndevices=", 18) == 0) ||
            (strncmp(arg[iarg], "--ndevices=", 11) == 0)))
        Impl::throw_runtime_exception(
            "Error: expecting an '=INT[,INT]' after command line argument "
            "'--ndevices/--kokkos-ndevices'. Raised by Kokkos::initialize(int "
            "narg, char* argc[]).");

      char* num1      = strchr(arg[iarg], '=') + 1;
      char* num2      = strpbrk(num1, ",");
      int num1_len    = num2 == NULL ? strlen(num1) : num2 - num1;
      char* num1_only = new char[num1_len + 1];
      strncpy(num1_only, num1, num1_len);
      num1_only[num1_len] = 0;

      if (!Impl::is_unsigned_int(num1_only) || (strlen(num1_only) == 0)) {
        Impl::throw_runtime_exception(
            "Error: expecting an integer number after command line argument "
            "'--kokkos-ndevices'. Raised by Kokkos::initialize(int narg, char* "
            "argc[]).");
      }
      if ((strncmp(arg[iarg], "--kokkos-ndevices", 17) == 0) ||
          !kokkos_ndevices_found)
        ndevices = atoi(num1_only);
      delete[] num1_only;

      if (num2 != NULL) {
        if ((!Impl::is_unsigned_int(num2 + 1)) || (strlen(num2) == 1))
          Impl::throw_runtime_exception(
              "Error: expecting an integer number after command line argument "
              "'--kokkos-ndevices=XX,'. Raised by Kokkos::initialize(int narg, "
              "char* argc[]).");

        if ((strncmp(arg[iarg], "--kokkos-ndevices", 17) == 0) ||
            !kokkos_ndevices_found)
          skip_device = atoi(num2 + 1);
      }

      // Remove the --kokkos-ndevices argument from the list but leave
      // --ndevices
      if (strncmp(arg[iarg], "--kokkos-ndevices", 17) == 0) {
        for (int k = iarg; k < narg - 1; k++) {
          arg[k] = arg[k + 1];
        }
        kokkos_ndevices_found = 1;
        narg--;
      } else {
        iarg++;
      }
    } else if (strcmp(arg[iarg], "--kokkos-disable-warnings") == 0) {
      disable_warnings = true;
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
    } else if ((strcmp(arg[iarg], "--kokkos-help") == 0) ||
               (strcmp(arg[iarg], "--help") == 0)) {
      std::cout << std::endl;
      std::cout << "-----------------------------------------------------------"
                   "---------------------"
                << std::endl;
      std::cout << "-------------Kokkos command line "
                   "arguments--------------------------------------"
                << std::endl;
      std::cout << "-----------------------------------------------------------"
                   "---------------------"
                << std::endl;
      std::cout << "The following arguments exist also without prefix 'kokkos' "
                   "(e.g. --help)."
                << std::endl;
      std::cout << "The prefixed arguments will be removed from the list by "
                   "Kokkos::initialize(),"
                << std::endl;
      std::cout << "the non-prefixed ones are not removed. Prefixed versions "
                   "take precedence over "
                << std::endl;
      std::cout << "non prefixed ones, and the last occurrence of an argument "
                   "overwrites prior"
                << std::endl;
      std::cout << "settings." << std::endl;
      std::cout << std::endl;
      std::cout << "--kokkos-help               : print this message"
                << std::endl;
      std::cout
          << "--kokkos-disable-warnings   : disable kokkos warning messages"
          << std::endl;
      std::cout
          << "--kokkos-threads=INT        : specify total number of threads or"
          << std::endl;
      std::cout << "                              number of threads per NUMA "
                   "region if "
                << std::endl;
      std::cout << "                              used in conjunction with "
                   "'--numa' option. "
                << std::endl;
      std::cout << "--kokkos-numa=INT           : specify number of NUMA "
                   "regions used by process."
                << std::endl;
      std::cout << "--kokkos-device=INT         : specify device id to be used "
                   "by Kokkos. "
                << std::endl;
      std::cout << "--kokkos-ndevices=INT[,INT] : used when running MPI jobs. "
                   "Specify number of"
                << std::endl;
      std::cout << "                              devices per node to be used. "
                   "Process to device"
                << std::endl;
      std::cout << "                              mapping happens by obtaining "
                   "the local MPI rank"
                << std::endl;
      std::cout << "                              and assigning devices "
                   "round-robin. The optional"
                << std::endl;
      std::cout << "                              second argument allows for "
                   "an existing device"
                << std::endl;
      std::cout << "                              to be ignored. This is most "
                   "useful on workstations"
                << std::endl;
      std::cout << "                              with multiple GPUs of which "
                   "one is used to drive"
                << std::endl;
      std::cout << "                              screen output." << std::endl;
      std::cout << std::endl;
      std::cout << "-----------------------------------------------------------"
                   "---------------------"
                << std::endl;
      std::cout << std::endl;

      // Remove the --kokkos-help argument from the list but leave --ndevices
      if (strcmp(arg[iarg], "--kokkos-help") == 0) {
        for (int k = iarg; k < narg - 1; k++) {
          arg[k] = arg[k + 1];
        }
        narg--;
      } else {
        iarg++;
      }
    } else
      iarg++;
  }

  // Read environment variables
  char* endptr;
  auto env_num_threads_str = std::getenv("KOKKOS_NUM_THREADS");
  if (env_num_threads_str != nullptr) {
    errno                = 0;
    auto env_num_threads = std::strtol(env_num_threads_str, &endptr, 10);
    if (endptr == env_num_threads_str)
      Impl::throw_runtime_exception(
          "Error: cannot convert KOKKOS_NUM_THREADS to an integer. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
    if (errno == ERANGE)
      Impl::throw_runtime_exception(
          "Error: KOKKOS_NUM_THREADS out of range of representable values by "
          "an integer. Raised by Kokkos::initialize(int narg, char* argc[]).");
    if ((num_threads != -1) && (env_num_threads != num_threads))
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-threads and "
          "KOKKOS_NUM_THREADS if both are set. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
    else
      num_threads = env_num_threads;
  }
  auto env_numa_str = std::getenv("KOKKOS_NUMA");
  if (env_numa_str != nullptr) {
    errno         = 0;
    auto env_numa = std::strtol(env_numa_str, &endptr, 10);
    if (endptr == env_numa_str)
      Impl::throw_runtime_exception(
          "Error: cannot convert KOKKOS_NUMA to an integer. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
    if (errno == ERANGE)
      Impl::throw_runtime_exception(
          "Error: KOKKOS_NUMA out of range of representable values by an "
          "integer. Raised by Kokkos::initialize(int narg, char* argc[]).");
    if ((numa != -1) && (env_numa != numa))
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-numa and KOKKOS_NUMA if "
          "both are set. Raised by Kokkos::initialize(int narg, char* "
          "argc[]).");
    else
      numa = env_numa;
  }
  auto env_device_str = std::getenv("KOKKOS_DEVICE_ID");
  if (env_device_str != nullptr) {
    errno           = 0;
    auto env_device = std::strtol(env_device_str, &endptr, 10);
    if (endptr == env_device_str)
      Impl::throw_runtime_exception(
          "Error: cannot convert KOKKOS_DEVICE_ID to an integer. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
    if (errno == ERANGE)
      Impl::throw_runtime_exception(
          "Error: KOKKOS_DEVICE_ID out of range of representable values by an "
          "integer. Raised by Kokkos::initialize(int narg, char* argc[]).");
    if ((device != -1) && (env_device != device))
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-device and "
          "KOKKOS_DEVICE_ID if both are set. Raised by Kokkos::initialize(int "
          "narg, char* argc[]).");
    else
      device = env_device;
  }
  auto env_rdevices_str = std::getenv("KOKKOS_RAND_DEVICES");
  auto env_ndevices_str = std::getenv("KOKKOS_NUM_DEVICES");
  if (env_ndevices_str != nullptr || env_rdevices_str != nullptr) {
    errno = 0;
    if (env_ndevices_str != nullptr && env_rdevices_str != nullptr) {
      Impl::throw_runtime_exception(
          "Error: cannot specify both KOKKOS_NUM_DEVICES and "
          "KOKKOS_RAND_DEVICES. "
          "Raised by Kokkos::initialize(int narg, char* argc[]).");
    }
    int rdevices = -1;
    if (env_ndevices_str != nullptr) {
      auto env_ndevices = std::strtol(env_ndevices_str, &endptr, 10);
      if (endptr == env_ndevices_str)
        Impl::throw_runtime_exception(
            "Error: cannot convert KOKKOS_NUM_DEVICES to an integer. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");
      if (errno == ERANGE)
        Impl::throw_runtime_exception(
            "Error: KOKKOS_NUM_DEVICES out of range of representable values by "
            "an integer. Raised by Kokkos::initialize(int narg, char* "
            "argc[]).");
      if ((ndevices != -1) && (env_ndevices != ndevices))
        Impl::throw_runtime_exception(
            "Error: expecting a match between --kokkos-ndevices and "
            "KOKKOS_NUM_DEVICES if both are set. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");
      else
        ndevices = env_ndevices;
    } else {  // you set KOKKOS_RAND_DEVICES
      auto env_rdevices = std::strtol(env_rdevices_str, &endptr, 10);
      if (endptr == env_ndevices_str)
        Impl::throw_runtime_exception(
            "Error: cannot convert KOKKOS_RAND_DEVICES to an integer. Raised "
            "by Kokkos::initialize(int narg, char* argc[]).");
      if (errno == ERANGE)
        Impl::throw_runtime_exception(
            "Error: KOKKOS_RAND_DEVICES out of range of representable values "
            "by an integer. Raised by Kokkos::initialize(int narg, char* "
            "argc[]).");
      else
        rdevices = env_rdevices;
    }
    // Skip device
    auto env_skip_device_str = std::getenv("KOKKOS_SKIP_DEVICE");
    if (env_skip_device_str != nullptr) {
      errno                = 0;
      auto env_skip_device = std::strtol(env_skip_device_str, &endptr, 10);
      if (endptr == env_skip_device_str)
        Impl::throw_runtime_exception(
            "Error: cannot convert KOKKOS_SKIP_DEVICE to an integer. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");
      if (errno == ERANGE)
        Impl::throw_runtime_exception(
            "Error: KOKKOS_SKIP_DEVICE out of range of representable values by "
            "an integer. Raised by Kokkos::initialize(int narg, char* "
            "argc[]).");
      if ((skip_device != 9999) && (env_skip_device != skip_device))
        Impl::throw_runtime_exception(
            "Error: expecting a match between --kokkos-ndevices and "
            "KOKKOS_SKIP_DEVICE if both are set. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");
      else
        skip_device = env_skip_device;
    }
    if (rdevices > 0) {
      if (skip_device > 0 && rdevices == 1)
        Impl::throw_runtime_exception(
            "Error: cannot KOKKOS_SKIP_DEVICE the only KOKKOS_RAND_DEVICE. "
            "Raised by Kokkos::initialize(int narg, char* argc[]).");

      std::srand(getpid());
      while (device < 0) {
        int test_device = std::rand() % rdevices;
        if (test_device != skip_device) device = test_device;
      }
    }
  }
  char* env_disablewarnings_str = std::getenv("KOKKOS_DISABLE_WARNINGS");
  if (env_disablewarnings_str != nullptr) {
    std::string env_str(env_disablewarnings_str);  // deep-copies string
    for (char& c : env_str) {
      c = toupper(c);
    }
    if ((env_str == "TRUE") || (env_str == "ON") || (env_str == "1"))
      disable_warnings = true;
    else if (disable_warnings)
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-disable-warnings and "
          "KOKKOS_DISABLE_WARNINGS if both are set. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
  }

  InitArguments arguments;
  arguments.num_threads      = num_threads;
  arguments.num_numa         = numa;
  arguments.device_id        = device;
  arguments.ndevices         = ndevices;
  arguments.skip_device      = skip_device;
  arguments.disable_warnings = disable_warnings;
  Impl::initialize_internal(arguments);
}

void initialize(const InitArguments& arguments) {
  Impl::initialize_internal(arguments);
}

void push_finalize_hook(std::function<void()> f) { finalize_hooks.push(f); }

void finalize() { Impl::finalize_internal(); }

void finalize_all() {
  enum { all_spaces = true };
  Impl::finalize_internal(all_spaces);
}

void fence() { Impl::fence_internal(); }

void print_configuration(std::ostream& out, const bool detail) {
  std::ostringstream msg;

  msg << "Compiler:" << std::endl;
#ifdef KOKKOS_COMPILER_APPLECC
  msg << "  KOKKOS_COMPILER_APPLECC: " << KOKKOS_COMPILER_APPLECC << std::endl;
#endif
#ifdef KOKKOS_COMPILER_CLANG
  msg << "  KOKKOS_COMPILER_CLANG: " << KOKKOS_COMPILER_CLANG << std::endl;
#endif
#ifdef KOKKOS_COMPILER_CRAYC
  msg << "  KOKKOS_COMPILER_CRAYC: " << KOKKOS_COMPILER_CRAYC << std::endl;
#endif
#ifdef KOKKOS_COMPILER_GNU
  msg << "  KOKKOS_COMPILER_GNU: " << KOKKOS_COMPILER_GNU << std::endl;
#endif
#ifdef KOKKOS_COMPILER_IBM
  msg << "  KOKKOS_COMPILER_IBM: " << KOKKOS_COMPILER_IBM << std::endl;
#endif
#ifdef KOKKOS_COMPILER_INTEL
  msg << "  KOKKOS_COMPILER_INTEL: " << KOKKOS_COMPILER_INTEL << std::endl;
#endif
#ifdef KOKKOS_COMPILER_NVCC
  msg << "  KOKKOS_COMPILER_NVCC: " << KOKKOS_COMPILER_NVCC << std::endl;
#endif
#ifdef KOKKOS_COMPILER_PGI
  msg << "  KOKKOS_COMPILER_PGI: " << KOKKOS_COMPILER_PGI << std::endl;
#endif

  msg << "Architecture:" << std::endl;
#ifdef KOKKOS_ENABLE_ISA_KNC
  msg << "  KOKKOS_ENABLE_ISA_KNC: yes" << std::endl;
#else
  msg << "  KOKKOS_ENABLE_ISA_KNC: no" << std::endl;
#endif
#ifdef KOKKOS_ENABLE_ISA_POWERPCLE
  msg << "  KOKKOS_ENABLE_ISA_POWERPCLE: yes" << std::endl;
#else
  msg << "  KOKKOS_ENABLE_ISA_POWERPCLE: no" << std::endl;
#endif
#ifdef KOKKOS_ENABLE_ISA_X86_64
  msg << "  KOKKOS_ENABLE_ISA_X86_64: yes" << std::endl;
#else
  msg << "  KOKKOS_ENABLE_ISA_X86_64: no" << std::endl;
#endif

  msg << "Devices:" << std::endl;
  msg << "  KOKKOS_ENABLE_CUDA: ";
#ifdef KOKKOS_ENABLE_CUDA
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_OPENMP: ";
#ifdef KOKKOS_ENABLE_OPENMP
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_HPX: ";
#ifdef KOKKOS_ENABLE_HPX
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_THREADS: ";
#ifdef KOKKOS_ENABLE_THREADS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_QTHREADS: ";
#ifdef KOKKOS_ENABLE_QTHREADS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_SERIAL: ";
#ifdef KOKKOS_ENABLE_SERIAL
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "Default Device:" << std::endl;
  msg << "  KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA: ";
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP: ";
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS: ";
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS: ";
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL: ";
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "Atomics:" << std::endl;
  msg << "  KOKKOS_ENABLE_CUDA_ATOMICS: ";
#ifdef KOKKOS_ENABLE_CUDA_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_GNU_ATOMICS: ";
#ifdef KOKKOS_ENABLE_GNU_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_INTEL_ATOMICS: ";
#ifdef KOKKOS_ENABLE_INTEL_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_OPENMP_ATOMICS: ";
#ifdef KOKKOS_ENABLE_OPENMP_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_WINDOWS_ATOMICS: ";
#ifdef KOKKOS_ENABLE_WINDOWS_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_SERIAL_ATOMICS: ";
#ifdef KOKKOS_ENABLE_SERIAL_ATOMICS
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "Vectorization:" << std::endl;
  msg << "  KOKKOS_ENABLE_PRAGMA_IVDEP: ";
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_PRAGMA_LOOPCOUNT: ";
#ifdef KOKKOS_ENABLE_PRAGMA_LOOPCOUNT
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_PRAGMA_SIMD: ";
#ifdef KOKKOS_ENABLE_PRAGMA_SIMD
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_PRAGMA_UNROLL: ";
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_PRAGMA_VECTOR: ";
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "Memory:" << std::endl;
  msg << "  KOKKOS_ENABLE_HBWSPACE: ";
#ifdef KOKKOS_ENABLE_HBWSPACE
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_INTEL_MM_ALLOC: ";
#ifdef KOKKOS_ENABLE_INTEL_MM_ALLOC
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_POSIX_MEMALIGN: ";
#ifdef KOKKOS_ENABLE_POSIX_MEMALIGN
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "Options:" << std::endl;
  msg << "  KOKKOS_ENABLE_ASM: ";
#ifdef KOKKOS_ENABLE_ASM
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CXX14: ";
#ifdef KOKKOS_ENABLE_CXX14
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CXX17: ";
#ifdef KOKKOS_ENABLE_CXX17
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CXX20: ";
#ifdef KOKKOS_ENABLE_CXX20
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK: ";
#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_HWLOC: ";
#ifdef KOKKOS_ENABLE_HWLOC
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_LIBRT: ";
#ifdef KOKKOS_ENABLE_LIBRT
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_MPI: ";
#ifdef KOKKOS_ENABLE_MPI
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_PROFILING: ";
#ifdef KOKKOS_ENABLE_PROFILING
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

#ifdef KOKKOS_ENABLE_CUDA
  msg << "Cuda Options:" << std::endl;
  msg << "  KOKKOS_ENABLE_CUDA_LAMBDA: ";
#ifdef KOKKOS_ENABLE_CUDA_LAMBDA
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CUDA_LDG_INTRINSIC: ";
#ifdef KOKKOS_ENABLE_CUDA_LDG_INTRINSIC
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE: ";
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CUDA_UVM: ";
#ifdef KOKKOS_ENABLE_CUDA_UVM
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CUSPARSE: ";
#ifdef KOKKOS_ENABLE_CUSPARSE
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif
  msg << "  KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA: ";
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

#endif

  msg << "\nRuntime Configuration:" << std::endl;
#ifdef KOKKOS_ENABLE_CUDA
  Cuda::print_configuration(msg, detail);
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  OpenMP::print_configuration(msg, detail);
#endif
#ifdef KOKKOS_ENABLE_HPX
  Experimental::HPX::print_configuration(msg, detail);
#endif
#if defined(KOKKOS_ENABLE_THREADS)
  Threads::print_configuration(msg, detail);
#endif
#ifdef KOKKOS_ENABLE_QTHREADS
  Qthreads::print_configuration(msg, detail);
#endif
#ifdef KOKKOS_ENABLE_SERIAL
  Serial::print_configuration(msg, detail);
#endif

  out << msg.str() << std::endl;
}

bool is_initialized() noexcept { return g_is_initialized; }

bool show_warnings() noexcept { return g_show_warnings; }

#ifdef KOKKOS_COMPILER_PGI
namespace Impl {
// Bizzarely, an extra jump instruction forces the PGI compiler to not have a
// bug related to (probably?) empty base optimization and/or aggregate
// construction.
void _kokkos_pgi_compiler_bug_workaround() {}
}  // end namespace Impl
#endif

}  // namespace Kokkos
