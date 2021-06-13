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
#include <impl/Kokkos_ExecSpaceInitializer.hpp>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stack>
#include <functional>
#include <list>
#include <cerrno>
#include <regex>
#ifndef _WIN32
#include <unistd.h>
#else
#include <Windows.h>
#endif

//----------------------------------------------------------------------------
namespace {
bool g_is_initialized = false;
bool g_show_warnings  = true;
bool g_tune_internals = false;
// When compiling with clang/LLVM and using the GNU (GCC) C++ Standard Library
// (any recent version between GCC 7.3 and GCC 9.2), std::deque SEGV's during
// the unwinding of the atexit(3C) handlers at program termination.  However,
// this bug is not observable when building with GCC.
// As an added bonus, std::list<T> provides constant insertion and
// deletion time complexity, which translates to better run-time performance. As
// opposed to std::deque<T> which does not provide the same constant time
// complexity for inserts/removals, since std::deque<T> is implemented as a
// segmented array.
using hook_function_type = std::function<void()>;
std::stack<hook_function_type, std::list<hook_function_type>> finalize_hooks;
}  // namespace

namespace Kokkos {
namespace Impl {
/**
 * The category is only used in printing, tools
 * get all metadata free of category
 */
using metadata_category_type = std::string;
using metadata_key_type      = std::string;
using metadata_value_type    = std::string;

std::map<metadata_category_type,
         std::map<metadata_key_type, metadata_value_type>>
    metadata_map;

void declare_configuration_metadata(const std::string& category,
                                    const std::string& key,
                                    const std::string& value) {
  metadata_map[category][key] = value;
}

ExecSpaceManager& ExecSpaceManager::get_instance() {
  static ExecSpaceManager space_initializer = {};
  return space_initializer;
}

void ExecSpaceManager::register_space_factory(
    const std::string name, std::unique_ptr<ExecSpaceInitializerBase> space) {
  exec_space_factory_list[name] = std::move(space);
}

void ExecSpaceManager::initialize_spaces(const Kokkos::InitArguments& args) {
  // Note: the names of the execution spaces, used as keys in the map, encode
  // the ordering of the initialization code from the old initializtion stuff.
  // Eventually, we may want to do something less brittle than this, but for now
  // we're just preserving compatibility with the old implementation.
  for (auto& to_init : exec_space_factory_list) {
    to_init.second->initialize(args);
  }
}

void ExecSpaceManager::finalize_spaces(const bool all_spaces) {
  for (auto& to_finalize : exec_space_factory_list) {
    to_finalize.second->finalize(all_spaces);
  }
}

void ExecSpaceManager::static_fence() {
  for (auto& to_fence : exec_space_factory_list) {
    to_fence.second->fence();
  }
}
void ExecSpaceManager::print_configuration(std::ostream& msg,
                                           const bool detail) {
  for (auto& to_print : exec_space_factory_list) {
    to_print.second->print_configuration(msg, detail);
  }
}

int get_ctest_gpu(const char* local_rank_str) {
  auto const* ctest_kokkos_device_type =
      std::getenv("CTEST_KOKKOS_DEVICE_TYPE");
  if (!ctest_kokkos_device_type) {
    return 0;
  }

  auto const* ctest_resource_group_count_str =
      std::getenv("CTEST_RESOURCE_GROUP_COUNT");
  if (!ctest_resource_group_count_str) {
    return 0;
  }

  // Make sure rank is within bounds of resource groups specified by CTest
  auto resource_group_count = std::stoi(ctest_resource_group_count_str);
  auto local_rank           = std::stoi(local_rank_str);
  if (local_rank >= resource_group_count) {
    std::ostringstream ss;
    ss << "Error: local rank " << local_rank
       << " is outside the bounds of resource groups provided by CTest. Raised"
       << " by Kokkos::Impl::get_ctest_gpu().";
    throw_runtime_exception(ss.str());
  }

  // Get the resource types allocated to this resource group
  std::ostringstream ctest_resource_group;
  ctest_resource_group << "CTEST_RESOURCE_GROUP_" << local_rank;
  std::string ctest_resource_group_name = ctest_resource_group.str();
  auto const* ctest_resource_group_str =
      std::getenv(ctest_resource_group_name.c_str());
  if (!ctest_resource_group_str) {
    std::ostringstream ss;
    ss << "Error: " << ctest_resource_group_name << " is not specified. Raised"
       << " by Kokkos::Impl::get_ctest_gpu().";
    throw_runtime_exception(ss.str());
  }

  // Look for the device type specified in CTEST_KOKKOS_DEVICE_TYPE
  bool found_device                        = false;
  std::string ctest_resource_group_cxx_str = ctest_resource_group_str;
  std::istringstream instream(ctest_resource_group_cxx_str);
  while (true) {
    std::string devName;
    std::getline(instream, devName, ',');
    if (devName == ctest_kokkos_device_type) {
      found_device = true;
      break;
    }
    if (instream.eof() || devName.length() == 0) {
      break;
    }
  }

  if (!found_device) {
    std::ostringstream ss;
    ss << "Error: device type '" << ctest_kokkos_device_type
       << "' not included in " << ctest_resource_group_name
       << ". Raised by Kokkos::Impl::get_ctest_gpu().";
    throw_runtime_exception(ss.str());
  }

  // Get the device ID
  std::string ctest_device_type_upper = ctest_kokkos_device_type;
  for (auto& c : ctest_device_type_upper) {
    c = std::toupper(c);
  }
  ctest_resource_group << "_" << ctest_device_type_upper;

  std::string ctest_resource_group_id_name = ctest_resource_group.str();
  auto resource_str = std::getenv(ctest_resource_group_id_name.c_str());
  if (!resource_str) {
    std::ostringstream ss;
    ss << "Error: " << ctest_resource_group_id_name
       << " is not specified. Raised by Kokkos::Impl::get_ctest_gpu().";
    throw_runtime_exception(ss.str());
  }

  auto const* comma = std::strchr(resource_str, ',');
  if (!comma || strncmp(resource_str, "id:", 3)) {
    std::ostringstream ss;
    ss << "Error: invalid value of " << ctest_resource_group_id_name << ": '"
       << resource_str << "'. Raised by Kokkos::Impl::get_ctest_gpu().";
    throw_runtime_exception(ss.str());
  }

  std::string id(resource_str + 3, comma - resource_str - 3);
  return std::stoi(id.c_str());
}

// function to extract gpu # from args
int get_gpu(const InitArguments& args) {
  int use_gpu        = args.device_id;
  const int ndevices = [](int num_devices) -> int {
    if (num_devices > 0) return num_devices;
#if defined(KOKKOS_ENABLE_CUDA)
    return Cuda::detect_device_count();
#elif defined(KOKKOS_ENABLE_HIP)
    return Experimental::HIP::detect_device_count();
#elif defined(KOKKOS_ENABLE_SYCL)
    return sycl::device::get_devices(sycl::info::device_type::gpu).size();
#else
    return num_devices;
#endif
  }(args.ndevices);
  const int skip_device = args.skip_device;

  // if the exact device is not set, but ndevices was given, assign round-robin
  // using on-node MPI rank
  if (use_gpu < 0) {
    auto const* local_rank_str =
        std::getenv("OMPI_COMM_WORLD_LOCAL_RANK");  // OpenMPI
    if (!local_rank_str)
      local_rank_str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK");  // MVAPICH2
    if (!local_rank_str)
      local_rank_str = std::getenv("SLURM_LOCALID");  // SLURM

    auto const* ctest_kokkos_device_type =
        std::getenv("CTEST_KOKKOS_DEVICE_TYPE");  // CTest
    auto const* ctest_resource_group_count_str =
        std::getenv("CTEST_RESOURCE_GROUP_COUNT");  // CTest
    if (ctest_kokkos_device_type && ctest_resource_group_count_str &&
        local_rank_str) {
      // Use the device assigned by CTest
      use_gpu = get_ctest_gpu(local_rank_str);
    } else if (ndevices > 0) {
      // Use the device assigned by the rank
      if (local_rank_str) {
        auto local_rank = std::stoi(local_rank_str);
        use_gpu         = local_rank % ndevices;
      } else {
        // user only gave use ndevices, but the MPI environment variable wasn't
        // set. start with GPU 0 at this point
        use_gpu = 0;
      }
    }
    // shift assignments over by one so no one is assigned to "skip_device"
    if (use_gpu >= skip_device) ++use_gpu;
  }
  return use_gpu;
}
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

void initialize_backends(const InitArguments& args) {
// This is an experimental setting
// For KNL in Flat mode this variable should be set, so that
// memkind allocates high bandwidth memory correctly.
#ifdef KOKKOS_ENABLE_HBWSPACE
  setenv("MEMKIND_HBW_NODES", "1", 0);
#endif

  Impl::ExecSpaceManager::get_instance().initialize_spaces(args);
}

void initialize_profiling(const InitArguments& args) {
  Kokkos::Profiling::initialize(args.tool_lib);
  if (args.tool_help) {
    if (!Kokkos::Tools::printHelp(args.tool_args)) {
      std::cerr << "Tool has not provided a help message" << std::endl;
    }
    g_is_initialized = true;
    ::Kokkos::finalize();
    std::exit(EXIT_SUCCESS);
  }
  Kokkos::Tools::parseArgs(args.tool_args);
  for (const auto& category_value : Kokkos::Impl::metadata_map) {
    for (const auto& key_value : category_value.second) {
      Kokkos::Tools::declareMetadata(key_value.first, key_value.second);
    }
  }
}

std::string version_string_from_int(int version_number) {
  std::stringstream str_builder;
  str_builder << version_number / 10000 << "." << (version_number % 10000) / 100
              << "." << version_number % 100;
  return str_builder.str();
}
void pre_initialize_internal(const InitArguments& args) {
  if (args.disable_warnings) g_show_warnings = false;
  if (args.tune_internals) g_tune_internals = true;
  declare_configuration_metadata("version_info", "Kokkos Version",
                                 version_string_from_int(KOKKOS_VERSION));
#ifdef KOKKOS_COMPILER_APPLECC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_APPLECC",
                                 std::to_string(KOKKOS_COMPILER_APPLECC));
  declare_configuration_metadata("tools_only", "compiler_family", "apple");
#endif
#ifdef KOKKOS_COMPILER_CLANG
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_CLANG",
                                 std::to_string(KOKKOS_COMPILER_CLANG));
  declare_configuration_metadata("tools_only", "compiler_family", "clang");
#endif
#ifdef KOKKOS_COMPILER_CRAYC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_CRAYC",
                                 std::to_string(KOKKOS_COMPILER_CRAYC));
  declare_configuration_metadata("tools_only", "compiler_family", "cray");
#endif
#ifdef KOKKOS_COMPILER_GNU
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_GNU",
                                 std::to_string(KOKKOS_COMPILER_GNU));
  declare_configuration_metadata("tools_only", "compiler_family", "gnu");
#endif
#ifdef KOKKOS_COMPILER_IBM
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_IBM",
                                 std::to_string(KOKKOS_COMPILER_IBM));
  declare_configuration_metadata("tools_only", "compiler_family", "ibm");
#endif
#ifdef KOKKOS_COMPILER_INTEL
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_INTEL",
                                 std::to_string(KOKKOS_COMPILER_INTEL));
  declare_configuration_metadata("tools_only", "compiler_family", "intel");
#endif
#ifdef KOKKOS_COMPILER_NVCC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_NVCC",
                                 std::to_string(KOKKOS_COMPILER_NVCC));
  declare_configuration_metadata("tools_only", "compiler_family", "nvcc");
#endif
#ifdef KOKKOS_COMPILER_PGI
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_PGI",
                                 std::to_string(KOKKOS_COMPILER_PGI));
  declare_configuration_metadata("tools_only", "compiler_family", "pgi");
#endif
#ifdef KOKKOS_COMPILER_MSVC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_MSVC",
                                 std::to_string(KOKKOS_COMPILER_MSVC));
  declare_configuration_metadata("tools_only", "compiler_family", "msvc");
#endif
#ifdef KOKKOS_ENABLE_ISA_KNC
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_KNC",
                                 "yes");
#else
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_KNC", "no");
#endif
#ifdef KOKKOS_ENABLE_ISA_POWERPCLE
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_POWERPCLE",
                                 "yes");
#else
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_POWERPCLE",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_ISA_X86_64
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_X86_64",
                                 "yes");
#else
  declare_configuration_metadata("architecture", "KOKKOS_ENABLE_ISA_X86_64",
                                 "no");
#endif

#ifdef KOKKOS_ENABLE_GNU_ATOMICS
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_GNU_ATOMICS", "yes");
#else
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_GNU_ATOMICS", "no");
#endif
#ifdef KOKKOS_ENABLE_INTEL_ATOMICS
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_INTEL_ATOMICS",
                                 "yes");
#else
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_INTEL_ATOMICS",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_WINDOWS_ATOMICS
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_WINDOWS_ATOMICS",
                                 "yes");
#else
  declare_configuration_metadata("atomics", "KOKKOS_ENABLE_WINDOWS_ATOMICS",
                                 "no");
#endif

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_IVDEP",
                                 "yes");
#else
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_IVDEP",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_LOOPCOUNT
  declare_configuration_metadata("vectorization",
                                 "KOKKOS_ENABLE_PRAGMA_LOOPCOUNT", "yes");
#else
  declare_configuration_metadata("vectorization",
                                 "KOKKOS_ENABLE_PRAGMA_LOOPCOUNT", "no");
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_SIMD
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_SIMD",
                                 "yes");
#else
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_SIMD",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_UNROLL",
                                 "yes");
#else
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_UNROLL",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_VECTOR",
                                 "yes");
#else
  declare_configuration_metadata("vectorization", "KOKKOS_ENABLE_PRAGMA_VECTOR",
                                 "no");
#endif

#ifdef KOKKOS_ENABLE_HBWSPACE
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_HBWSPACE", "yes");
#else
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_HBWSPACE", "no");
#endif
#ifdef KOKKOS_ENABLE_INTEL_MM_ALLOC
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_INTEL_MM_ALLOC",
                                 "yes");
#else
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_INTEL_MM_ALLOC",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_POSIX_MEMALIGN
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_POSIX_MEMALIGN",
                                 "yes");
#else
  declare_configuration_metadata("memory", "KOKKOS_ENABLE_POSIX_MEMALIGN",
                                 "no");
#endif

#ifdef KOKKOS_ENABLE_ASM
  declare_configuration_metadata("options", "KOKKOS_ENABLE_ASM", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_ASM", "no");
#endif
#ifdef KOKKOS_ENABLE_CXX14
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX14", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX14", "no");
#endif
#ifdef KOKKOS_ENABLE_CXX17
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX17", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX17", "no");
#endif
#ifdef KOKKOS_ENABLE_CXX20
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX20", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX20", "no");
#endif
#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  declare_configuration_metadata("options", "KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK",
                                 "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK",
                                 "no");
#endif
#ifdef KOKKOS_ENABLE_HWLOC
  declare_configuration_metadata("options", "KOKKOS_ENABLE_HWLOC", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_HWLOC", "no");
#endif
#ifdef KOKKOS_ENABLE_LIBRT
  declare_configuration_metadata("options", "KOKKOS_ENABLE_LIBRT", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_LIBRT", "no");
#endif
#ifdef KOKKOS_ENABLE_MPI
  declare_configuration_metadata("options", "KOKKOS_ENABLE_MPI", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_MPI", "no");
#endif
  declare_configuration_metadata("architecture", "Default Device",
                                 typeid(Kokkos::DefaultExecutionSpace).name());
}

void post_initialize_internal(const InitArguments& args) {
  initialize_profiling(args);
  g_is_initialized = true;
}

void initialize_internal(const InitArguments& args) {
  pre_initialize_internal(args);
  initialize_backends(args);
  post_initialize_internal(args);
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

  Kokkos::Profiling::finalize();

  Impl::ExecSpaceManager::get_instance().finalize_spaces(all_spaces);

  g_is_initialized = false;
  g_show_warnings  = true;
  g_tune_internals = false;
}

void fence_internal() { Impl::ExecSpaceManager::get_instance().static_fence(); }

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
  *value = std::stoi(number);
  if (!okay) {
    std::ostringstream ss;
    ss << "Error: expecting an '=INT' after command line argument '" << expected
       << "'";
    ss << ". Raised by Kokkos::initialize(int narg, char* argc[]).";
    Impl::throw_runtime_exception(ss.str());
  }
  return true;
}

bool check_str_arg(char const* arg, char const* expected, std::string& value) {
  if (!check_arg(arg, expected)) return false;
  std::size_t arg_len = std::strlen(arg);
  std::size_t exp_len = std::strlen(expected);
  bool okay           = true;
  if (arg_len == exp_len || arg[exp_len] != '=') okay = false;
  char const* remain = arg + exp_len + 1;
  value              = remain;
  if (!okay) {
    std::ostringstream ss;
    ss << "Error: expecting an '=STRING' after command line argument '"
       << expected << "'";
    ss << ". Raised by Kokkos::initialize(int narg, char* argc[]).";
    Impl::throw_runtime_exception(ss.str());
  }
  return true;
}

void warn_deprecated_command_line_argument(std::string deprecated,
                                           std::string valid) {
  std::cerr
      << "Warning: command line argument '" << deprecated
      << "' is deprecated. Use '" << valid
      << "' instead. Raised by Kokkos::initialize(int narg, char* argc[])."
      << std::endl;
}

unsigned get_process_id() {
#ifdef _WIN32
  return unsigned(GetCurrentProcessId());
#else
  return unsigned(getpid());
#endif
}

void parse_command_line_arguments(int& narg, char* arg[],
                                  InitArguments& arguments) {
  auto& num_threads      = arguments.num_threads;
  auto& numa             = arguments.num_numa;
  auto& device           = arguments.device_id;
  auto& ndevices         = arguments.ndevices;
  auto& skip_device      = arguments.skip_device;
  auto& disable_warnings = arguments.disable_warnings;
  auto& tune_internals   = arguments.tune_internals;
  auto& tool_help        = arguments.tool_help;
  auto& tool_args        = arguments.tool_args;
  auto& tool_lib         = arguments.tool_lib;

  bool kokkos_threads_found  = false;
  bool kokkos_numa_found     = false;
  bool kokkos_device_found   = false;
  bool kokkos_ndevices_found = false;

  int iarg = 0;

  while (iarg < narg) {
    if (check_int_arg(arg[iarg], "--kokkos-threads", &num_threads)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_threads_found = true;
      narg--;
    } else if (!kokkos_threads_found &&
               check_int_arg(arg[iarg], "--threads", &num_threads)) {
      iarg++;
    } else if (check_int_arg(arg[iarg], "--kokkos-numa", &numa)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_numa_found = true;
      narg--;
    } else if (!kokkos_numa_found &&
               check_int_arg(arg[iarg], "--numa", &numa)) {
      iarg++;
    } else if (check_int_arg(arg[iarg], "--kokkos-device-id", &device) ||
               check_int_arg(arg[iarg], "--kokkos-device", &device)) {
      if (check_arg(arg[iarg], "--kokkos-device")) {
        warn_deprecated_command_line_argument("--kokkos-device",
                                              "--kokkos-device-id");
      }
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      kokkos_device_found = true;
      narg--;
    } else if (!kokkos_device_found &&
               (check_int_arg(arg[iarg], "--device-id", &device) ||
                check_int_arg(arg[iarg], "--device", &device))) {
      if (check_arg(arg[iarg], "--device")) {
        warn_deprecated_command_line_argument("--device", "--device-id");
      }
      iarg++;
    } else if (check_arg(arg[iarg], "--kokkos-num-devices") ||
               check_arg(arg[iarg], "--num-devices") ||
               check_arg(arg[iarg], "--kokkos-ndevices") ||
               check_arg(arg[iarg], "--ndevices")) {
      if (check_arg(arg[iarg], "--ndevices")) {
        warn_deprecated_command_line_argument("--ndevices", "--num-devices");
      }
      if (check_arg(arg[iarg], "--kokkos-ndevices")) {
        warn_deprecated_command_line_argument("--kokkos-ndevices",
                                              "--kokkos-num-devices");
      }
      // Find the number of device (expecting --device=XX)
      if (!((strncmp(arg[iarg], "--kokkos-num-devices=", 21) == 0) ||
            (strncmp(arg[iarg], "--num-ndevices=", 14) == 0) ||
            (strncmp(arg[iarg], "--kokkos-ndevices=", 18) == 0) ||
            (strncmp(arg[iarg], "--ndevices=", 11) == 0)))
        throw_runtime_exception(
            "Error: expecting an '=INT[,INT]' after command line argument "
            "'--num-devices/--kokkos-num-devices'. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");

      char* num1      = strchr(arg[iarg], '=') + 1;
      char* num2      = strpbrk(num1, ",");
      int num1_len    = num2 == nullptr ? strlen(num1) : num2 - num1;
      char* num1_only = new char[num1_len + 1];
      strncpy(num1_only, num1, num1_len);
      num1_only[num1_len] = '\0';

      if (!is_unsigned_int(num1_only) || (strlen(num1_only) == 0)) {
        throw_runtime_exception(
            "Error: expecting an integer number after command line argument "
            "'--kokkos-numdevices'. Raised by "
            "Kokkos::initialize(int narg, char* argc[]).");
      }
      if (check_arg(arg[iarg], "--kokkos-num-devices") ||
          check_arg(arg[iarg], "--kokkos-ndevices") || !kokkos_ndevices_found)
        ndevices = std::stoi(num1_only);
      delete[] num1_only;

      if (num2 != nullptr) {
        if ((!is_unsigned_int(num2 + 1)) || (strlen(num2) == 1))
          throw_runtime_exception(
              "Error: expecting an integer number after command line argument "
              "'--kokkos-num-devices=XX,'. Raised by "
              "Kokkos::initialize(int narg, char* argc[]).");

        if (check_arg(arg[iarg], "--kokkos-num-devices") ||
            check_arg(arg[iarg], "--kokkos-ndevices") || !kokkos_ndevices_found)
          skip_device = std::stoi(num2 + 1);
      }

      // Remove the --kokkos-num-devices argument from the list but leave
      // --num-devices
      if (check_arg(arg[iarg], "--kokkos-num-devices") ||
          check_arg(arg[iarg], "--kokkos-ndevices")) {
        for (int k = iarg; k < narg - 1; k++) {
          arg[k] = arg[k + 1];
        }
        kokkos_ndevices_found = true;
        narg--;
      } else {
        iarg++;
      }
    } else if (check_arg(arg[iarg], "--kokkos-disable-warnings")) {
      disable_warnings = true;
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
    } else if (check_arg(arg[iarg], "--kokkos-tune-internals")) {
      tune_internals = true;
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
    } else if (check_str_arg(arg[iarg], "--kokkos-tools-library", tool_lib)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
    } else if (check_str_arg(arg[iarg], "--kokkos-tools-args", tool_args)) {
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
      // strip any leading and/or trailing quotes if they were retained in the
      // string because this will very likely cause parsing issues for tools.
      // If the quotes are retained (via bypassing the shell):
      //    <EXE> --kokkos-tools-args="-c my example"
      // would be tokenized as:
      //    "<EXE>" "\"-c" "my" "example\""
      // instead of:
      //    "<EXE>" "-c" "my" "example"
      if (!tool_args.empty()) {
        if (tool_args.front() == '"') tool_args = tool_args.substr(1);
        if (tool_args.back() == '"')
          tool_args = tool_args.substr(0, tool_args.length() - 1);
      }
      // add the name of the executable to the beginning
      if (narg > 0) tool_args = std::string(arg[0]) + " " + tool_args;
    } else if (check_arg(arg[iarg], "--kokkos-tools-help")) {
      tool_help = true;
      for (int k = iarg; k < narg - 1; k++) {
        arg[k] = arg[k + 1];
      }
      narg--;
    } else if (check_arg(arg[iarg], "--kokkos-help") ||
               check_arg(arg[iarg], "--help")) {
      auto const help_message = R"(
      --------------------------------------------------------------------------------
      -------------Kokkos command line arguments--------------------------------------
      --------------------------------------------------------------------------------
      The following arguments exist also without prefix 'kokkos' (e.g. --help).
      The prefixed arguments will be removed from the list by Kokkos::initialize(),
      the non-prefixed ones are not removed. Prefixed versions take precedence over
      non prefixed ones, and the last occurrence of an argument overwrites prior
      settings.

      --kokkos-help                  : print this message
      --kokkos-disable-warnings      : disable kokkos warning messages
      --kokkos-tune-internals        : allow Kokkos to autotune policies and declare
                                       tuning features through the tuning system. If
                                       left off, Kokkos uses heuristics
      --kokkos-threads=INT           : specify total number of threads or
                                       number of threads per NUMA region if
                                       used in conjunction with '--numa' option.
      --kokkos-numa=INT              : specify number of NUMA regions used by process.
      --kokkos-device-id=INT         : specify device id to be used by Kokkos.
      --kokkos-num-devices=INT[,INT] : used when running MPI jobs. Specify number of
                                       devices per node to be used. Process to device
                                       mapping happens by obtaining the local MPI rank
                                       and assigning devices round-robin. The optional
                                       second argument allows for an existing device
                                       to be ignored. This is most useful on workstations
                                       with multiple GPUs of which one is used to drive
                                       screen output.
      --kokkos-tools-library         : Equivalent to KOKKOS_PROFILE_LIBRARY environment
                                       variable. Must either be full path to library or
                                       name of library if the path is present in the
                                       runtime library search path (e.g. LD_LIBRARY_PATH)
      --kokkos-tools-help            : Query the (loaded) kokkos-tool for its command-line
                                       option support (which should then be passed via
                                       --kokkos-tools-args="...")
      --kokkos-tools-args=STR        : A single (quoted) string of options which will be
                                       whitespace delimited and passed to the loaded
                                       kokkos-tool as command-line arguments. E.g.
                                       `<EXE> --kokkos-tools-args="-c input.txt"` will
                                       pass `<EXE> -c input.txt` as argc/argv to tool
      --------------------------------------------------------------------------------
)";
      std::cout << help_message << std::endl;

      // Remove the --kokkos-help argument from the list but leave --help
      if (check_arg(arg[iarg], "--kokkos-help")) {
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
  if (tool_args.empty() && narg > 0) tool_args = arg[0];
}

void parse_environment_variables(InitArguments& arguments) {
  auto& num_threads      = arguments.num_threads;
  auto& numa             = arguments.num_numa;
  auto& device           = arguments.device_id;
  auto& ndevices         = arguments.ndevices;
  auto& skip_device      = arguments.skip_device;
  auto& disable_warnings = arguments.disable_warnings;
  auto& tune_internals   = arguments.tune_internals;
  auto& tool_lib         = arguments.tool_lib;
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

      std::srand(get_process_id());
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
    const auto _rc = std::regex_constants::icase | std::regex_constants::egrep;
    const auto _re = std::regex("^(true|on|yes|[1-9])$", _rc);
    if (std::regex_match(env_str, _re))
      disable_warnings = true;
    else if (disable_warnings)
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-disable-warnings and "
          "KOKKOS_DISABLE_WARNINGS if both are set. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
  }
  char* env_tuneinternals_str = std::getenv("KOKKOS_TUNE_INTERNALS");
  if (env_tuneinternals_str != nullptr) {
    std::string env_str(env_tuneinternals_str);  // deep-copies string
    for (char& c : env_str) {
      c = toupper(c);
    }
    if ((env_str == "TRUE") || (env_str == "ON") || (env_str == "1"))
      tune_internals = true;
    else if (tune_internals)
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-tune-internals and "
          "KOKKOS_TUNE_INTERNALS if both are set. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
  }
  auto env_tool_lib = std::getenv("KOKKOS_PROFILE_LIBRARY");
  if (env_tool_lib != nullptr) {
    if (!tool_lib.empty() && std::string(env_tool_lib) != tool_lib)
      Impl::throw_runtime_exception(
          "Error: expecting a match between --kokkos-tools-library and "
          "KOKKOS_PROFILE_LIBRARY if both are set. Raised by "
          "Kokkos::initialize(int narg, char* argc[]).");
    else
      tool_lib = env_tool_lib;
  }
}

}  // namespace

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

void initialize(int& narg, char* arg[]) {
  InitArguments arguments;
  Impl::parse_command_line_arguments(narg, arg, arguments);
  Impl::parse_environment_variables(arguments);
  Impl::initialize_internal(arguments);
}

void initialize(InitArguments arguments) {
  Impl::parse_environment_variables(arguments);
  Impl::initialize_internal(arguments);
}

namespace Impl {

void pre_initialize(const InitArguments& args) {
  pre_initialize_internal(args);
}

void post_initialize(const InitArguments& args) {
  post_initialize_internal(args);
}

}  // namespace Impl

void push_finalize_hook(std::function<void()> f) { finalize_hooks.push(f); }

void finalize() { Impl::finalize_internal(); }

void finalize_all() {
  enum : bool { all_spaces = true };
  Impl::finalize_internal(all_spaces);
}

void fence() { Impl::fence_internal(); }

void print_helper(std::ostringstream& out,
                  const std::map<std::string, std::string>& print_me) {
  for (const auto& kv : print_me) {
    out << kv.first << ": " << kv.second << '\n';
  }
}

void print_configuration(std::ostream& out, const bool detail) {
  std::ostringstream msg;

  print_helper(msg, Kokkos::Impl::metadata_map["version_info"]);

  msg << "Compiler:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["compiler_version"]);

  msg << "Architecture:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["architecture"]);

  msg << "Atomics:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["atomics"]);

  msg << "Vectorization:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["vectorization"]);

  msg << "Memory:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["memory"]);

  msg << "Options:" << std::endl;
  print_helper(msg, Kokkos::Impl::metadata_map["options"]);

  Impl::ExecSpaceManager::get_instance().print_configuration(msg, detail);

  out << msg.str() << std::endl;
}

bool is_initialized() noexcept { return g_is_initialized; }

bool show_warnings() noexcept { return g_show_warnings; }
bool tune_internals() noexcept { return g_tune_internals; }

#ifdef KOKKOS_COMPILER_PGI
namespace Impl {
// Bizzarely, an extra jump instruction forces the PGI compiler to not have a
// bug related to (probably?) empty base optimization and/or aggregate
// construction.
void _kokkos_pgi_compiler_bug_workaround() {}
}  // end namespace Impl
#endif

}  // namespace Kokkos
