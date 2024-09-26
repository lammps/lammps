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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Command_Line_Parsing.hpp>
#include <impl/Kokkos_ParseCommandLineArgumentsAndEnvironmentVariables.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stack>
#include <functional>
#include <list>
#include <cerrno>
#include <random>
#include <regex>
#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#endif

//----------------------------------------------------------------------------
namespace {
bool g_is_initialized = false;
bool g_is_finalized   = false;
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

void combine(Kokkos::InitializationSettings& out,
             Kokkos::InitializationSettings const& in) {
#define KOKKOS_IMPL_COMBINE_SETTING(NAME) \
  if (in.has_##NAME()) {                  \
    out.set_##NAME(in.get_##NAME());      \
  }                                       \
  static_assert(true, "no-op to require trailing semicolon")
  KOKKOS_IMPL_COMBINE_SETTING(num_threads);
  KOKKOS_IMPL_COMBINE_SETTING(map_device_id_by);
  KOKKOS_IMPL_COMBINE_SETTING(device_id);
  KOKKOS_IMPL_COMBINE_SETTING(disable_warnings);
  KOKKOS_IMPL_COMBINE_SETTING(print_configuration);
  KOKKOS_IMPL_COMBINE_SETTING(tune_internals);
  KOKKOS_IMPL_COMBINE_SETTING(tools_help);
  KOKKOS_IMPL_COMBINE_SETTING(tools_libs);
  KOKKOS_IMPL_COMBINE_SETTING(tools_args);
#undef KOKKOS_IMPL_COMBINE_SETTING
}

void combine(Kokkos::InitializationSettings& out,
             Kokkos::Tools::InitArguments const& in) {
  using Kokkos::Tools::InitArguments;
  if (in.help != InitArguments::PossiblyUnsetOption::unset) {
    out.set_tools_help(in.help == InitArguments::PossiblyUnsetOption::on);
  }
  if (in.lib != InitArguments::unset_string_option) {
    out.set_tools_libs(in.lib);
  }
  if (in.args != InitArguments::unset_string_option) {
    out.set_tools_args(in.args);
  }
}

void combine(Kokkos::Tools::InitArguments& out,
             Kokkos::InitializationSettings const& in) {
  using Kokkos::Tools::InitArguments;
  if (in.has_tools_help()) {
    out.help = in.get_tools_help() ? InitArguments::PossiblyUnsetOption::on
                                   : InitArguments::PossiblyUnsetOption::off;
  }
  if (in.has_tools_libs()) {
    out.lib = in.get_tools_libs();
  }
  if (in.has_tools_args()) {
    out.args = in.get_tools_args();
  }
}

int get_device_count() {
#if defined(KOKKOS_ENABLE_CUDA)
  int count;
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&count));
  return count;
#elif defined(KOKKOS_ENABLE_HIP)
  int count;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&count));
  return count;
#elif defined(KOKKOS_ENABLE_SYCL)
  return Kokkos::Experimental::Impl::get_sycl_devices().size();
#elif defined(KOKKOS_ENABLE_OPENACC)
  return acc_get_num_devices(
      Kokkos::Experimental::Impl::OpenACC_Traits::dev_type);
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  return omp_get_num_devices();
#else
  Kokkos::abort("implementation bug");
  return -1;
#endif
}

unsigned get_process_id() {
#ifdef _WIN32
  return unsigned(GetCurrentProcessId());
#else
  return unsigned(getpid());
#endif
}

bool is_valid_num_threads(int x) { return x > 0; }

bool is_valid_device_id(int x) { return x >= 0; }

bool is_valid_map_device_id_by(std::string const& x) {
  return x == "mpi_rank" || x == "random";
}

}  // namespace

std::vector<int> const& Kokkos::Impl::get_visible_devices() {
  static auto devices = get_visible_devices(get_device_count());
  return devices;
}

[[nodiscard]] int Kokkos::device_id() noexcept {
#if defined(KOKKOS_ENABLE_CUDA)
  int device = Cuda().cuda_device();
#elif defined(KOKKOS_ENABLE_HIP)
  int device = HIP().hip_device();
#elif defined(KOKKOS_ENABLE_OPENACC)
  int device = Experimental::OpenACC().acc_device_number();
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  int device = omp_get_default_device();  // FIXME_OPENMPTARGET
#elif defined(KOKKOS_ENABLE_SYCL)
  int device = Experimental::Impl::SYCLInternal::m_syclDev;
#else
  int device = -1;
  return device;
#endif
  auto const& visible_devices = Impl::get_visible_devices();
  for (std::size_t i = 0; i < visible_devices.size(); ++i) {
    if (visible_devices[i] == device) {
      return i;
    }
  }
  Kokkos::abort("Unexpected error: cannot determine device id");
  return -1;
}

[[nodiscard]] int Kokkos::num_devices() noexcept {
  if constexpr (std::is_same_v<DefaultExecutionSpace,
                               DefaultHostExecutionSpace>) {
    return -1;  // no GPU backend enabled
  } else {
    return Impl::get_visible_devices().size();
  }
}

[[nodiscard]] int Kokkos::num_threads() noexcept {
  return DefaultHostExecutionSpace().concurrency();
}

Kokkos::Impl::ExecSpaceManager& Kokkos::Impl::ExecSpaceManager::get_instance() {
  static ExecSpaceManager space_initializer = {};
  return space_initializer;
}

void Kokkos::Impl::ExecSpaceManager::register_space_factory(
    const std::string name, std::unique_ptr<ExecSpaceBase> space) {
  exec_space_factory_list[name] = std::move(space);
}

void Kokkos::Impl::ExecSpaceManager::initialize_spaces(
    const InitializationSettings& settings) {
  // Note: the names of the execution spaces, used as keys in the map, encode
  // the ordering of the initialization code from the old initialization stuff.
  // Eventually, we may want to do something less brittle than this, but for now
  // we're just preserving compatibility with the old implementation.
  for (auto& to_init : exec_space_factory_list) {
    to_init.second->initialize(settings);
  }
}

void Kokkos::Impl::ExecSpaceManager::finalize_spaces() {
  for (auto& to_finalize : exec_space_factory_list) {
    to_finalize.second->finalize();
  }
}

void Kokkos::Impl::ExecSpaceManager::static_fence(const std::string& name) {
  for (auto& to_fence : exec_space_factory_list) {
    to_fence.second->static_fence(name);
  }
}
void Kokkos::Impl::ExecSpaceManager::print_configuration(std::ostream& os,
                                                         bool verbose) {
  for (auto const& to_print : exec_space_factory_list) {
    to_print.second->print_configuration(os, verbose);
  }
}

int Kokkos::Impl::get_ctest_gpu(int local_rank) {
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
  assert(local_rank >= 0);
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

std::vector<int> Kokkos::Impl::get_visible_devices(int device_count) {
  std::vector<int> visible_devices;
  char* env_visible_devices = std::getenv("KOKKOS_VISIBLE_DEVICES");
  if (env_visible_devices) {
    std::stringstream ss(env_visible_devices);
    for (int i; ss >> i;) {
      visible_devices.push_back(i);
      if (ss.peek() == ',') ss.ignore();
    }
    for (auto id : visible_devices) {
      if (id < 0) {
        ss << "Error: Invalid device id '" << id
           << "' in environment variable 'KOKKOS_VISIBLE_DEVICES="
           << env_visible_devices << "'."
           << " Device id cannot be negative!"
           << " Raised by Kokkos::initialize().\n";
      }
      if (id >= device_count) {
        ss << "Error: Invalid device id '" << id
           << "' in environment variable 'KOKKOS_VISIBLE_DEVICES="
           << env_visible_devices << "'."
           << " Device id must be smaller than the number of GPUs available"
           << " for execution '" << device_count << "'!"
           << " Raised by Kokkos::initialize().\n";
      }
    }
  } else {
    for (int i = 0; i < device_count; ++i) {
      visible_devices.push_back(i);
    }
  }
  if (visible_devices.empty()) {
    Kokkos::abort(
        "Error: no GPU available for execution.\n"
        " Raised by Kokkos::initialize().\n");
  }
  return visible_devices;
}

std::optional<int> Kokkos::Impl::get_gpu(
    const InitializationSettings& settings) {
  std::vector<int> visible_devices = get_visible_devices(get_device_count());
  int const num_devices            = visible_devices.size();
  // device_id is provided
  if (settings.has_device_id()) {
    int const id = settings.get_device_id();
    if (id < 0) {
      std::stringstream ss;
      ss << "Error: Requested GPU with invalid id '" << id << "'."
         << " Device id cannot be negative!"
         << " Raised by Kokkos::initialize().\n";
      Kokkos::abort(ss.str().c_str());
    }
    if (id >= num_devices) {
      std::stringstream ss;
      ss << "Error: Requested GPU with id '" << id << "' but only "
         << num_devices << "GPU(s) available!"
         << " Raised by Kokkos::initialize().\n";
      Kokkos::abort(ss.str().c_str());
    }
    return visible_devices[settings.get_device_id()];
  }

  // either random or round-robin assignment based on local MPI rank
  if (settings.has_map_device_id_by() &&
      !is_valid_map_device_id_by(settings.get_map_device_id_by())) {
    std::stringstream ss;
    ss << "Error: map_device_id_by setting '" << settings.get_map_device_id_by()
       << "' is not recognized."
       << " Raised by Kokkos::initialize().\n";
    Kokkos::abort(ss.str().c_str());
  }

  if (settings.has_map_device_id_by() &&
      settings.get_map_device_id_by() == "random") {
    std::default_random_engine gen(get_process_id());
    std::uniform_int_distribution<int> distribution(0, num_devices - 1);
    return visible_devices[distribution(gen)];
  }

  // either map_device_id_by is not specified or it is mpi_rank
  if (settings.has_map_device_id_by() &&
      settings.get_map_device_id_by() != "mpi_rank") {
    Kokkos::abort("implementation bug");
  }

  int const mpi_local_rank = mpi_local_rank_on_node();

  // if unable to detect local MPI rank return nullopt to delegate device
  // selection to the backend
  if (mpi_local_rank < 0) {
    if (settings.has_map_device_id_by()) {
      std::cerr << "Warning: unable to detect local MPI rank."
                << " Falling back to the first GPU available for execution."
                << " Raised by Kokkos::initialize()." << std::endl;
    }
    return std::nullopt;
  }

  // use device assigned by CTest when resource allocation is activated
  if (std::getenv("CTEST_KOKKOS_DEVICE_TYPE") &&
      std::getenv("CTEST_RESOURCE_GROUP_COUNT")) {
    return get_ctest_gpu(mpi_local_rank);
  }

  return visible_devices[mpi_local_rank % visible_devices.size()];
}

namespace {

void initialize_backends(const Kokkos::InitializationSettings& settings) {
  Kokkos::Impl::ExecSpaceManager::get_instance().initialize_spaces(settings);
}

void initialize_profiling(const Kokkos::Tools::InitArguments& args) {
  auto initialization_status =
      Kokkos::Tools::Impl::initialize_tools_subsystem(args);
  if (initialization_status.result ==
      Kokkos::Tools::Impl::InitializationStatus::InitializationResult::
          help_request) {
    g_is_initialized = true;
    ::Kokkos::finalize();
    std::exit(EXIT_SUCCESS);
  } else if (initialization_status.result ==
             Kokkos::Tools::Impl::InitializationStatus::InitializationResult::
                 success) {
    Kokkos::Tools::parseArgs(args.args);
    for (const auto& category_value : metadata_map) {
      for (const auto& key_value : category_value.second) {
        Kokkos::Tools::declareMetadata(key_value.first, key_value.second);
      }
    }
  } else {
    std::cerr << "Error initializing Kokkos Tools subsystem" << std::endl;
    g_is_initialized = true;
    ::Kokkos::finalize();
    std::exit(EXIT_FAILURE);
  }
}

std::string version_string_from_int(int version_number) {
  std::stringstream str_builder;
  str_builder << version_number / 10000 << "." << (version_number % 10000) / 100
              << "." << version_number % 100;
  return str_builder.str();
}

void pre_initialize_internal(const Kokkos::InitializationSettings& settings) {
  if (settings.has_disable_warnings() && settings.get_disable_warnings())
    g_show_warnings = false;
  if (settings.has_tune_internals() && settings.get_tune_internals())
    g_tune_internals = true;
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
#ifdef KOKKOS_COMPILER_INTEL
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_INTEL",
                                 std::to_string(KOKKOS_COMPILER_INTEL));
  declare_configuration_metadata("tools_only", "compiler_family", "intel");
#endif
#ifdef KOKKOS_COMPILER_INTEL_LLVM
  declare_configuration_metadata("compiler_version",
                                 "KOKKOS_COMPILER_INTEL_LLVM",
                                 std::to_string(KOKKOS_COMPILER_INTEL_LLVM));
  declare_configuration_metadata("tools_only", "compiler_family", "intel_llvm");
#endif
#ifdef KOKKOS_COMPILER_NVCC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_NVCC",
                                 std::to_string(KOKKOS_COMPILER_NVCC));
  declare_configuration_metadata("tools_only", "compiler_family", "nvcc");
#endif
#ifdef KOKKOS_COMPILER_NVHPC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_NVHPC",
                                 std::to_string(KOKKOS_COMPILER_NVHPC));
  declare_configuration_metadata("tools_only", "compiler_family", "pgi");
#endif
#ifdef KOKKOS_COMPILER_MSVC
  declare_configuration_metadata("compiler_version", "KOKKOS_COMPILER_MSVC",
                                 std::to_string(KOKKOS_COMPILER_MSVC));
  declare_configuration_metadata("tools_only", "compiler_family", "msvc");
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

#ifdef KOKKOS_ENABLE_ASM
  declare_configuration_metadata("options", "KOKKOS_ENABLE_ASM", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_ASM", "no");
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
#ifdef KOKKOS_ENABLE_CXX23
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX23", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX23", "no");
#endif
#ifdef KOKKOS_ENABLE_CXX26
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX26", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_CXX26", "no");
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
#ifdef KOKKOS_ENABLE_LIBDL
  declare_configuration_metadata("options", "KOKKOS_ENABLE_LIBDL", "yes");
#else
  declare_configuration_metadata("options", "KOKKOS_ENABLE_LIBDL", "no");
#endif

  declare_configuration_metadata("architecture", "Default Device",
                                 typeid(Kokkos::DefaultExecutionSpace).name());

#if defined(KOKKOS_ARCH_A64FX)
  declare_configuration_metadata("architecture", "CPU architecture", "A64FX");
#elif defined(KOKKOS_ARCH_AMDAVX)
  declare_configuration_metadata("architecture", "CPU architecture", "AMDAVX");
#elif defined(KOKKOS_ARCH_ARMV80)
  declare_configuration_metadata("architecture", "CPU architecture", "ARMV80");
#elif defined(KOKKOS_ARCH_ARMV81)
  declare_configuration_metadata("architecture", "CPU architecture", "ARMV81");
#elif defined(KOKKOS_ARCH_ARMV8_THUNDERX)
  declare_configuration_metadata("architecture", "CPU architecture",
                                 "ARMV8_THUNDERX");
#elif defined(KOKKOS_ARCH_ARMV8_THUNDERX2)
  declare_configuration_metadata("architecture", "CPU architecture",
                                 "ARMV8_THUNDERX2");
#elif defined(KOKKOS_ARCH_BDW)
  declare_configuration_metadata("architecture", "CPU architecture", "BDW");
#elif defined(KOKKOS_ARCH_HSW)
  declare_configuration_metadata("architecture", "CPU architecture", "HSW");
#elif defined(KOKKOS_ARCH_ICL)
  declare_configuration_metadata("architecture", "CPU architecture", "ICL");
#elif defined(KOKKOS_ARCH_ICX)
  declare_configuration_metadata("architecture", "CPU architecture", "ICX");
#elif defined(KOKKOS_ARCH_KNC)
  declare_configuration_metadata("architecture", "CPU architecture", "KNC");
#elif defined(KOKKOS_ARCH_KNL)
  declare_configuration_metadata("architecture", "CPU architecture", "KNL");
#elif defined(KOKKOS_ARCH_NATIVE)
  declare_configuration_metadata("architecture", "CPU architecture", "NATIVE");
#elif defined(KOKKOS_ARCH_POWER8)
  declare_configuration_metadata("architecture", "CPU architecture", "POWER8");
#elif defined(KOKKOS_ARCH_POWER9)
  declare_configuration_metadata("architecture", "CPU architecture", "POWER9");
#elif defined(KOKKOS_ARCH_SKL)
  declare_configuration_metadata("architecture", "CPU architecture", "SKL");
#elif defined(KOKKOS_ARCH_SKX)
  declare_configuration_metadata("architecture", "CPU architecture", "SKX");
#elif defined(KOKKOS_ARCH_SNB)
  declare_configuration_metadata("architecture", "CPU architecture", "SNB");
#elif defined(KOKKOS_ARCH_SPR)
  declare_configuration_metadata("architecture", "CPU architecture", "SPR");
#elif defined(KOKKOS_ARCH_AMD_ZEN)
  declare_configuration_metadata("architecture", "CPU architecture", "AMD_ZEN");
#elif defined(KOKKOS_ARCH_AMD_ZEN2)
  declare_configuration_metadata("architecture", "CPU architecture",
                                 "AMD_ZEN2");
#elif defined(KOKKOS_ARCH_AMD_ZEN3)
  declare_configuration_metadata("architecture", "CPU architecture",
                                 "AMD_ZEN3");
#elif defined(KOKKOS_ARCH_RISCV_SG2042)
  declare_configuration_metadata("architecture", "CPU architecture",
                                 "SG2042 (RISC-V)")
#else
  declare_configuration_metadata("architecture", "CPU architecture", "none");
#endif

#if defined(KOKKOS_ARCH_INTEL_GEN)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_GEN");
#elif defined(KOKKOS_ARCH_INTEL_DG1)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_DG1");
#elif defined(KOKKOS_ARCH_INTEL_GEN9)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_GEN9");
#elif defined(KOKKOS_ARCH_INTEL_GEN11)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_GEN11");
#elif defined(KOKKOS_ARCH_INTEL_GEN12LP)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_GEN12LP");
#elif defined(KOKKOS_ARCH_INTEL_XEHP)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_XEHP");
#elif defined(KOKKOS_ARCH_INTEL_PVC)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "INTEL_PVC");

#elif defined(KOKKOS_ARCH_KEPLER30)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "KEPLER30");
#elif defined(KOKKOS_ARCH_KEPLER32)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "KEPLER32");
#elif defined(KOKKOS_ARCH_KEPLER35)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "KEPLER35");
#elif defined(KOKKOS_ARCH_KEPLER37)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "KELPER37");
#elif defined(KOKKOS_ARCH_MAXWELL50)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "MAXWELL50");
#elif defined(KOKKOS_ARCH_MAXWELL52)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "MAXWELL52");
#elif defined(KOKKOS_ARCH_MAXWELL53)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "MAXWELL53");
#elif defined(KOKKOS_ARCH_PASCAL60)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "PASCAL60");
#elif defined(KOKKOS_ARCH_PASCAL61)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "PASCAL61");
#elif defined(KOKKOS_ARCH_VOLTA70)
  declare_configuration_metadata("architecture", "GPU architecture", "VOLTA70");
#elif defined(KOKKOS_ARCH_VOLTA72)
  declare_configuration_metadata("architecture", "GPU architecture", "VOLTA72");
#elif defined(KOKKOS_ARCH_TURING75)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "TURING75");
#elif defined(KOKKOS_ARCH_AMPERE80)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMPERE80");
#elif defined(KOKKOS_ARCH_AMPERE86)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMPERE86");
#elif defined(KOKKOS_ARCH_ADA89)
  declare_configuration_metadata("architecture", "GPU architecture", "ADA89");
#elif defined(KOKKOS_ARCH_HOPPER90)
      declare_configuration_metadata("architecture", "GPU architecture",
                                     "HOPPER90");
#elif defined(KOKKOS_ARCH_AMD_GFX906)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMD_GFX906");
#elif defined(KOKKOS_ARCH_AMD_GFX908)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMD_GFX908");
#elif defined(KOKKOS_ARCH_AMD_GFX90A)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMD_GFX90A");
#elif defined(KOKKOS_ARCH_AMD_GFX1030)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMD_GFX1030");
#elif defined(KOKKOS_ARCH_AMD_GFX1100)
  declare_configuration_metadata("architecture", "GPU architecture",
                                 "AMD_GFX1100");

#else
  declare_configuration_metadata("architecture", "GPU architecture", "none");
#endif

#ifdef KOKKOS_IMPL_32BIT
  declare_configuration_metadata("architecture", "platform", "32bit");
#else
  declare_configuration_metadata("architecture", "platform", "64bit");
#endif
}

void post_initialize_internal(const Kokkos::InitializationSettings& settings) {
  Kokkos::Tools::InitArguments tools_init_arguments;
  combine(tools_init_arguments, settings);
  initialize_profiling(tools_init_arguments);
  g_is_initialized = true;
  if (settings.has_print_configuration() &&
      settings.get_print_configuration()) {
    ::Kokkos::print_configuration(std::cout);
  }
}

void initialize_internal(const Kokkos::InitializationSettings& settings) {
  // The tool initialization is only called in post_initialize_internal.
  // Pausing tools here, so that if someone has set callbacks programmatically
  // these callbacks are not called inside the backend initialization, before
  // the tool initialization happened.
  Kokkos::Tools::Experimental::pause_tools();
  pre_initialize_internal(settings);
  initialize_backends(settings);
  Kokkos::Tools::Experimental::resume_tools();
  post_initialize_internal(settings);
}

// declared noexcept such that std::terminate is called if any of the registered
// function throws
void call_registered_finalize_hook_functions() noexcept {
  while (!finalize_hooks.empty()) {
    auto const& func = finalize_hooks.top();
    func();
    finalize_hooks.pop();
  }
}

void pre_finalize_internal() {
  call_registered_finalize_hook_functions();
  Kokkos::Profiling::finalize();
}

void post_finalize_internal() {
  g_is_initialized = false;
  g_is_finalized   = true;
  g_show_warnings  = true;
  g_tune_internals = false;
}

void fence_internal(const std::string& name) {
  Kokkos::Impl::ExecSpaceManager::get_instance().static_fence(name);
}

void print_help_message() {
  auto const help_message = R"(
--------------------------------------------------------------------------------
-------------Kokkos command line arguments--------------------------------------
--------------------------------------------------------------------------------
This program is using Kokkos.  You can use the following command line flags to
control its behavior:

Kokkos Core Options:
  --kokkos-help                  : print this message
  --kokkos-disable-warnings      : disable kokkos warning messages
  --kokkos-print-configuration   : print configuration
  --kokkos-tune-internals        : allow Kokkos to autotune policies and declare
                                   tuning features through the tuning system. If
                                   left off, Kokkos uses heuristics
  --kokkos-num-threads=INT       : specify total number of threads to use for
                                   parallel regions on the host.
  --kokkos-device-id=INT         : specify device id to be used by Kokkos.
  --kokkos-map-device-id-by=(random|mpi_rank)
                                 : strategy to select device-id automatically from
                                   available devices.
                                   - random:   choose a random device from available.
                                   - mpi_rank: choose device-id based on a round robin
                                               assignment of local MPI ranks.
                                               Works with OpenMPI, MVAPICH, SLURM, and
                                               derived implementations.

Kokkos Tools Options:
  --kokkos-tools-libs=STR        : Specify which of the tools to use. Must either
                                   be full path to library or name of library if the
                                   path is present in the runtime library search path
                                   (e.g. LD_LIBRARY_PATH)
  --kokkos-tools-help            : Query the (loaded) kokkos-tool for its command-line
                                   option support (which should then be passed via
                                   --kokkos-tools-args="...")
  --kokkos-tools-args=STR        : A single (quoted) string of options which will be
                                   whitespace delimited and passed to the loaded
                                   kokkos-tool as command-line arguments. E.g.
                                   `<EXE> --kokkos-tools-args="-c input.txt"` will
                                   pass `<EXE> -c input.txt` as argc/argv to tool

Except for --kokkos[-tools]-help, you can alternatively set the corresponding
environment variable of a flag (all letters in upper-case and underscores
instead of hyphens). For example, to disable warning messages, you can either
specify --kokkos-disable-warnings or set the KOKKOS_DISABLE_WARNINGS
environment variable to yes.

Join us on Slack, visit https://kokkosteam.slack.com
Report bugs to https://github.com/kokkos/kokkos/issues
--------------------------------------------------------------------------------
)";
  std::cout << help_message << std::endl;
}

}  // namespace

void Kokkos::Impl::parse_command_line_arguments(
    int& argc, char* argv[], InitializationSettings& settings) {
  Tools::InitArguments tools_init_arguments;
  combine(tools_init_arguments, settings);
  Tools::Impl::parse_command_line_arguments(argc, argv, tools_init_arguments);
  combine(settings, tools_init_arguments);

  int num_threads;
  int device_id;
  std::string map_device_id_by;
  bool disable_warnings;
  bool print_configuration;
  bool tune_internals;

  bool help_flag = false;

  int iarg = 0;
  while (iarg < argc) {
    bool remove_flag = false;

    if (check_arg_int(argv[iarg], "--kokkos-num-threads", num_threads)) {
      if (!is_valid_num_threads(num_threads)) {
        std::stringstream ss;
        ss << "Error: command line argument '" << argv[iarg] << "' is invalid."
           << " The number of threads must be greater than or equal to one."
           << " Raised by Kokkos::initialize().\n";
        Kokkos::abort(ss.str().c_str());
      }
      settings.set_num_threads(num_threads);
      remove_flag = true;
    } else if (check_arg_int(argv[iarg], "--kokkos-device-id", device_id)) {
      if (!is_valid_device_id(device_id)) {
        std::stringstream ss;
        ss << "Error: command line argument '" << argv[iarg] << "' is invalid."
           << " The device id must be greater than or equal to zero."
           << " Raised by Kokkos::initialize().\n";
        Kokkos::abort(ss.str().c_str());
      }
      settings.set_device_id(device_id);
      remove_flag = true;
    } else if (check_arg_bool(argv[iarg], "--kokkos-disable-warnings",
                              disable_warnings)) {
      settings.set_disable_warnings(disable_warnings);
      remove_flag = true;
    } else if (check_arg_bool(argv[iarg], "--kokkos-print-configuration",
                              print_configuration)) {
      settings.set_print_configuration(print_configuration);
      remove_flag = true;
    } else if (check_arg_bool(argv[iarg], "--kokkos-tune-internals",
                              tune_internals)) {
      settings.set_tune_internals(tune_internals);
      remove_flag = true;
    } else if (check_arg(argv[iarg], "--kokkos-help") ||
               check_arg(argv[iarg], "--help")) {
      help_flag   = true;
      remove_flag = std::string(argv[iarg]).find("--kokkos-") == 0;
    } else if (check_arg_str(argv[iarg], "--kokkos-map-device-id-by",
                             map_device_id_by)) {
      if (!is_valid_map_device_id_by(map_device_id_by)) {
        std::stringstream ss;
        ss << "Warning: command line argument '--kokkos-map-device-id-by="
           << map_device_id_by << "' is not recognized."
           << " Raised by Kokkos::initialize().\n";
        Kokkos::abort(ss.str().c_str());
      }
      settings.set_map_device_id_by(map_device_id_by);
      remove_flag = true;
    } else if (std::regex_match(argv[iarg],
                                std::regex("-?-kokkos.*", std::regex::egrep))) {
      warn_not_recognized_command_line_argument(argv[iarg]);
    }

    if (remove_flag) {
      // Shift the remainder of the argv list by one.  Note that argv has
      // (argc + 1) arguments, the last one always being nullptr.  The following
      // loop moves the trailing nullptr element as well
      for (int k = iarg; k < argc; ++k) {
        argv[k] = argv[k + 1];
      }
      argc--;
    } else {
      iarg++;
    }
  }

  if (help_flag) {
    print_help_message();
  }

  if ((tools_init_arguments.args ==
       Kokkos::Tools::InitArguments::unset_string_option) &&
      argc > 0) {
    settings.set_tools_args(argv[0]);
  }
}

void Kokkos::Impl::parse_environment_variables(
    InitializationSettings& settings) {
  Tools::InitArguments tools_init_arguments;
  combine(tools_init_arguments, settings);
  auto init_result =
      Tools::Impl::parse_environment_variables(tools_init_arguments);
  if (init_result.result ==
      Tools::Impl::InitializationStatus::environment_argument_mismatch) {
    Impl::throw_runtime_exception(init_result.error_message);
  }
  combine(settings, tools_init_arguments);

  int num_threads;
  if (check_env_int("KOKKOS_NUM_THREADS", num_threads)) {
    if (!is_valid_num_threads(num_threads)) {
      std::stringstream ss;
      ss << "Error: environment variable 'KOKKOS_NUM_THREADS=" << num_threads
         << "' is invalid."
         << " The number of threads must be greater than or equal to one."
         << " Raised by Kokkos::initialize().\n";
      Kokkos::abort(ss.str().c_str());
    }
    settings.set_num_threads(num_threads);
  }
  int device_id;
  if (check_env_int("KOKKOS_DEVICE_ID", device_id)) {
    if (!is_valid_device_id(device_id)) {
      std::stringstream ss;
      ss << "Error: environment variable 'KOKKOS_DEVICE_ID" << device_id
         << "' is invalid."
         << " The device id must be greater than or equal to zero."
         << " Raised by Kokkos::initialize().\n";
      Kokkos::abort(ss.str().c_str());
    }
    settings.set_device_id(device_id);
  }
  bool disable_warnings;
  if (check_env_bool("KOKKOS_DISABLE_WARNINGS", disable_warnings)) {
    settings.set_disable_warnings(disable_warnings);
  }
  bool print_configuration;
  if (check_env_bool("KOKKOS_PRINT_CONFIGURATION", print_configuration)) {
    settings.set_print_configuration(print_configuration);
  }
  bool tune_internals;
  if (check_env_bool("KOKKOS_TUNE_INTERNALS", tune_internals)) {
    settings.set_tune_internals(tune_internals);
  }
  char const* map_device_id_by = std::getenv("KOKKOS_MAP_DEVICE_ID_BY");
  if (map_device_id_by != nullptr) {
    if (std::getenv("KOKKOS_DEVICE_ID")) {
      std::cerr << "Warning: environment variable KOKKOS_MAP_DEVICE_ID_BY"
                << "ignored since KOKKOS_DEVICE_ID is specified."
                << " Raised by Kokkos::initialize()." << std::endl;
    }
    if (!is_valid_map_device_id_by(map_device_id_by)) {
      std::stringstream ss;
      ss << "Warning: environment variable 'KOKKOS_MAP_DEVICE_ID_BY="
         << map_device_id_by << "' is not recognized."
         << " Raised by Kokkos::initialize().\n";
      Kokkos::abort(ss.str().c_str());
    }
    settings.set_map_device_id_by(map_device_id_by);
  }
}

//----------------------------------------------------------------------------
namespace {
bool kokkos_initialize_was_called() {
  return Kokkos::is_initialized() || Kokkos::is_finalized();
}
bool kokkos_finalize_was_called() { return Kokkos::is_finalized(); }
}  // namespace

void Kokkos::initialize(int& argc, char* argv[]) {
  if (kokkos_initialize_was_called()) {
    Kokkos::abort(
        "Error: Kokkos::initialize() has already been called."
        " Kokkos can be initialized at most once.\n");
  }
  InitializationSettings settings;
  Impl::parse_environment_variables(settings);
  Impl::parse_command_line_arguments(argc, argv, settings);
  initialize_internal(settings);
}

void Kokkos::initialize(InitializationSettings const& settings) {
  if (kokkos_initialize_was_called()) {
    Kokkos::abort(
        "Error: Kokkos::initialize() has already been called."
        " Kokkos can be initialized at most once.\n");
  }
  InitializationSettings tmp;
  Impl::parse_environment_variables(tmp);
  combine(tmp, settings);
  initialize_internal(tmp);
}

void Kokkos::Impl::pre_initialize(const InitializationSettings& settings) {
  pre_initialize_internal(settings);
}

void Kokkos::Impl::post_initialize(const InitializationSettings& settings) {
  post_initialize_internal(settings);
}

void Kokkos::Impl::pre_finalize() { pre_finalize_internal(); }

void Kokkos::Impl::post_finalize() { post_finalize_internal(); }

void Kokkos::push_finalize_hook(std::function<void()> f) {
  finalize_hooks.push(f);
}

void Kokkos::finalize() {
  if (!kokkos_initialize_was_called()) {
    Kokkos::abort(
        "Error: Kokkos::finalize() may only be called after Kokkos has been "
        "initialized.\n");
  }
  if (kokkos_finalize_was_called()) {
    Kokkos::abort("Error: Kokkos::finalize() has already been called.\n");
  }
  pre_finalize_internal();
  Impl::ExecSpaceManager::get_instance().finalize_spaces();
  post_finalize_internal();
}

#ifdef KOKKOS_COMPILER_INTEL
void Kokkos::fence() { fence("Kokkos::fence: Unnamed Global Fence"); }
#endif
void Kokkos::fence(const std::string& name) { fence_internal(name); }

namespace {
void print_helper(std::ostream& os,
                  const std::map<std::string, std::string>& print_me) {
  for (const auto& kv : print_me) {
    os << "  " << kv.first << ": " << kv.second << '\n';
  }
}
}  // namespace

void Kokkos::print_configuration(std::ostream& os, bool verbose) {
  print_helper(os, metadata_map["version_info"]);

  os << "Compiler:\n";
  print_helper(os, metadata_map["compiler_version"]);

  os << "Architecture:\n";
  print_helper(os, metadata_map["architecture"]);

  os << "Atomics:\n";
  print_helper(os, metadata_map["atomics"]);

  os << "Vectorization:\n";
  print_helper(os, metadata_map["vectorization"]);

  os << "Memory:\n";
  print_helper(os, metadata_map["memory"]);

  os << "Options:\n";
  print_helper(os, metadata_map["options"]);

  Impl::ExecSpaceManager::get_instance().print_configuration(os, verbose);
}

[[nodiscard]] bool Kokkos::is_initialized() noexcept {
  return g_is_initialized;
}

[[nodiscard]] bool Kokkos::is_finalized() noexcept { return g_is_finalized; }

bool Kokkos::show_warnings() noexcept { return g_show_warnings; }

bool Kokkos::tune_internals() noexcept { return g_tune_internals; }
