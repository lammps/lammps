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

#ifndef KOKKOS_CORE_HPP
#define KOKKOS_CORE_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

//----------------------------------------------------------------------------
// Include the execution space header files for the enabled execution spaces.

#include <Kokkos_Core_fwd.hpp>

#include <KokkosCore_Config_DeclareBackend.hpp>

#include <Kokkos_Half.hpp>
#include <Kokkos_AnonymousSpace.hpp>
#include <Kokkos_LogicalSpaces.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_MinMaxClamp.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_MathematicalSpecialFunctions.hpp>
#include <Kokkos_MemoryPool.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_Tuners.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_CopyViews.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <functional>
#include <iosfwd>
#include <memory>
#include <vector>

//----------------------------------------------------------------------------

namespace Kokkos {

void initialize(int& argc, char* argv[]);

void initialize(
    InitializationSettings const& settings = InitializationSettings());

namespace Impl {

void pre_initialize(const InitializationSettings& settings);

void post_initialize(const InitializationSettings& settings);

void declare_configuration_metadata(const std::string& category,
                                    const std::string& key,
                                    const std::string& value);

}  // namespace Impl

KOKKOS_ATTRIBUTE_NODISCARD bool is_initialized() noexcept;
KOKKOS_ATTRIBUTE_NODISCARD bool is_finalized() noexcept;

bool show_warnings() noexcept;
bool tune_internals() noexcept;

/** \brief  Finalize the spaces that were initialized via Kokkos::initialize */
void finalize();

/**
 * \brief Push a user-defined function to be called in
 *   Kokkos::finalize, before any Kokkos state is finalized.
 *
 * \warning Only call this after Kokkos::initialize, but before
 *   Kokkos::finalize.
 *
 * This function is the Kokkos analog to std::atexit.  If you call
 * this with a function f, then your function will get called when
 * Kokkos::finalize is called.  Specifically, it will be called BEFORE
 * Kokkos does any finalization.  This means that all execution
 * spaces, memory spaces, etc. that were initialized will still be
 * initialized when your function is called.
 *
 * Just like std::atexit, if you call push_finalize_hook in sequence
 * with multiple functions (f, g, h), Kokkos::finalize will call them
 * in reverse order (h, g, f), as if popping a stack.  Furthermore,
 * just like std::atexit, if any of your functions throws but does not
 * catch an exception, Kokkos::finalize will call std::terminate.
 */
void push_finalize_hook(std::function<void()> f);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
/** \brief  Finalize all known execution spaces */
KOKKOS_DEPRECATED void finalize_all();
#endif

/** \brief Print "Bill of Materials" */
void print_configuration(std::ostream& os, bool verbose = false);

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const std::string& arg_alloc_label,
                           const size_t arg_alloc_size) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), arg_alloc_label, arg_alloc_size);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const size_t arg_alloc_size) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), "no-label", arg_alloc_size);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void kokkos_free(void* arg_alloc) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::deallocate_tracked(
      arg_alloc);
}

template <class Space = Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_realloc(void* arg_alloc, const size_t arg_alloc_size) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::reallocate_tracked(
      arg_alloc, arg_alloc_size);
}

}  // namespace Kokkos

namespace Kokkos {

/** \brief  ScopeGuard
 *  Some user scope issues have been identified with some Kokkos::finalize
 * calls; ScopeGuard aims to correct these issues.
 *
 *  Two requirements for ScopeGuard:
 *     if Kokkos::is_initialized() in the constructor, don't call
 * Kokkos::initialize or Kokkos::finalize it is not copyable or assignable
 */
namespace Impl {

inline std::string scopeguard_correct_usage() {
  return std::string(
      "Do instead:\n"
      "  std::unique_ptr<Kokkos::ScopeGuard> guard =\n"
      "    !Kokkos::is_initialized() && !Kokkos::is_finalized()?\n"
      "    new ScopeGuard(argc,argv) : nullptr;\n");
}

inline std::string scopeguard_create_while_initialized_warning() {
  return std::string(
             "Kokkos Error: Creating a ScopeGuard while Kokkos is initialized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

inline std::string scopeguard_create_after_finalize_warning() {
  return std::string(
             "Kokkos Error: Creating a ScopeGuard after Kokkos was finalized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

inline std::string scopeguard_destruct_after_finalize_warning() {
  return std::string(
             "Kokkos Error: Destroying a ScopeGuard after Kokkos was finalized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

}  // namespace Impl

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
class KOKKOS_ATTRIBUTE_NODISCARD ScopeGuard {
 public:
#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  KOKKOS_ATTRIBUTE_NODISCARD
#endif
  ScopeGuard(int& argc, char* argv[]) {
    sg_init = false;
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (is_initialized()) {
      std::cerr << Impl::scopeguard_create_while_initialized_warning()
                << std::endl;
    }
    if (is_finalized()) {
      std::cerr << Impl::scopeguard_create_after_finalize_warning()
                << std::endl;
    }
#endif
    if (!is_initialized()) {
      initialize(argc, argv);
      sg_init = true;
    }
  }

#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  KOKKOS_ATTRIBUTE_NODISCARD
#endif
  explicit ScopeGuard(
      const InitializationSettings& settings = InitializationSettings()) {
    sg_init = false;
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (is_initialized()) {
      std::cerr << Impl::scopeguard_create_while_initialized_warning()
                << std::endl;
    }
    if (is_finalized()) {
      std::cerr << Impl::scopeguard_create_after_finalize_warning()
                << std::endl;
    }
#endif
    if (!is_initialized()) {
      initialize(settings);
      sg_init = true;
    }
  }

  ~ScopeGuard() {
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (is_finalized()) {
      std::cerr << Impl::scopeguard_destruct_after_finalize_warning()
                << std::endl;
    }
#endif
    if (is_initialized() && sg_init) {
      finalize();
    }
  }

 private:
  bool sg_init;

 public:
  ScopeGuard& operator=(const ScopeGuard&) = delete;
  ScopeGuard& operator=(ScopeGuard&&) = delete;
  ScopeGuard(const ScopeGuard&)       = delete;
  ScopeGuard(ScopeGuard&&)            = delete;
};

#else  // ifndef KOKKOS_ENABLE_DEPRECATED_CODE3

class KOKKOS_ATTRIBUTE_NODISCARD ScopeGuard {
 public:
#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  KOKKOS_ATTRIBUTE_NODISCARD
#endif
  ScopeGuard(int& argc, char* argv[]) {
    if (is_initialized()) {
      Kokkos::abort(
          Impl::scopeguard_create_while_initialized_warning().c_str());
    }
    if (is_finalized()) {
      Kokkos::abort(Impl::scopeguard_create_after_finalize_warning().c_str());
    }
    initialize(argc, argv);
  }

#if defined(__has_cpp_attribute) && __has_cpp_attribute(nodiscard) >= 201907
  KOKKOS_ATTRIBUTE_NODISCARD
#endif
  ScopeGuard(
      const InitializationSettings& settings = InitializationSettings()) {
    if (is_initialized()) {
      Kokkos::abort(
          Impl::scopeguard_create_while_initialized_warning().c_str());
    }
    if (is_finalized()) {
      Kokkos::abort(Impl::scopeguard_create_after_finalize_warning().c_str());
    }
    initialize(settings);
  }

  ~ScopeGuard() {
    if (is_finalized()) {
      Kokkos::abort(Impl::scopeguard_destruct_after_finalize_warning().c_str());
    }
    finalize();
  }

  ScopeGuard& operator=(const ScopeGuard&) = delete;
  ScopeGuard& operator=(ScopeGuard&&) = delete;
  ScopeGuard(const ScopeGuard&)       = delete;
  ScopeGuard(ScopeGuard&&)            = delete;
};
#endif

}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
// Partitioning an Execution Space: expects space and integer arguments for
// relative weight
//   Customization point for backends
//   Default behavior is to return the passed in instance
template <class ExecSpace, class... Args>
std::vector<ExecSpace> partition_space(ExecSpace space, Args...) {
  static_assert(is_execution_space<ExecSpace>::value,
                "Kokkos Error: partition_space expects an Execution Space as "
                "first argument");
#ifdef __cpp_fold_expressions
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");
#endif
  std::vector<ExecSpace> instances(sizeof...(Args));
  for (int s = 0; s < int(sizeof...(Args)); s++) instances[s] = space;
  return instances;
}

template <class ExecSpace, class T>
std::vector<ExecSpace> partition_space(ExecSpace space,
                                       std::vector<T>& weights) {
  static_assert(is_execution_space<ExecSpace>::value,
                "Kokkos Error: partition_space expects an Execution Space as "
                "first argument");
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  std::vector<ExecSpace> instances(weights.size());
  for (int s = 0; s < int(weights.size()); s++) instances[s] = space;
  return instances;
}
}  // namespace Experimental
}  // namespace Kokkos

#include <Kokkos_Crs.hpp>
#include <Kokkos_WorkGraphPolicy.hpp>
// Including this in Kokkos_Parallel_Reduce.hpp led to a circular dependency
// because Kokkos::Sum is used in Kokkos_Combined_Reducer.hpp and the default.
// The real answer is to finally break up Kokkos_Parallel_Reduce.hpp into
// smaller parts...
#include <impl/Kokkos_Combined_Reducer.hpp>
// Yet another workaround to deal with circular dependency issues because the
// implementation of the RAII wrapper is using Kokkos::single.
#include <Kokkos_AcquireUniqueTokenImpl.hpp>

// Specializations required after core definitions
#include <KokkosCore_Config_PostInclude.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#endif
