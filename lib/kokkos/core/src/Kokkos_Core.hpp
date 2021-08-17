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

//----------------------------------------------------------------------------
// Include the execution space header files for the enabled execution spaces.

#include <Kokkos_Core_fwd.hpp>

// Fundamental type description for half precision
// Should not rely on other backend infrastructure
#include <Kokkos_Half.hpp>
#include <KokkosCore_Config_DeclareBackend.hpp>

#include <Kokkos_AnonymousSpace.hpp>
#include <Kokkos_LogicalSpaces.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
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
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>

//----------------------------------------------------------------------------

namespace Kokkos {

struct InitArguments {
  int num_threads;
  int num_numa;
  int device_id;
  int ndevices;
  int skip_device;
  bool disable_warnings;
  bool tune_internals;
  bool tool_help        = false;
  std::string tool_lib  = {};
  std::string tool_args = {};

  InitArguments(int nt = -1, int nn = -1, int dv = -1, bool dw = false,
                bool ti = false)
      : num_threads{nt},
        num_numa{nn},
        device_id{dv},
        ndevices{-1},
        skip_device{9999},
        disable_warnings{dw},
        tune_internals{ti} {}
};

namespace Impl {

/* ExecSpaceManager - Responsible for initializing all of the registered
 * backends. Backends are registered using the register_space_initializer()
 * function which should be called from a global context so that it is called
 * prior to initialize_spaces() which is called from Kokkos::initialize()
 */
class ExecSpaceManager {
  std::map<std::string, std::unique_ptr<ExecSpaceInitializerBase>>
      exec_space_factory_list;

 public:
  ExecSpaceManager() = default;

  void register_space_factory(std::string name,
                              std::unique_ptr<ExecSpaceInitializerBase> ptr);
  void initialize_spaces(const Kokkos::InitArguments& args);
  void finalize_spaces(const bool all_spaces);
  void static_fence();
  void print_configuration(std::ostream& msg, const bool detail);
  static ExecSpaceManager& get_instance();
};

template <class SpaceInitializerType>
int initialize_space_factory(std::string name) {
  auto space_ptr = std::make_unique<SpaceInitializerType>();
  ExecSpaceManager::get_instance().register_space_factory(name,
                                                          std::move(space_ptr));
  return 1;
}

}  // namespace Impl
void initialize(int& narg, char* arg[]);

void initialize(InitArguments args = InitArguments());

namespace Impl {

void pre_initialize(const InitArguments& args);

void post_initialize(const InitArguments& args);

void declare_configuration_metadata(const std::string& category,
                                    const std::string& key,
                                    const std::string& value);

}  // namespace Impl

bool is_initialized() noexcept;

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

/** \brief  Finalize all known execution spaces */
void finalize_all();

void fence();

/** \brief Print "Bill of Materials" */
void print_configuration(std::ostream&, const bool detail = false);

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/* Allocate memory from a memory space.
 * The allocation is tracked in Kokkos memory tracking system, so
 * leaked memory can be identified.
 */
template <class Space = typename Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const std::string& arg_alloc_label,
                           const size_t arg_alloc_size) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), arg_alloc_label, arg_alloc_size);
}

template <class Space = typename Kokkos::DefaultExecutionSpace::memory_space>
inline void* kokkos_malloc(const size_t arg_alloc_size) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::allocate_tracked(
      MemorySpace(), "no-label", arg_alloc_size);
}

template <class Space = typename Kokkos::DefaultExecutionSpace::memory_space>
inline void kokkos_free(void* arg_alloc) {
  using MemorySpace = typename Space::memory_space;
  return Impl::SharedAllocationRecord<MemorySpace>::deallocate_tracked(
      arg_alloc);
}

template <class Space = typename Kokkos::DefaultExecutionSpace::memory_space>
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

class ScopeGuard {
 public:
  ScopeGuard(int& narg, char* arg[]) {
    sg_init = false;
    if (!Kokkos::is_initialized()) {
      initialize(narg, arg);
      sg_init = true;
    }
  }

  ScopeGuard(const InitArguments& args = InitArguments()) {
    sg_init = false;
    if (!Kokkos::is_initialized()) {
      initialize(args);
      sg_init = true;
    }
  }

  ~ScopeGuard() {
    if (Kokkos::is_initialized() && sg_init) {
      finalize();
    }
  }

  // private:
  bool sg_init;

  ScopeGuard& operator=(const ScopeGuard&) = delete;
  ScopeGuard(const ScopeGuard&)            = delete;
};

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

// Specializations requires after core definitions
#include <KokkosCore_Config_PostInclude.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
