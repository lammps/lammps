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

#ifndef KOKKOS_EXEC_SPACE_MANAGER_HPP
#define KOKKOS_EXEC_SPACE_MANAGER_HPP

#include <impl/Kokkos_InitializationSettings.hpp>
#include <Kokkos_DetectionIdiom.hpp>
#include <Kokkos_Concepts.hpp>

#include <iosfwd>
#include <map>
#include <string>

namespace {

template <class T>
using public_member_types_t = std::enable_if_t<
    Kokkos::is_execution_space<typename T::execution_space>::value &&
    Kokkos::is_memory_space<typename T::memory_space>::value &&
    Kokkos::is_device<typename T::device_type>::value &&
    Kokkos::is_array_layout<typename T::array_layout>::value &&
    std::is_integral<typename T::size_type>::value &&
    Kokkos::is_memory_space<typename T::scratch_memory_space>::value>;

template <class T>
using print_configuration_t = std::enable_if_t<
    std::is_void<decltype(std::declval<T const&>().print_configuration(
        std::declval<std::ostream&>()))>::value &&
    std::is_void<decltype(std::declval<T const&>().print_configuration(
        std::declval<std::ostream&>(), false))>::value>;

template <class T>
using initialize_finalize_t = std::enable_if_t<
    std::is_void<decltype(T::impl_initialize(
        std::declval<Kokkos::InitializationSettings const&>()))>::value &&
    std::is_void<decltype(T::impl_finalize())>::value>;

template <class T>
using fence_t = std::enable_if_t<
    std::is_void<decltype(std::declval<T const&>().fence())>::value &&
    std::is_void<decltype(std::declval<T const&>().fence("name"))>::value &&
    std::is_void<decltype(T::impl_static_fence("name"))>::value>;

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, "")  // FIXME C++17

template <class ExecutionSpace>
constexpr bool check_valid_execution_space() {
  using Kokkos::is_detected;
  STATIC_ASSERT(std::is_default_constructible<ExecutionSpace>::value);
  STATIC_ASSERT(is_detected<public_member_types_t, ExecutionSpace>::value);
  STATIC_ASSERT(is_detected<print_configuration_t, ExecutionSpace>::value);
  STATIC_ASSERT(is_detected<initialize_finalize_t, ExecutionSpace>::value);
  STATIC_ASSERT(is_detected<fence_t, ExecutionSpace>::value);
#ifndef KOKKOS_ENABLE_HPX  // FIXME_HPX
  STATIC_ASSERT(sizeof(ExecutionSpace) <= 2 * sizeof(void*));
#endif
  return true;
}

#undef STATIC_ASSERT

}  // namespace

namespace Kokkos {
namespace Impl {

struct ExecSpaceBase {
  virtual void initialize(InitializationSettings const&)           = 0;
  virtual void finalize()                                          = 0;
  virtual void static_fence(std::string const&)                    = 0;
  virtual void print_configuration(std::ostream& os, bool verbose) = 0;
  virtual ~ExecSpaceBase()                                         = default;
};

template <class ExecutionSpace>
struct ExecSpaceDerived : ExecSpaceBase {
  static_assert(check_valid_execution_space<ExecutionSpace>(), "");
  void initialize(InitializationSettings const& settings) final {
    ExecutionSpace::impl_initialize(settings);
  }
  void finalize() final { ExecutionSpace::impl_finalize(); }
  void static_fence(std::string const& label) final {
    ExecutionSpace::impl_static_fence(label);
  }
  void print_configuration(std::ostream& os, bool verbose) final {
    ExecutionSpace().print_configuration(os, verbose);
  }
};

/* ExecSpaceManager - Responsible for initializing all the registered
 * backends. Backends are registered using the register_space_initializer()
 * function which should be called from a global context so that it is called
 * prior to initialize_spaces() which is called from Kokkos::initialize()
 */
class ExecSpaceManager {
  std::map<std::string, std::unique_ptr<ExecSpaceBase>> exec_space_factory_list;
  ExecSpaceManager() = default;

 public:
  void register_space_factory(std::string name,
                              std::unique_ptr<ExecSpaceBase> ptr);
  void initialize_spaces(const Kokkos::InitializationSettings& settings);
  void finalize_spaces();
  void static_fence(const std::string&);
  void print_configuration(std::ostream& os, bool verbose);
  static ExecSpaceManager& get_instance();
};

template <class ExecutionSpace>
int initialize_space_factory(std::string name) {
  auto space_ptr = std::make_unique<ExecSpaceDerived<ExecutionSpace>>();
  ExecSpaceManager::get_instance().register_space_factory(name,
                                                          std::move(space_ptr));
  return 1;
}

}  // namespace Impl
}  // namespace Kokkos

#endif
