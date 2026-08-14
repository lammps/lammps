// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_EXEC_SPACE_MANAGER_HPP
#define KOKKOS_EXEC_SPACE_MANAGER_HPP

#include <impl/Kokkos_InitializationSettings.hpp>
#include <Kokkos_DetectionIdiom.hpp>
#include <Kokkos_Concepts.hpp>

#include <concepts>
#include <iosfwd>
#include <map>
#include <string>
#include <utility>

namespace Kokkos::Impl {

template <class T>
using public_member_types_t = std::enable_if_t<
    Kokkos::is_execution_space_v<typename T::execution_space> &&
    Kokkos::is_memory_space_v<typename T::memory_space> &&
    Kokkos::is_device_v<typename T::device_type> &&
    Kokkos::is_array_layout_v<typename T::array_layout> &&
    std::unsigned_integral<typename T::size_type> &&
    Kokkos::is_memory_space_v<typename T::scratch_memory_space>>;

template <class T>
using print_configuration_t = std::enable_if_t<
    std::is_void_v<decltype(std::declval<T const&>().print_configuration(
        std::declval<std::ostream&>()))> &&
    std::is_void_v<decltype(std::declval<T const&>().print_configuration(
        std::declval<std::ostream&>(), false))>>;

template <class T>
using initialize_finalize_t = std::enable_if_t<
    std::is_void_v<decltype(T::impl_initialize(
        std::declval<Kokkos::InitializationSettings const&>()))> &&
    std::is_void_v<decltype(T::impl_finalize())>>;

template <class T>
using fence_t = std::enable_if_t<
    std::is_void_v<decltype(std::declval<T const&>().fence())> &&
    std::is_void_v<decltype(std::declval<T const&>().fence("name"))> &&
    std::is_void_v<decltype(T::impl_static_fence("name"))>>;

template <class T>
using concurrency_t = std::enable_if_t<
    std::is_same_v<int, decltype(std::declval<T const&>().concurrency())>>;

template <class T>
constexpr bool check_is_semiregular() {
  static_assert(std::is_default_constructible_v<T>);
  static_assert(std::is_copy_constructible_v<T>);
  static_assert(std::is_move_constructible_v<T>);
  static_assert(std::is_copy_assignable_v<T>);
  static_assert(std::is_move_assignable_v<T>);
  static_assert(std::is_destructible_v<T>);
  return true;
}

template <class T>
using equal_to_t =
    decltype(std::declval<T const&>() == std::declval<T const&>());

template <class T>
using not_equal_to_t =
    decltype(std::declval<T const&>() != std::declval<T const&>());

template <class T>
constexpr bool check_is_equality_comparable() {
  using Kokkos::is_detected_exact_v;
  static_assert(is_detected_exact_v<bool, equal_to_t, T>);
  static_assert(is_detected_exact_v<bool, not_equal_to_t, T>);
  return true;
}

template <class T>
constexpr bool check_is_regular() {
  static_assert(check_is_semiregular<T>() && check_is_equality_comparable<T>());
  return true;
}

template <class ExecutionSpace>
constexpr bool check_valid_execution_space() {
  using Kokkos::is_detected_v;
  static_assert(std::is_default_constructible_v<ExecutionSpace>);
  static_assert(is_detected_v<public_member_types_t, ExecutionSpace>);
  static_assert(is_detected_v<print_configuration_t, ExecutionSpace>);
  static_assert(is_detected_v<initialize_finalize_t, ExecutionSpace>);
  static_assert(is_detected_v<fence_t, ExecutionSpace>);
  static_assert(is_detected_v<concurrency_t, ExecutionSpace>);
  static_assert(sizeof(ExecutionSpace) <= 2 * sizeof(void*));
  return true;
}

struct ExecSpaceBase {
  virtual void initialize(InitializationSettings const&)           = 0;
  virtual void finalize()                                          = 0;
  virtual void static_fence(std::string const&)                    = 0;
  virtual void print_configuration(std::ostream& os, bool verbose) = 0;
  virtual ~ExecSpaceBase()                                         = default;
};

template <class ExecutionSpace>
struct ExecSpaceDerived : ExecSpaceBase {
  static_assert(check_valid_execution_space<ExecutionSpace>());
  static_assert(check_is_regular<ExecutionSpace>());
  void initialize(InitializationSettings const& settings) override final {
    ExecutionSpace::impl_initialize(settings);
  }
  void finalize() override final { ExecutionSpace::impl_finalize(); }
  void static_fence(std::string const& label) override final {
    ExecutionSpace::impl_static_fence(label);
  }
  void print_configuration(std::ostream& os, bool verbose) override final {
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

}  // namespace Kokkos::Impl

#endif
