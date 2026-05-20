// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CHECK_USAGE_HPP
#define KOKKOS_CHECK_USAGE_HPP

#include <Kokkos_Abort.hpp>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_InitializeFinalize.hpp>
#include <impl/Kokkos_Utilities.hpp>

#include <sstream>
#include <type_traits>

// FIXME: Obtain file and line number information via std::source_location
// (since C++20) which requires GCC 11 etc.

namespace Kokkos {

template <typename... Args>
class RangePolicy;

template <typename... Args>
struct MDRangePolicy;

template <typename... Args>
class TeamPolicy;

namespace Impl {

template <typename T>
const char* fetch_policy_name(const T&) {
  if ((std::is_same_v<T, std::size_t>) ||
      (Kokkos::Impl::is_specialization_of_v<T, Kokkos::RangePolicy>))
    return "RangePolicy";
  else if (Kokkos::Impl::is_specialization_of_v<T, Kokkos::MDRangePolicy>)
    return "MDRangePolicy";
  else if (Kokkos::Impl::is_specialization_of_v<T, Kokkos::TeamPolicy>)
    return "TeamPolicy";
  else
    return "Probably nested execution policy";
}

// The types CheckUsage and UsageRequires are for specifying run-time checks
// of correct usage of Kokkos functions and other constructs.
template <typename... T>
struct CheckUsage;

struct UsageRequires {
  struct isInitialized {};
  struct isNotFinalized {};
  struct insideExecEnv {};
};

template <>
struct CheckUsage<UsageRequires::isInitialized> {
  template <typename T>
  static void check(const char* func_name, const T& exec_policy,
                    const char* meta_data = "no-label") {
    if (!Kokkos::is_initialized()) {
      std::stringstream ss;
      ss << "Kokkos ERROR: attempting to call " << func_name << "() "
         << "**before** Kokkos::initialize() was called."
         << " Concerns " << meta_data << " with exec policy "
         << fetch_policy_name(exec_policy) << ".";
      Kokkos::abort(ss.str().c_str());
    }
  }
};

template <>
struct CheckUsage<UsageRequires::isNotFinalized> {
  template <typename T>
  static void check(const char* func_name, const T& exec_policy,
                    const char* meta_data = "no-label") {
    if (Kokkos::is_finalized()) {
      std::stringstream ss;
      ss << "Kokkos ERROR: attempting to call " << func_name << "() "
         << "**after** Kokkos::finalize() was called."
         << " Concerns " << meta_data << " with exec policy "
         << fetch_policy_name(exec_policy) << ".";
      Kokkos::abort(ss.str().c_str());
    }
  }
};

template <>
struct CheckUsage<UsageRequires::insideExecEnv> {
  template <typename T>
  static void check(const char* func_name, const T& exec_policy,
                    const char* meta_data = "no-label") {
    Kokkos::Impl::CheckUsage<
        Kokkos::Impl::UsageRequires::isNotFinalized>::check(func_name,
                                                            exec_policy,
                                                            meta_data);
    Kokkos::Impl::CheckUsage<Kokkos::Impl::UsageRequires::isInitialized>::check(
        func_name, exec_policy, meta_data);
  }
};

// NOLINTBEGIN(bugprone-exception-escape)
inline void check_execution_space_constructor_precondition(
    char const* name) noexcept {
  if (Kokkos::is_finalized()) {
    std::stringstream err;
    err << "Kokkos ERROR: " << name
        << " execution space is being constructed"
           " after finalize() has been called";
    Kokkos::abort(err.str().c_str());
  }
  if (!Kokkos::is_initialized()) {
    std::stringstream err;
    err << "Kokkos ERROR: " << name
        << " execution space is being constructed"
           " before initialize() has been called";
    Kokkos::abort(err.str().c_str());
  }
}

inline void check_execution_space_destructor_precondition(
    char const* name) noexcept {
  if (Kokkos::is_finalized()) {
    std::stringstream err;
    err << "Kokkos ERROR: " << name
        << " execution space is being destructed"
           " after finalize() has been called";
    Kokkos::abort(err.str().c_str());
  }
}
// NOLINTEND(bugprone-exception-escape)

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CHECK_USAGE_HPP
