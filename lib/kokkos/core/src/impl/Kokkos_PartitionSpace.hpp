// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_PARTITION_SPACE_HPP
#define KOKKOS_PARTITION_SPACE_HPP

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PARTITION_SPACE
#endif

#include <Kokkos_Concepts.hpp>

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

namespace Kokkos::Experimental::Impl {

// Customization point for backends. Default behavior is to return the passed
// in instance, ignoring weights
template <class ExecSpace, class T>
std::vector<ExecSpace> impl_partition_space(const ExecSpace& base_instance,
                                            const std::vector<T>& weights) {
  std::vector<ExecSpace> instances;
  instances.reserve(weights.size());
  std::generate_n(std::back_inserter(instances), weights.size(),
                  [&base_instance]() { return base_instance; });

  return instances;
}

}  // namespace Kokkos::Experimental::Impl

namespace Kokkos::Experimental {

// Partitioning an Execution Space
// Input:
//   - Base execution space
//   - integer arguments for relative weight, either input per weight or vector
//   of weights
// Ouput:
//   - Array (or vector) of execution spaces partitioned based on weights
template <class ExecSpace, class... Args>
std::array<ExecSpace, sizeof...(Args)> partition_space(
    ExecSpace const& base_instance, Args... args) {
  static_assert(is_execution_space<ExecSpace>::value,
                "Kokkos Error: partition_space expects an Execution Space as "
                "first argument");
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");

  // Get vector of instances from backend specific impl
  std::vector<std::common_type_t<Args...>> weights = {args...};
  auto instances_vec = Impl::impl_partition_space(base_instance, weights);

  // Convert to std::array and return
  std::array<ExecSpace, sizeof...(Args)> instances;
  std::copy(instances_vec.begin(), instances_vec.end(), instances.begin());
  return instances;
}

template <class ExecSpace, class T>
std::vector<ExecSpace> partition_space(ExecSpace const& base_instance,
                                       std::vector<T> const& weights) {
  static_assert(is_execution_space<ExecSpace>::value,
                "Kokkos Error: partition_space expects an Execution Space as "
                "first argument");
  static_assert(
      std::is_arithmetic_v<T>,
      "Kokkos Error: partitioning arguments must be integers or floats");

  // Return vector of instances from backend specific impl
  return Impl::impl_partition_space(base_instance, weights);
}

}  // namespace Kokkos::Experimental

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PARTITION_SPACE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_PARTITION_SPACE
#endif

#endif
