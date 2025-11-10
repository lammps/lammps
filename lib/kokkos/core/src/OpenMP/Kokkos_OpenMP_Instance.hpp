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

#ifndef KOKKOS_OPENMP_INSTANCE_HPP
#define KOKKOS_OPENMP_INSTANCE_HPP

#include <Kokkos_Macros.hpp>
#if !defined(_OPENMP) && !defined(__CUDA_ARCH__) && \
    !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
#error \
    "You enabled Kokkos OpenMP support without enabling OpenMP in the compiler!"
#endif

#include <OpenMP/Kokkos_OpenMP.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>

#include <Kokkos_Atomic.hpp>

#include <impl/Kokkos_ConcurrentBitset.hpp>

#include <omp.h>

#include <mutex>
#include <numeric>
#include <type_traits>
#include <vector>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

class OpenMPInternal;

struct OpenMPTraits {
  static constexpr int MAX_THREAD_COUNT = 512;
};

class OpenMPInternal {
 private:
  OpenMPInternal(int arg_pool_size)
      : m_pool_size{arg_pool_size}, m_level{omp_get_level()}, m_pool() {
    // guard pushing to all_instances
    {
      std::scoped_lock lock(all_instances_mutex);
      all_instances.push_back(this);
    }
  }

  OpenMPInternal()                                 = delete;
  OpenMPInternal(const OpenMPInternal&)            = delete;
  OpenMPInternal& operator=(const OpenMPInternal&) = delete;

  static int get_current_max_threads() noexcept;

  bool m_initialized = false;

  int m_pool_size;
  int m_level;

  HostThreadTeamData* m_pool[OpenMPTraits::MAX_THREAD_COUNT];

 public:
  friend class Kokkos::OpenMP;

  static OpenMPInternal& singleton();

  void initialize(int thread_cound);

  void finalize();

  void clear_thread_data();

  static int max_hardware_threads() noexcept;

  int thread_pool_size() const { return m_pool_size; }

  void resize_thread_data(size_t pool_reduce_bytes, size_t team_reduce_bytes,
                          size_t team_shared_bytes, size_t thread_local_bytes);

  HostThreadTeamData* get_thread_data() const noexcept {
    return m_pool[m_level == omp_get_level() ? 0 : omp_get_thread_num()];
  }

  HostThreadTeamData* get_thread_data(int i) const noexcept {
    return m_pool[i];
  }

  int get_level() const { return m_level; }

  bool is_initialized() const { return m_initialized; }

  bool verify_is_initialized(const char* const label) const;

  void print_configuration(std::ostream& s) const;

  std::mutex m_instance_mutex;

  static std::vector<OpenMPInternal*> all_instances;
  static std::mutex all_instances_mutex;
};

inline bool execute_in_serial(OpenMP const& space = OpenMP()) {
// The default value returned by `omp_get_max_active_levels` with gcc version
// lower than 11.1.0 is 2147483647 instead of 1.
#if (!defined(KOKKOS_COMPILER_GNU) || KOKKOS_COMPILER_GNU >= 1110) && \
    _OPENMP >= 201511
  bool is_nested = omp_get_max_active_levels() > 1;
#else
  bool is_nested = static_cast<bool>(omp_get_nested());
#endif
  return (space.impl_internal_space_instance()->get_level() < omp_get_level() &&
          !(is_nested && (omp_get_level() == 1)));
}

}  // namespace Impl

namespace Experimental::Impl {
// Calculate pool sizes for partitioned OpenMP spaces
template <typename T>
inline std::vector<int> calculate_omp_pool_sizes(
    OpenMP const& main_instance, std::vector<T> const& weights) {
  static_assert(
      std::is_arithmetic_v<T>,
      "Kokkos Error: partitioning arguments must be integers or floats");
  if (weights.size() == 0) {
    Kokkos::abort("Kokkos::abort: Partition weights vector is empty.");
  }
  std::vector<int> pool_sizes(weights.size());
  double total_weight = std::accumulate(weights.begin(), weights.end(), 0.);
  int const main_pool_size =
      main_instance.impl_internal_space_instance()->thread_pool_size();

  int resources_left = main_pool_size;
  for (unsigned int i = 0; i < weights.size() - 1; ++i) {
    int instance_pool_size = (weights[i] / total_weight) * main_pool_size;
    if (instance_pool_size == 0) {
      Kokkos::abort("Kokkos::abort: Instance has no resource allocated to it");
    }
    pool_sizes[i] = instance_pool_size;
    resources_left -= instance_pool_size;
  }
  // Last instance get all resources left
  if (resources_left <= 0) {
    Kokkos::abort(
        "Kokkos::abort: Partition not enough resources left to create the last "
        "instance.");
  }
  pool_sizes[weights.size() - 1] = resources_left;

  return pool_sizes;
}

// Create new OpenMP instances with pool sizes relative to input weights
template <class T>
std::vector<OpenMP> impl_partition_space(const OpenMP& base_instance,
                                         const std::vector<T>& weights) {
  const auto pool_sizes =
      Impl::calculate_omp_pool_sizes(base_instance, weights);

  std::vector<OpenMP> instances;
  instances.reserve(pool_sizes.size());
  for (size_t i = 0; i < pool_sizes.size(); ++i) {
    instances.emplace_back(OpenMP(pool_sizes[i]));
  }

  return instances;
}
}  // namespace Experimental::Impl
}  // namespace Kokkos

#endif
