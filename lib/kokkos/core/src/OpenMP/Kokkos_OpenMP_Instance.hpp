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

  ~OpenMPInternal() { clear_thread_data(); }

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

namespace Experimental {
namespace Impl {
// Partitioning an Execution Space: expects space and integer arguments for
// relative weight
template <typename T>
inline std::vector<OpenMP> create_OpenMP_instances(
    OpenMP const& main_instance, std::vector<T> const& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");
  if (weights.size() == 0) {
    Kokkos::abort("Kokkos::abort: Partition weights vector is empty.");
  }
  std::vector<OpenMP> instances(weights.size());
  double total_weight = std::accumulate(weights.begin(), weights.end(), 0.);
  int const main_pool_size =
      main_instance.impl_internal_space_instance()->thread_pool_size();

  int resources_left = main_pool_size;
  for (unsigned int i = 0; i < weights.size() - 1; ++i) {
    int instance_pool_size = (weights[i] / total_weight) * main_pool_size;
    if (instance_pool_size == 0) {
      Kokkos::abort("Kokkos::abort: Instance has no resource allocated to it");
    }
    instances[i] = OpenMP(instance_pool_size);
    resources_left -= instance_pool_size;
  }
  // Last instance get all resources left
  if (resources_left <= 0) {
    Kokkos::abort(
        "Kokkos::abort: Partition not enough resources left to create the last "
        "instance.");
  }
  instances[weights.size() - 1] = OpenMP(resources_left);

  return instances;
}
}  // namespace Impl

template <typename... Args>
std::vector<OpenMP> partition_space(OpenMP const& main_instance, Args... args) {
  // Unpack the arguments and create the weight vector. Note that if not all of
  // the types are the same, you will get a narrowing warning.
  std::vector<std::common_type_t<Args...>> const weights = {args...};
  return Impl::create_OpenMP_instances(main_instance, weights);
}

template <typename T>
std::vector<OpenMP> partition_space(OpenMP const& main_instance,
                                    std::vector<T> const& weights) {
  return Impl::create_OpenMP_instances(main_instance, weights);
}
}  // namespace Experimental
}  // namespace Kokkos

#endif
