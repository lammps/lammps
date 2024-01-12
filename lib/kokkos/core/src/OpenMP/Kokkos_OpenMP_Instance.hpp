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

inline bool execute_in_serial(OpenMP const& space = OpenMP()) {
  return (OpenMP::in_parallel(space) &&
          !(omp_get_nested() && (omp_get_level() == 1)));
}

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

class OpenMPInternal;

inline int g_openmp_hardware_max_threads = 1;

struct OpenMPTraits {
  static constexpr int MAX_THREAD_COUNT = 512;
};

class OpenMPInternal {
 private:
  OpenMPInternal(int arg_pool_size)
      : m_pool_size{arg_pool_size}, m_level{omp_get_level()}, m_pool() {}

  ~OpenMPInternal() { clear_thread_data(); }

  static int get_current_max_threads() noexcept;

  bool m_initialized = false;

  int m_pool_size;
  int m_level;
  int m_pool_mutex = 0;

  HostThreadTeamData* m_pool[OpenMPTraits::MAX_THREAD_COUNT];

 public:
  friend class Kokkos::OpenMP;

  static OpenMPInternal& singleton();

  void initialize(int thread_cound);

  void finalize();

  void clear_thread_data();

  int thread_pool_size() const { return m_pool_size; }

  // Acquire lock used to protect access to m_pool
  void acquire_lock();

  // Release lock used to protect access to m_pool
  void release_lock();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  static void validate_partition_impl(const int nthreads, int& num_partitions,
                                      int& partition_size);
#endif

  void resize_thread_data(size_t pool_reduce_bytes, size_t team_reduce_bytes,
                          size_t team_shared_bytes, size_t thread_local_bytes);

  HostThreadTeamData* get_thread_data() const noexcept {
    return m_pool[m_level == omp_get_level() ? 0 : omp_get_thread_num()];
  }

  HostThreadTeamData* get_thread_data(int i) const noexcept {
    return m_pool[i];
  }

  bool is_initialized() const { return m_initialized; }

  bool verify_is_initialized(const char* const label) const;

  void print_configuration(std::ostream& s) const;
};

}  // namespace Impl

namespace Experimental {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
template <>
class MasterLock<OpenMP> {
 public:
  void lock() { omp_set_lock(&m_lock); }
  void unlock() { omp_unset_lock(&m_lock); }
  bool try_lock() { return static_cast<bool>(omp_test_lock(&m_lock)); }

  KOKKOS_DEPRECATED MasterLock() { omp_init_lock(&m_lock); }
  ~MasterLock() { omp_destroy_lock(&m_lock); }

  MasterLock(MasterLock const&) = delete;
  MasterLock(MasterLock&&)      = delete;
  MasterLock& operator=(MasterLock const&) = delete;
  MasterLock& operator=(MasterLock&&) = delete;

 private:
  omp_lock_t m_lock;
};
#endif

}  // namespace Experimental

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
  instances[weights.size() - 1] = resources_left;

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

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
template <typename F>
KOKKOS_DEPRECATED void OpenMP::partition_master(F const& f, int num_partitions,
                                                int partition_size) {
#if _OPENMP >= 201511
  if (omp_get_max_active_levels() > 1) {
#else
  if (omp_get_nested()) {
#endif
    using Exec = Impl::OpenMPInternal;

    Exec* prev_instance = &Impl::OpenMPInternal::singleton();

    Exec::validate_partition_impl(prev_instance->m_pool_size, num_partitions,
                                  partition_size);

    OpenMP::memory_space space;

#pragma omp parallel num_threads(num_partitions)
    {
      Exec thread_local_instance(partition_size);
      Impl::t_openmp_instance = &thread_local_instance;

      size_t pool_reduce_bytes  = 32 * partition_size;
      size_t team_reduce_bytes  = 32 * partition_size;
      size_t team_shared_bytes  = 1024 * partition_size;
      size_t thread_local_bytes = 1024;

      thread_local_instance.resize_thread_data(
          pool_reduce_bytes, team_reduce_bytes, team_shared_bytes,
          thread_local_bytes);

      omp_set_num_threads(partition_size);
      f(omp_get_thread_num(), omp_get_num_threads());
      Impl::t_openmp_instance = nullptr;
    }
  } else {
    // nested openmp not enabled
    f(0, 1);
  }
}
#endif

}  // namespace Kokkos

#endif
