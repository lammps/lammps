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

#include <Kokkos_OpenMP.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>

#include <Kokkos_Atomic.hpp>

#include <Kokkos_UniqueToken.hpp>
#include <impl/Kokkos_ConcurrentBitset.hpp>

#include <omp.h>

#include <mutex>
#include <numeric>
#include <type_traits>
#include <vector>

namespace Kokkos {
namespace Impl {

class OpenMPInternal;

inline int g_openmp_hardware_max_threads = 1;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
// FIXME_OPENMP we can remove this after we remove partition_master
inline thread_local OpenMPInternal* t_openmp_instance = nullptr;
#endif

struct OpenMPTraits {
  static int constexpr MAX_THREAD_COUNT = 512;
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
inline bool OpenMP::impl_is_initialized() noexcept {
  return Impl::OpenMPInternal::singleton().is_initialized();
}

inline bool OpenMP::in_parallel(OpenMP const& exec_space) noexcept {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  return (
      (exec_space.impl_internal_space_instance()->m_level < omp_get_level()) &&
      (!Impl::t_openmp_instance ||
       Impl::t_openmp_instance->m_level < omp_get_level()));
#else
  return exec_space.impl_internal_space_instance()->m_level < omp_get_level();
#endif
}

inline int OpenMP::impl_thread_pool_size(OpenMP const& exec_space) noexcept {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  return OpenMP::in_parallel(exec_space)
             ? omp_get_num_threads()
             : (Impl::t_openmp_instance
                    ? Impl::t_openmp_instance->m_pool_size
                    : exec_space.impl_internal_space_instance()->m_pool_size);
#else
  return OpenMP::in_parallel(exec_space)
             ? omp_get_num_threads()
             : exec_space.impl_internal_space_instance()->m_pool_size;
#endif
}

inline int OpenMP::impl_thread_pool_rank() noexcept {
  // FIXME_OPENMP Can we remove this when removing partition_master? It's only
  // used in one partition_master test
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  KOKKOS_IF_ON_HOST(
      (return Impl::t_openmp_instance ? 0 : omp_get_thread_num();))
#else
  KOKKOS_IF_ON_HOST((return omp_get_thread_num();))
#endif

  KOKKOS_IF_ON_DEVICE((return -1;))
}

inline void OpenMP::impl_static_fence(std::string const& name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::OpenMP>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      []() {});
}

inline bool OpenMP::is_asynchronous(OpenMP const& /*instance*/) noexcept {
  return false;
}

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

template <>
class UniqueToken<OpenMP, UniqueTokenScope::Instance> {
 private:
  using buffer_type = Kokkos::View<uint32_t*, Kokkos::HostSpace>;
  int m_count;
  buffer_type m_buffer_view;
  uint32_t volatile* m_buffer;

 public:
  using execution_space = OpenMP;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept
      : m_count(::Kokkos::OpenMP::impl_thread_pool_size()),
        m_buffer_view(buffer_type()),
        m_buffer(nullptr) {}

  UniqueToken(size_type max_size, execution_space const& = execution_space())
      : m_count(max_size),
        m_buffer_view("UniqueToken::m_buffer_view",
                      ::Kokkos::Impl::concurrent_bitset::buffer_bound(m_count)),
        m_buffer(m_buffer_view.data()) {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept {
    KOKKOS_IF_ON_HOST((return m_count;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief acquire value such that 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST(
        (if (m_count >= ::Kokkos::OpenMP::impl_thread_pool_size()) return ::
             Kokkos::OpenMP::impl_thread_pool_rank();
         const ::Kokkos::pair<int, int> result =
             ::Kokkos::Impl::concurrent_bitset::acquire_bounded(
                 m_buffer, m_count, ::Kokkos::Impl::clock_tic() % m_count);

         if (result.first < 0) {
           ::Kokkos::abort(
               "UniqueToken<OpenMP> failure to acquire tokens, no tokens "
               "available");
         }

         return result.first;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int i) const noexcept {
    KOKKOS_IF_ON_HOST(
        (if (m_count < ::Kokkos::OpenMP::impl_thread_pool_size()) {
          ::Kokkos::Impl::concurrent_bitset::release(m_buffer, i);
        }))

    KOKKOS_IF_ON_DEVICE(((void)i;))
  }
};

template <>
class UniqueToken<OpenMP, UniqueTokenScope::Global> {
 public:
  using execution_space = OpenMP;
  using size_type       = int;

  /// \brief create object size for concurrency on the given instance
  ///
  /// This object should not be shared between instances
  UniqueToken(execution_space const& = execution_space()) noexcept {}

  /// \brief upper bound for acquired values, i.e. 0 <= value < size()
  KOKKOS_INLINE_FUNCTION
  int size() const noexcept {
    KOKKOS_IF_ON_HOST((return Kokkos::Impl::g_openmp_hardware_max_threads;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief acquire value such that 0 <= value < size()
  // FIXME this is wrong when using nested parallelism. In that case multiple
  // threads have the same thread ID.
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((return omp_get_thread_num();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

}  // namespace Experimental

inline int OpenMP::impl_thread_pool_size(int depth, OpenMP const& exec_space) {
  return depth < 2 ? impl_thread_pool_size(exec_space) : 1;
}

KOKKOS_INLINE_FUNCTION
int OpenMP::impl_hardware_thread_id() noexcept {
  KOKKOS_IF_ON_HOST((return omp_get_thread_num();))

  KOKKOS_IF_ON_DEVICE((return -1;))
}

inline int OpenMP::impl_max_hardware_threads() noexcept {
  return Impl::g_openmp_hardware_max_threads;
}

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
                                    std::vector<T>& weights) {
  return Impl::create_OpenMP_instances(main_instance, weights);
}
}  // namespace Experimental
}  // namespace Kokkos

#endif
