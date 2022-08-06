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

#ifndef KOKKOS_OPENMPEXEC_HPP
#define KOKKOS_OPENMPEXEC_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_OPENMP)

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPExec;

extern int g_openmp_hardware_max_threads;

extern __thread int t_openmp_hardware_id;
extern __thread OpenMPExec* t_openmp_instance;

//----------------------------------------------------------------------------
/** \brief  Data for OpenMP thread execution */

class OpenMPExec {
 public:
  friend class Kokkos::OpenMP;

  enum { MAX_THREAD_COUNT = 512 };

  void clear_thread_data();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  KOKKOS_DEPRECATED static void validate_partition(const int nthreads,
                                                   int& num_partitions,
                                                   int& partition_size) {
    validate_partition_impl(nthreads, num_partitions, partition_size);
  }
  static void validate_partition_impl(const int nthreads, int& num_partitions,
                                      int& partition_size);
#endif

 private:
  OpenMPExec(int arg_pool_size)
      : m_pool_size{arg_pool_size}, m_level{omp_get_level()}, m_pool() {}

  ~OpenMPExec() { clear_thread_data(); }

  int m_pool_size;
  int m_level;

  HostThreadTeamData* m_pool[MAX_THREAD_COUNT];

 public:
  static void verify_is_master(const char* const);

  void resize_thread_data(size_t pool_reduce_bytes, size_t team_reduce_bytes,
                          size_t team_shared_bytes, size_t thread_local_bytes);

  inline HostThreadTeamData* get_thread_data() const noexcept {
    return m_pool[m_level == omp_get_level() ? 0 : omp_get_thread_num()];
  }

  inline HostThreadTeamData* get_thread_data(int i) const noexcept {
    return m_pool[i];
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline bool OpenMP::impl_is_initialized() noexcept {
  return Impl::t_openmp_instance != nullptr;
}

inline bool OpenMP::in_parallel(OpenMP const&) noexcept {
  // t_openmp_instance is only non-null on a master thread
  return !Impl::t_openmp_instance ||
         Impl::t_openmp_instance->m_level < omp_get_level();
}

inline int OpenMP::impl_thread_pool_size() noexcept {
  return OpenMP::in_parallel() ? omp_get_num_threads()
                               : Impl::t_openmp_instance->m_pool_size;
}

KOKKOS_INLINE_FUNCTION
int OpenMP::impl_thread_pool_rank() noexcept {
  KOKKOS_IF_ON_HOST(
      (return Impl::t_openmp_instance ? 0 : omp_get_thread_num();))

  KOKKOS_IF_ON_DEVICE((return -1;))
}

inline void OpenMP::impl_static_fence(OpenMP const& /**instance*/,
                                      const std::string& name) noexcept {
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
    using Exec = Impl::OpenMPExec;

    Exec* prev_instance = Impl::t_openmp_instance;

    Exec::validate_partition_impl(prev_instance->m_pool_size, num_partitions,
                                  partition_size);

    OpenMP::memory_space space;

#pragma omp parallel num_threads(num_partitions)
    {
      void* ptr = nullptr;
      try {
        ptr = space.allocate(sizeof(Exec));
      } catch (
          Kokkos::Experimental::RawMemoryAllocationFailure const& failure) {
        // For now, just rethrow the error message the existing way
        Kokkos::Impl::throw_runtime_exception(failure.get_error_message());
      }

      Impl::t_openmp_instance = new (ptr) Exec(partition_size);

      size_t pool_reduce_bytes  = 32 * partition_size;
      size_t team_reduce_bytes  = 32 * partition_size;
      size_t team_shared_bytes  = 1024 * partition_size;
      size_t thread_local_bytes = 1024;

      Impl::t_openmp_instance->resize_thread_data(
          pool_reduce_bytes, team_reduce_bytes, team_shared_bytes,
          thread_local_bytes);

      omp_set_num_threads(partition_size);
      f(omp_get_thread_num(), omp_get_num_threads());

      Impl::t_openmp_instance->~Exec();
      space.deallocate(Impl::t_openmp_instance, sizeof(Exec));
      Impl::t_openmp_instance = nullptr;
    }

    Impl::t_openmp_instance = prev_instance;
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
  KOKKOS_INLINE_FUNCTION
  int acquire() const noexcept {
    KOKKOS_IF_ON_HOST((return Kokkos::Impl::t_openmp_hardware_id;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /// \brief release a value acquired by generate
  KOKKOS_INLINE_FUNCTION
  void release(int) const noexcept {}
};

}  // namespace Experimental

inline int OpenMP::impl_thread_pool_size(int depth) {
  return depth < 2 ? impl_thread_pool_size() : 1;
}

KOKKOS_INLINE_FUNCTION
int OpenMP::impl_hardware_thread_id() noexcept {
  KOKKOS_IF_ON_HOST((return Impl::t_openmp_hardware_id;))

  KOKKOS_IF_ON_DEVICE((return -1;))
}

inline int OpenMP::impl_max_hardware_threads() noexcept {
  return Impl::g_openmp_hardware_max_threads;
}

}  // namespace Kokkos

#endif
#endif /* #ifndef KOKKOS_OPENMPEXEC_HPP */
