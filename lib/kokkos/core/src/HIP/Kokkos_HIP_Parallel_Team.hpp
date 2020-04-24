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

#ifndef KOKKO_HIP_PARALLEL_TEAM_HPP
#define KOKKO_HIP_PARALLEL_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_Locks.hpp>
#include <HIP/Kokkos_HIP_Team.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>

namespace Kokkos {
namespace Impl {
template <typename... Properties>
class TeamPolicyInternal<Kokkos::Experimental::HIP, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  template <typename ExecSpace, typename... OtherProperties>
  friend class TeamPolicyInternal;

 private:
  static int constexpr MAX_WARP = 8;

  typename traits::execution_space m_space;
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  int m_team_scratch_size[2];
  int m_thread_scratch_size[2];
  int m_chunk_size;

 public:
  using execution_space = Kokkos::Experimental::HIP;

  template <class... OtherProperties>
  TeamPolicyInternal(TeamPolicyInternal<OtherProperties...> const& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    m_space                  = p.m_space;
  }

  TeamPolicyInternal& operator=(TeamPolicyInternal const& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    m_space                  = p.m_space;
    return *this;
  }

  template <typename FunctorType>
  int team_size_max(FunctorType const& f, ParallelForTag const&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...> >;
    hipFuncAttributes attr = ::Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type,
        typename traits::launch_bounds>::get_hip_func_attributes();
    int const block_size = ::Kokkos::Experimental::Impl::hip_get_max_block_size<
        FunctorType, typename traits::launch_bounds>(
        space().impl_internal_space_instance(), attr, f,
        static_cast<size_t>(vector_length()),
        static_cast<size_t>(team_scratch_size(0)) + 2 * sizeof(double),
        static_cast<size_t>(thread_scratch_size(0)) + sizeof(double));
    return block_size / vector_length();
  }

  template <typename FunctorType>
  int team_size_recommended(FunctorType const& f, ParallelForTag const&) const {
    typedef Impl::ParallelFor<FunctorType, TeamPolicy<Properties...> >
        closure_type;
    hipFuncAttributes attr = ::Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type,
        typename traits::launch_bounds>::get_hip_func_attributes();
    int const block_size = ::Kokkos::Experimental::Impl::hip_get_opt_block_size<
        FunctorType, typename traits::launch_bounds>(
        space().impl_internal_space_instance(), attr, f,
        static_cast<size_t>(vector_length()),
        static_cast<size_t>(team_scratch_size(0)) + 2 * sizeof(double),
        static_cast<size_t>(thread_scratch_size(0)) + sizeof(double));
    return block_size / vector_length();
  }

  static int vector_length_max() {
    return ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
  }

  static int verify_requested_vector_length(int requested_vector_length) {
    int test_vector_length =
        std::min(requested_vector_length, vector_length_max());

    // Allow only power-of-two vector_length
    if (!(is_integral_power_of_two(test_vector_length))) {
      int test_pow2 = 1;
      for (int i = 0; i < 5; i++) {
        test_pow2 = test_pow2 << 1;
        if (test_pow2 > test_vector_length) {
          break;
        }
      }
      test_vector_length = test_pow2 >> 1;
    }

    return test_vector_length;
  }

  static int scratch_size_max(int level) {
    return (
        level == 0 ? 1024 * 40 :  // FIXME_HIP arbitrarily setting this to 48kB
            20 * 1024 * 1024);    // FIXME_HIP arbitrarily setting this to 20MB
  }

  int vector_length() const { return m_vector_length; }

  int team_size() const { return m_team_size; }

  int league_size() const { return m_league_size; }

  int scratch_size(int level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  int team_scratch_size(int level) const { return m_team_scratch_size[level]; }

  int thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  typename traits::execution_space space() const { return m_space; }

  TeamPolicyInternal()
      : m_space(typename traits::execution_space()),
        m_league_size(0),
        m_team_size(-1),
        m_vector_length(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request, int vector_length_request = 1)
      : m_space(space_),
        m_league_size(league_size_),
        m_team_size(team_size_request),
        m_vector_length(verify_requested_vector_length(vector_length_request)),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {
    // Make sure league size is permissable
    if (league_size_ >=
        static_cast<int>(
            ::Kokkos::Experimental::Impl::hip_internal_maximum_grid_count()))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on HIP execution "
          "space.");

    // Make sure total block size is permissable
    if (m_team_size * m_vector_length > 1024) {
      Impl::throw_runtime_exception(
          std::string("Kokkos::TeamPolicy< HIP > the team size is too large. "
                      "Team size x vector length must be smaller than 1024."));
    }
  }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : m_space(space_),
        m_league_size(league_size_),
        m_team_size(-1),
        m_vector_length(verify_requested_vector_length(vector_length_request)),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {
    // Make sure league size is permissable
    if (league_size_ >=
        static_cast<int>(
            ::Kokkos::Experimental::Impl::hip_internal_maximum_grid_count()))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on HIP execution "
          "space.");
  }

  TeamPolicyInternal(int league_size_, int team_size_request,
                     int vector_length_request = 1)
      : m_space(typename traits::execution_space()),
        m_league_size(league_size_),
        m_team_size(team_size_request),
        m_vector_length(verify_requested_vector_length(vector_length_request)),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {
    // Make sure league size is permissable
    if (league_size_ >=
        static_cast<int>(
            ::Kokkos::Experimental::Impl::hip_internal_maximum_grid_count()))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on HIP execution "
          "space.");

    // Make sure total block size is permissable
    if (m_team_size * m_vector_length > 1024) {
      Impl::throw_runtime_exception(
          std::string("Kokkos::TeamPolicy< HIP > the team size is too large. "
                      "Team size x vector length must be smaller than 1024."));
    }
  }

  TeamPolicyInternal(int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : m_space(typename traits::execution_space()),
        m_league_size(league_size_),
        m_team_size(-1),
        m_vector_length(verify_requested_vector_length(vector_length_request)),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(::Kokkos::Experimental::Impl::HIPTraits::WarpSize) {
    // Make sure league size is permissable
    if (league_size_ >=
        static_cast<int>(
            ::Kokkos::Experimental::Impl::hip_internal_maximum_grid_count()))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on HIP execution "
          "space.");
  }

  int chunk_size() const { return m_chunk_size; }

  TeamPolicyInternal& set_chunk_size(typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  TeamPolicyInternal& set_scratch_size(int level,
                                       PerTeamValue const& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  TeamPolicyInternal& set_scratch_size(int level,
                                       PerThreadValue const& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  TeamPolicyInternal& set_scratch_size(int level, PerTeamValue const& per_team,
                                       PerThreadValue const& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  using member_type = Kokkos::Impl::HIPTeamMember;

 protected:
  template <class ClosureType, class FunctorType, class BlockSizeCallable>
  int internal_team_size_common(const FunctorType& f,
                                BlockSizeCallable&& block_size_callable) const {
    using closure_type = ClosureType;
    using functor_value_traits =
        Impl::FunctorValueTraits<FunctorType, typename traits::work_tag>;

    hipFuncAttributes attr = ::Kokkos::Experimental::Impl::HIPParallelLaunch<
        closure_type,
        typename traits::launch_bounds>::get_hip_func_attributes();
    const int block_size = std::forward<BlockSizeCallable>(block_size_callable)(
        space().impl_internal_space_instance(), attr, f,
        static_cast<size_t>(vector_length()),
        static_cast<size_t>(team_scratch_size(0)) + 2 * sizeof(double),
        static_cast<size_t>(thread_scratch_size(0)) + sizeof(double) +
            ((functor_value_traits::StaticValueSize != 0)
                 ? 0
                 : functor_value_traits::value_size(f)));
    KOKKOS_ASSERT(block_size > 0);

    // Currently we require Power-of-2 team size for reductions.
    int p2 = 1;
    while (p2 <= block_size) p2 *= 2;
    p2 /= 2;
    return p2 / vector_length();
  }

  template <class ClosureType, class FunctorType>
  int internal_team_size_max(const FunctorType& f) const {
    return internal_team_size_common<ClosureType>(
        f, ::Kokkos::Experimental::Impl::hip_get_max_block_size<
               FunctorType, typename traits::launch_bounds>);
  }

  template <class ClosureType, class FunctorType>
  int internal_team_size_recommended(const FunctorType& f) const {
    return internal_team_size_common<ClosureType>(
        f, ::Kokkos::Experimental::Impl::hip_get_opt_block_size<
               FunctorType, typename traits::launch_bounds>);
  }
};

struct HIPLockArrays {
  std::int32_t* atomic  = nullptr;
  std::int32_t* scratch = nullptr;
  std::int32_t n        = 0;
};

template <typename FunctorType, typename... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Experimental::HIP> {
 public:
  using Policy = TeamPolicyInternal<Kokkos::Experimental::HIP, Properties...>;
  using functor_type = FunctorType;
  using size_type    = ::Kokkos::Experimental::HIP::size_type;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  // Algorithmic constraints: hipBlockDim_y is a power of two AND hipBlockDim_y
  // == hipBlockDim_z == 1 shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]

  FunctorType const m_functor;
  Policy const m_policy;
  size_type const m_league_size;
  int m_team_size;
  size_type const m_vector_size;
  int m_shmem_begin;
  int m_shmem_size;
  void* m_scratch_ptr[2];
  int m_scratch_size[2];
  mutable HIPLockArrays hip_lock_arrays;

  template <typename TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_team(const Member& member) const {
    m_functor(member);
  }

  template <typename TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_team(const Member& member) const {
    m_functor(TagType(), member);
  }

 public:
  __device__ inline void operator()(void) const {
    // Iterate this block through the league
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      __shared__ int64_t base_thread_id;
      if (hipThreadIdx_x == 0 && hipThreadIdx_y == 0) {
        threadid = (hipBlockIdx_x * hipBlockDim_z + hipThreadIdx_z) %
                   (hip_lock_arrays.n / (hipBlockDim_x * hipBlockDim_y));
        threadid *= hipBlockDim_x * hipBlockDim_y;
        int done = 0;
        while (!done) {
          done = (0 == atomicCAS(&hip_lock_arrays.scratch[threadid], 0, 1));
          if (!done) {
            threadid += hipBlockDim_x * hipBlockDim_y;
            if (int64_t(threadid + hipBlockDim_x * hipBlockDim_y) >=
                int64_t(hip_lock_arrays.n))
              threadid = 0;
          }
        }
        base_thread_id = threadid;
      }
      __syncthreads();
      threadid = base_thread_id;
    }

    int const int_league_size = static_cast<int>(m_league_size);
    for (int league_rank = hipBlockIdx_x; league_rank < int_league_size;
         league_rank += hipGridDim_x) {
      this->template exec_team<WorkTag>(typename Policy::member_type(
          ::Kokkos::Experimental::kokkos_impl_hip_shared_memory<void>(),
          m_shmem_begin, m_shmem_size,
          static_cast<void*>(
              static_cast<char*>(m_scratch_ptr[1]) +
              ptrdiff_t(threadid / (hipBlockDim_x * hipBlockDim_y)) *
                  m_scratch_size[1]),
          m_scratch_size[1], league_rank, m_league_size));
    }
    if (m_scratch_size[1] > 0) {
      __syncthreads();
      if (hipThreadIdx_x == 0 && hipThreadIdx_y == 0)
        hip_lock_arrays.scratch[threadid] = 0;
    }
  }

  inline void execute() const {
    HIP_SAFE_CALL(hipMalloc(
        &hip_lock_arrays.atomic,
        sizeof(std::int32_t) * (KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1)));
    HIP_SAFE_CALL(hipMalloc(
        &hip_lock_arrays.scratch,
        sizeof(std::int32_t) * (::Kokkos::Experimental::HIP::concurrency())));
    HIP_SAFE_CALL(hipMemset(
        hip_lock_arrays.scratch, 0,
        sizeof(std::int32_t) * (::Kokkos::Experimental::HIP::concurrency())));
    hip_lock_arrays.n = ::Kokkos::Experimental::HIP::concurrency();

    int64_t const shmem_size_total = m_shmem_begin + m_shmem_size;
    dim3 const grid(static_cast<int>(m_league_size), 1, 1);
    dim3 const block(static_cast<int>(m_vector_size),
                     static_cast<int>(m_team_size), 1);

    ::Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
        *this, grid, block, shmem_size_total,
        m_policy.space().impl_internal_space_instance(),
        true);  // copy to device and execute

    if (hip_lock_arrays.atomic) {
      HIP_SAFE_CALL(hipFree(hip_lock_arrays.atomic));
      hip_lock_arrays.atomic = nullptr;
    }
    if (hip_lock_arrays.scratch) {
      HIP_SAFE_CALL(hipFree(hip_lock_arrays.scratch));
      hip_lock_arrays.scratch = nullptr;
    }
    hip_lock_arrays.n = 0;
  }

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.vector_length()) {
    hipFuncAttributes attr = ::Kokkos::Experimental::Impl::HIPParallelLaunch<
        ParallelFor, LaunchBounds>::get_hip_func_attributes();
    m_team_size =
        m_team_size >= 0
            ? m_team_size
            : ::Kokkos::Experimental::Impl::hip_get_opt_block_size<
                  FunctorType, LaunchBounds>(
                  m_policy.space().impl_internal_space_instance(), attr,
                  m_functor, m_vector_size, m_policy.team_scratch_size(0),
                  m_policy.thread_scratch_size(0)) /
                  m_vector_size;

    m_shmem_begin = (sizeof(double) * (m_team_size + 2));
    m_shmem_size =
        (m_policy.scratch_size(0, m_team_size) +
         FunctorTeamShmemSize<FunctorType>::value(m_functor, m_team_size));
    m_scratch_size[0] = m_policy.scratch_size(0, m_team_size);
    m_scratch_size[1] = m_policy.scratch_size(1, m_team_size);

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    m_scratch_ptr[0] = nullptr;
    m_scratch_ptr[1] =
        m_team_size <= 0
            ? nullptr
            : ::Kokkos::Experimental::Impl::hip_resize_scratch_space(
                  static_cast<ptrdiff_t>(m_scratch_size[1]) *
                  static_cast<ptrdiff_t>(
                      ::Kokkos::Experimental::HIP::concurrency() /
                      (m_team_size * m_vector_size)));

    int const shmem_size_total = m_shmem_begin + m_shmem_size;
    if (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
        shmem_size_total) {
      printf(
          "%i %i\n",
          m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock,
          shmem_size_total);
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > insufficient shared memory"));
    }

    if (static_cast<int>(m_team_size) >
        static_cast<int>(
            ::Kokkos::Experimental::Impl::hip_get_max_block_size<FunctorType,
                                                                 LaunchBounds>(
                m_policy.space().impl_internal_space_instance(), attr,
                arg_functor, arg_policy.vector_length(),
                arg_policy.team_scratch_size(0),
                arg_policy.thread_scratch_size(0)) /
            arg_policy.vector_length())) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > requested too large team size."));
    }
  }
};
}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
