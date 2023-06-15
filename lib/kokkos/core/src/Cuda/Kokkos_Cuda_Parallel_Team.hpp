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

#ifndef KOKKOS_CUDA_PARALLEL_TEAM_HPP
#define KOKKOS_CUDA_PARALLEL_TEAM_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdint>

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <Cuda/Kokkos_Cuda_Locks.hpp>
#include <Cuda/Kokkos_Cuda_Team.hpp>
#include <Kokkos_MinMaxClamp.hpp>
#include <Kokkos_Vectorization.hpp>

#include <impl/Kokkos_Tools.hpp>
#include <typeinfo>

#include <impl/KokkosExp_IterateTileGPU.hpp>

namespace Kokkos {

extern bool show_warnings() noexcept;

namespace Impl {

template <class... Properties>
class TeamPolicyInternal<Kokkos::Cuda, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  template <class ExecSpace, class... OtherProperties>
  friend class TeamPolicyInternal;

 private:
  static constexpr int MAX_WARP = 8;

  typename traits::execution_space m_space;
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];
  int m_chunk_size;
  bool m_tune_team;
  bool m_tune_vector;

 public:
  //! Execution space of this execution policy
  using execution_space = Kokkos::Cuda;

  template <class... OtherProperties>
  TeamPolicyInternal(const TeamPolicyInternal<OtherProperties...>& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    m_space                  = p.m_space;
    m_tune_team              = p.m_tune_team;
    m_tune_vector            = p.m_tune_vector;
  }

  //----------------------------------------

  template <class FunctorType>
  int team_size_max(const FunctorType& f, const ParallelForTag&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...>>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes();
    int block_size =
        Kokkos::Impl::cuda_get_max_block_size<FunctorType,
                                              typename traits::launch_bounds>(
            space().impl_internal_space_instance(), attr, f,
            (size_t)impl_vector_length(),
            (size_t)team_scratch_size(0) + 2 * sizeof(double),
            (size_t)thread_scratch_size(0) + sizeof(double));
    return block_size / impl_vector_length();
  }

  template <class FunctorType>
  inline int team_size_max(const FunctorType& f,
                           const ParallelReduceTag&) const {
    using functor_analysis_type =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                              TeamPolicyInternal, FunctorType>;
    using reducer_type = typename Impl::ParallelReduceReturnValue<
        void, typename functor_analysis_type::value_type,
        FunctorType>::reducer_type;
    using closure_type =
        Impl::ParallelReduce<FunctorType, TeamPolicy<Properties...>,
                             reducer_type>;
    return internal_team_size_max<closure_type>(f);
  }

  template <class FunctorType, class ReducerType>
  inline int team_size_max(const FunctorType& f, const ReducerType& /*r*/,
                           const ParallelReduceTag&) const {
    using closure_type =
        Impl::ParallelReduce<FunctorType, TeamPolicy<Properties...>,
                             ReducerType>;
    return internal_team_size_max<closure_type>(f);
  }

  template <class FunctorType>
  int team_size_recommended(const FunctorType& f, const ParallelForTag&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...>>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes();
    const int block_size =
        Kokkos::Impl::cuda_get_opt_block_size<FunctorType,
                                              typename traits::launch_bounds>(
            space().impl_internal_space_instance(), attr, f,
            (size_t)impl_vector_length(),
            (size_t)team_scratch_size(0) + 2 * sizeof(double),
            (size_t)thread_scratch_size(0) + sizeof(double));
    return block_size / impl_vector_length();
  }

  template <class FunctorType>
  inline int team_size_recommended(const FunctorType& f,
                                   const ParallelReduceTag&) const {
    using functor_analysis_type =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                              TeamPolicyInternal, FunctorType>;
    using reducer_type = typename Impl::ParallelReduceReturnValue<
        void, typename functor_analysis_type::value_type,
        FunctorType>::reducer_type;
    using closure_type =
        Impl::ParallelReduce<FunctorType, TeamPolicy<Properties...>,
                             reducer_type>;
    return internal_team_size_recommended<closure_type>(f);
  }

  template <class FunctorType, class ReducerType>
  int team_size_recommended(const FunctorType& f, const ReducerType&,
                            const ParallelReduceTag&) const {
    using closure_type =
        Impl::ParallelReduce<FunctorType, TeamPolicy<Properties...>,
                             ReducerType>;
    return internal_team_size_recommended<closure_type>(f);
  }

  inline static int vector_length_max() { return Impl::CudaTraits::WarpSize; }

  inline static int verify_requested_vector_length(
      int requested_vector_length) {
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

  inline static int scratch_size_max(int level) {
    return (
        level == 0 ? 1024 * 40 :  // 48kB is the max for CUDA, but we need some
                                  // for team_member.reduce etc.
            20 * 1024 *
                1024);  // arbitrarily setting this to 20MB, for a Volta V100
                        // that would give us about 3.2GB for 2 teams per SM
  }

  //----------------------------------------

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  KOKKOS_DEPRECATED inline int vector_length() const {
    return impl_vector_length();
  }
#endif
  inline int impl_vector_length() const { return m_vector_length; }
  inline int team_size() const { return m_team_size; }
  inline int league_size() const { return m_league_size; }
  inline bool impl_auto_team_size() const { return m_tune_team; }
  inline bool impl_auto_vector_length() const { return m_tune_vector; }
  inline void impl_set_team_size(size_t team_size) { m_team_size = team_size; }
  inline void impl_set_vector_length(size_t vector_length) {
    m_vector_length = vector_length;
  }
  size_t scratch_size(int level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }
  size_t team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }
  size_t thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  const typename traits::execution_space& space() const { return m_space; }

  TeamPolicyInternal()
      : m_space(typename traits::execution_space()),
        m_league_size(0),
        m_team_size(-1),
        m_vector_length(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(Impl::CudaTraits::WarpSize),
        m_tune_team(false),
        m_tune_vector(false) {}

  /** \brief  Specify league size, specify team size, specify vector length */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request, int vector_length_request = 1)
      : m_space(space_),
        m_league_size(league_size_),
        m_team_size(team_size_request),
        m_vector_length(
            (vector_length_request > 0)
                ? verify_requested_vector_length(vector_length_request)
                : verify_requested_vector_length(1)),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(Impl::CudaTraits::WarpSize),
        m_tune_team(bool(team_size_request <= 0)),
        m_tune_vector(bool(vector_length_request <= 0)) {
    // Make sure league size is permissible
    if (league_size_ >= int(Impl::cuda_internal_maximum_grid_count()[0]))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on Cuda execution "
          "space.");

    // Make sure total block size is permissible
    if (m_team_size * m_vector_length >
        int(Impl::CudaTraits::MaxHierarchicalParallelism)) {
      Impl::throw_runtime_exception(
          std::string("Kokkos::TeamPolicy< Cuda > the team size is too large. "
                      "Team size x vector length must be smaller than 1024."));
    }
  }

  /** \brief  Specify league size, request team size, specify vector length */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : TeamPolicyInternal(space_, league_size_, -1, vector_length_request) {}

  /** \brief  Specify league size, request team size and vector length */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     const Kokkos::AUTO_t& /* vector_length_request */
                     )
      : TeamPolicyInternal(space_, league_size_, -1, -1) {}

  /** \brief  Specify league size, specify team size, request vector length */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request, const Kokkos::AUTO_t&)
      : TeamPolicyInternal(space_, league_size_, team_size_request, -1) {}

  TeamPolicyInternal(int league_size_, int team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request) {}

  TeamPolicyInternal(int league_size_, const Kokkos::AUTO_t& team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request)

  {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(int league_size_, const Kokkos::AUTO_t& team_size_request,
                     const Kokkos::AUTO_t& vector_length_request)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(int league_size_, int team_size_request,
                     const Kokkos::AUTO_t& vector_length_request)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request) {}

  inline int chunk_size() const { return m_chunk_size; }

  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal& set_chunk_size(
      typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level,
                                              const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerTeamValue& per_team,
      const PerThreadValue& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  using member_type = Kokkos::Impl::CudaTeamMember;

 protected:
  template <class ClosureType, class FunctorType, class BlockSizeCallable>
  int internal_team_size_common(const FunctorType& f,
                                BlockSizeCallable&& block_size_callable) const {
    using closure_type = ClosureType;
    using Interface =
        typename Impl::DeduceFunctorPatternInterface<ClosureType>::type;
    using Analysis =
        Impl::FunctorAnalysis<Interface, typename ClosureType::Policy,
                              FunctorType>;

    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes();
    const int block_size = std::forward<BlockSizeCallable>(block_size_callable)(
        space().impl_internal_space_instance(), attr, f,
        (size_t)impl_vector_length(),
        (size_t)team_scratch_size(0) + 2 * sizeof(double),
        (size_t)thread_scratch_size(0) + sizeof(double) +
            ((Analysis::StaticValueSize != 0) ? 0 : Analysis::value_size(f)));
    KOKKOS_ASSERT(block_size > 0);

    // Currently we require Power-of-2 team size for reductions.
    int p2 = 1;
    while (p2 <= block_size) p2 *= 2;
    p2 /= 2;
    return p2 / impl_vector_length();
  }

  template <class ClosureType, class FunctorType>
  int internal_team_size_max(const FunctorType& f) const {
    return internal_team_size_common<ClosureType>(
        f,
        Kokkos::Impl::cuda_get_max_block_size<FunctorType,
                                              typename traits::launch_bounds>);
  }

  template <class ClosureType, class FunctorType>
  int internal_team_size_recommended(const FunctorType& f) const {
    return internal_team_size_common<ClosureType>(
        f,
        Kokkos::Impl::cuda_get_opt_block_size<FunctorType,
                                              typename traits::launch_bounds>);
  }
};

__device__ inline int64_t cuda_get_scratch_index(Cuda::size_type league_size,
                                                 int32_t* scratch_locks) {
  int64_t threadid = 0;
  __shared__ int64_t base_thread_id;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int64_t const wraparound_len = Kokkos::max(
        int64_t(1), Kokkos::min(int64_t(league_size),
                                (int64_t(g_device_cuda_lock_arrays.n)) /
                                    (blockDim.x * blockDim.y)));
    threadid = (blockIdx.x * blockDim.z + threadIdx.z) % wraparound_len;
    threadid *= blockDim.x * blockDim.y;
    int done = 0;
    while (!done) {
      done = (0 == atomicCAS(&scratch_locks[threadid], 0, 1));
      if (!done) {
        threadid += blockDim.x * blockDim.y;
        if (int64_t(threadid + blockDim.x * blockDim.y) >=
            wraparound_len * blockDim.x * blockDim.y)
          threadid = 0;
      }
    }
    base_thread_id = threadid;
  }
  __syncthreads();
  threadid = base_thread_id;
  return threadid;
}

__device__ inline void cuda_release_scratch_index(int32_t* scratch_locks,
                                                  int64_t threadid) {
  __syncthreads();
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    scratch_locks[threadid] = 0;
  }
}

template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Cuda> {
 public:
  using Policy = TeamPolicy<Properties...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

 public:
  using functor_type = FunctorType;
  using size_type    = Cuda::size_type;

 private:
  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y ==
  // blockDim.z == 1 shared memory utilization:
  //
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor;
  const Policy m_policy;
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;
  int m_shmem_begin;
  int m_shmem_size;
  void* m_scratch_ptr[2];
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const Member& member) const {
    m_functor(member);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const Member& member) const {
    m_functor(TagType(), member);
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  __device__ inline void operator()() const {
    // Iterate this block through the league
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = cuda_get_scratch_index(m_league_size, m_scratch_locks);
    }

    const int int_league_size = (int)m_league_size;
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<WorkTag>(typename Policy::member_type(
          kokkos_impl_cuda_shared_memory<void>(), m_shmem_begin, m_shmem_size,
          (void*)(((char*)m_scratch_ptr[1]) +
                  ptrdiff_t(threadid / (blockDim.x * blockDim.y)) *
                      m_scratch_size[1]),
          m_scratch_size[1], league_rank, m_league_size));
    }
    if (m_scratch_size[1] > 0) {
      cuda_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  inline void execute() const {
    const int64_t shmem_size_total = m_shmem_begin + m_shmem_size;
    dim3 grid(int(m_league_size), 1, 1);
    const dim3 block(int(m_vector_size), int(m_team_size), 1);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (Kokkos::Impl::CudaInternal::cuda_use_serial_execution()) {
      grid = dim3(1, 1, 1);
    }
#endif

    CudaParallelLaunch<ParallelFor, LaunchBounds>(
        *this, grid, block, shmem_size_total,
        m_policy.space()
	    .impl_internal_space_instance());  // copy to device and execute
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelFor,
                           LaunchBounds>::get_cuda_func_attributes();
    m_team_size =
        m_team_size >= 0
            ? m_team_size
            : Kokkos::Impl::cuda_get_opt_block_size<FunctorType, LaunchBounds>(
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
    m_scratch_locks =
        m_policy.space().impl_internal_space_instance()->m_scratch_locks;

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    m_scratch_ptr[0] = nullptr;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      auto scratch_ptr_id =
          m_policy.space()
              .impl_internal_space_instance()
              ->resize_team_scratch_space(
                  static_cast<std::int64_t>(m_scratch_size[1]) *
                  (std::min(
                      static_cast<std::int64_t>(Cuda::concurrency() /
                                                (m_team_size * m_vector_size)),
                      static_cast<std::int64_t>(m_league_size))));
      m_scratch_ptr[1]  = scratch_ptr_id.first;
      m_scratch_pool_id = scratch_ptr_id.second;
    }

    const int shmem_size_total = m_shmem_begin + m_shmem_size;
    if (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
        shmem_size_total) {
      printf(
          "%i %i\n",
          m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock,
          shmem_size_total);
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
    }

    if (int(m_team_size) >
        int(Kokkos::Impl::cuda_get_max_block_size<FunctorType, LaunchBounds>(
                m_policy.space().impl_internal_space_instance(), attr,
                arg_functor, arg_policy.impl_vector_length(),
                arg_policy.team_scratch_size(0),
                arg_policy.thread_scratch_size(0)) /
            arg_policy.impl_vector_length())) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< Cuda > requested too large team size."));
    }
  }

  ~ParallelFor() {
    if (m_scratch_pool_id >= 0) {
      m_policy.space()
          .impl_internal_space_instance()
          ->m_team_scratch_pool[m_scratch_pool_id] = 0;
    }
  }
};

template <class FunctorType, class ReducerType, class... Properties>
class ParallelReduce<FunctorType, Kokkos::TeamPolicy<Properties...>,
                     ReducerType, Kokkos::Cuda> {
 public:
  using Policy = TeamPolicy<Properties...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using Analysis =
      Kokkos::Impl::FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy,
                                    ReducerTypeFwd>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;
  using value_type     = typename Analysis::value_type;

 public:
  using functor_type = FunctorType;
  using size_type    = Cuda::size_type;
  using reducer_type = ReducerType;

  static constexpr bool UseShflReduction =
      (true && (Analysis::StaticValueSize != 0));

 private:
  struct ShflReductionTag {};
  struct SHMEMReductionTag {};

  // Algorithmic constraints: blockDim.y is a power of two AND blockDim.y ==
  // blockDim.z == 1 shared memory utilization:
  //
  //  [ global reduce space ]
  //  [ team   reduce space ]
  //  [ team   shared space ]
  //

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type* m_unified_space;
  size_type m_team_begin;
  size_type m_shmem_begin;
  size_type m_shmem_size;
  void* m_scratch_ptr[2];
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const Member& member, reference_type update) const {
    m_functor(member, update);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const Member& member, reference_type update) const {
    m_functor(TagType(), member, update);
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  __device__ inline void operator()() const {
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = cuda_get_scratch_index(m_league_size, m_scratch_locks);
    }

    using ReductionTag = std::conditional_t<UseShflReduction, ShflReductionTag,
                                            SHMEMReductionTag>;
    run(ReductionTag{}, threadid);
    if (m_scratch_size[1] > 0) {
      cuda_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  __device__ inline void run(SHMEMReductionTag&, const int& threadid) const {
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    const integral_nonzero_constant<size_type, Analysis::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(Analysis::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    reference_type value =
        final_reducer.init(kokkos_impl_cuda_shared_memory<size_type>() +
                           threadIdx.y * word_count.value);

    // Iterate this block through the league
    const int int_league_size = (int)m_league_size;
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<WorkTag>(
          Member(kokkos_impl_cuda_shared_memory<char>() + m_team_begin,
                 m_shmem_begin, m_shmem_size,
                 (void*)(((char*)m_scratch_ptr[1]) +
                         ptrdiff_t(threadid / (blockDim.x * blockDim.y)) *
                             m_scratch_size[1]),
                 m_scratch_size[1], league_rank, m_league_size),
          value);
    }

    // Reduce with final value at blockDim.y - 1 location.
    bool zero_length        = m_league_size == 0;
    bool do_final_reduction = true;
    if (!zero_length)
      do_final_reduction = cuda_single_inter_block_reduce_scan<false>(
          final_reducer, blockIdx.x, gridDim.x,
          kokkos_impl_cuda_shared_memory<size_type>(), m_scratch_space,
          m_scratch_flags);

    if (do_final_reduction) {
      // This is the final block with the final result at the final threads'
      // location

      size_type* const shared = kokkos_impl_cuda_shared_memory<size_type>() +
                                (blockDim.y - 1) * word_count.value;
      size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        final_reducer.final(reinterpret_cast<value_type*>(shared));
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  __device__ inline void run(ShflReductionTag, const int& threadid) const {
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    value_type value;
    final_reducer.init(&value);

    // Iterate this block through the league
    const int int_league_size = (int)m_league_size;
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<WorkTag>(
          Member(kokkos_impl_cuda_shared_memory<char>() + m_team_begin,
                 m_shmem_begin, m_shmem_size,
                 (void*)(((char*)m_scratch_ptr[1]) +
                         ptrdiff_t(threadid / (blockDim.x * blockDim.y)) *
                             m_scratch_size[1]),
                 m_scratch_size[1], league_rank, m_league_size),
          value);
    }

    pointer_type const result =
        m_result_ptr_device_accessible
            ? m_result_ptr
            : (pointer_type)(m_unified_space ? m_unified_space
                                             : m_scratch_space);

    value_type init;
    final_reducer.init(&init);

    if (int_league_size == 0) {
      final_reducer.final(&value);
      *result = value;
    } else if (Impl::cuda_inter_block_reduction(value, init, final_reducer,
                                                m_scratch_space, result,
                                                m_scratch_flags, blockDim.y)) {
      const unsigned id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        final_reducer.final(&value);
        *result = value;
      }
    }
  }

  inline void execute() {
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));

    const bool is_empty_range  = m_league_size == 0 || m_team_size == 0;
    const bool need_device_set = Analysis::has_init_member_function ||
                                 Analysis::has_final_member_function ||
                                 !m_result_ptr_host_accessible ||
#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
                                 Policy::is_graph_kernel::value ||
#endif
                                 !std::is_same<ReducerType, InvalidType>::value;
    if (!is_empty_range || need_device_set) {
      const int block_count = std::max(
          1u, UseShflReduction ? std::min(m_league_size, size_type(1024 * 32))
                               : std::min(int(m_league_size), m_team_size));

      m_scratch_space = cuda_internal_scratch_space(
          m_policy.space(), Analysis::value_size(ReducerConditional::select(
                                m_functor, m_reducer)) *
                                block_count);
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type));
      m_unified_space = cuda_internal_scratch_unified(
          m_policy.space(), Analysis::value_size(ReducerConditional::select(
                                m_functor, m_reducer)));

      dim3 block(m_vector_size, m_team_size, 1);
      dim3 grid(block_count, 1, 1);
      const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

      if (is_empty_range
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
          || Kokkos::Impl::CudaInternal::cuda_use_serial_execution()
#endif
      ) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }

      CudaParallelLaunch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem_size_total,
          m_policy.space()
	      .impl_internal_space_instance());  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().fence(
            "Kokkos::Impl::ParallelReduce<Cuda, TeamPolicy>::execute: Result "
            "Not Device Accessible");

        if (m_result_ptr) {
          if (m_unified_space) {
            const int count = Analysis::value_count(
                ReducerConditional::select(m_functor, m_reducer));
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = Analysis::value_size(
                ReducerConditional::select(m_functor, m_reducer));
            DeepCopy<HostSpace, CudaSpace>(m_result_ptr, m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        final_reducer.init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(
      const FunctorType& arg_functor, const Policy& arg_policy,
      const ViewType& arg_result,
      std::enable_if_t<Kokkos::is_view<ViewType>::value, void*> = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::CudaSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_unified_space(nullptr),
        m_team_begin(0),
        m_shmem_begin(0),
        m_shmem_size(0),
        m_scratch_ptr{nullptr, nullptr},
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelReduce,
                           LaunchBounds>::get_cuda_func_attributes();
    m_team_size =
        m_team_size >= 0
            ? m_team_size
            : Kokkos::Impl::cuda_get_opt_block_size<FunctorType, LaunchBounds>(
                  m_policy.space().impl_internal_space_instance(), attr,
                  m_functor, m_vector_size, m_policy.team_scratch_size(0),
                  m_policy.thread_scratch_size(0)) /
                  m_vector_size;

    m_team_begin =
        UseShflReduction
            ? 0
            : cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                        WorkTag>(arg_functor,
                                                                 m_team_size);
    m_shmem_begin = sizeof(double) * (m_team_size + 2);
    m_shmem_size =
        m_policy.scratch_size(0, m_team_size) +
        FunctorTeamShmemSize<FunctorType>::value(arg_functor, m_team_size);
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks =
        m_policy.space().impl_internal_space_instance()->m_scratch_locks;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      auto scratch_ptr_id =
          m_policy.space()
              .impl_internal_space_instance()
              ->resize_team_scratch_space(
                  static_cast<std::int64_t>(m_scratch_size[1]) *
                  (std::min(
                      static_cast<std::int64_t>(Cuda::concurrency() /
                                                (m_team_size * m_vector_size)),
                      static_cast<std::int64_t>(m_league_size))));
      m_scratch_ptr[1]  = scratch_ptr_id.first;
      m_scratch_pool_id = scratch_ptr_id.second;
    }

    // The global parallel_reduce does not support vector_length other than 1 at
    // the moment
    if ((arg_policy.impl_vector_length() > 1) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a vector length of "
          "greater than 1 is not currently supported for CUDA for dynamic "
          "sized reduction types.");

    if ((m_team_size < 32) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller "
          "than 32 is not currently supported with CUDA for dynamic sized "
          "reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

    if (!Kokkos::Impl::is_integral_power_of_two(m_team_size) &&
        !UseShflReduction) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    if (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
        shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too much "
                      "L0 scratch memory"));
    }

    if (int(m_team_size) >
        arg_policy.team_size_max(m_functor, m_reducer, ParallelReduceTag())) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too "
                      "large team size."));
    }
  }

  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::CudaSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_unified_space(nullptr),
        m_team_begin(0),
        m_shmem_begin(0),
        m_shmem_size(0),
        m_scratch_ptr{nullptr, nullptr},
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelReduce,
                           LaunchBounds>::get_cuda_func_attributes();

    // Valid team size not provided, deduce team size
    m_team_size =
        m_team_size >= 0
            ? m_team_size
            : Kokkos::Impl::cuda_get_opt_block_size<FunctorType, LaunchBounds>(
                  m_policy.space().impl_internal_space_instance(), attr,
                  m_functor, m_vector_size, m_policy.team_scratch_size(0),
                  m_policy.thread_scratch_size(0)) /
                  m_vector_size;

    m_team_begin =
        UseShflReduction
            ? 0
            : cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                        WorkTag>(arg_functor,
                                                                 m_team_size);
    m_shmem_begin = sizeof(double) * (m_team_size + 2);
    m_shmem_size =
        m_policy.scratch_size(0, m_team_size) +
        FunctorTeamShmemSize<FunctorType>::value(arg_functor, m_team_size);
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks =
        m_policy.space().impl_internal_space_instance()->m_scratch_locks;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      auto scratch_ptr_id =
          m_policy.space()
              .impl_internal_space_instance()
              ->resize_team_scratch_space(
                  static_cast<std::int64_t>(m_scratch_size[1]) *
                  (std::min(
                      static_cast<std::int64_t>(Cuda::concurrency() /
                                                (m_team_size * m_vector_size)),
                      static_cast<std::int64_t>(m_league_size))));
      m_scratch_ptr[1]  = scratch_ptr_id.first;
      m_scratch_pool_id = scratch_ptr_id.second;
    }

    // The global parallel_reduce does not support vector_length other than 1 at
    // the moment
    if ((arg_policy.impl_vector_length() > 1) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a vector length of "
          "greater than 1 is not currently supported for CUDA for dynamic "
          "sized reduction types.");

    if ((m_team_size < 32) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller "
          "than 32 is not currently supported with CUDA for dynamic sized "
          "reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

    if ((!Kokkos::Impl::is_integral_power_of_two(m_team_size) &&
         !UseShflReduction) ||
        m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
            shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    size_type team_size_max =
        Kokkos::Impl::cuda_get_max_block_size<FunctorType, LaunchBounds>(
            m_policy.space().impl_internal_space_instance(), attr, m_functor,
            m_vector_size, m_policy.team_scratch_size(0),
            m_policy.thread_scratch_size(0)) /
        m_vector_size;

    if ((int)m_team_size > (int)team_size_max) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too "
                      "large team size."));
    }
  }

  ~ParallelReduce() {
    if (m_scratch_pool_id >= 0) {
      m_policy.space()
          .impl_internal_space_instance()
          ->m_team_scratch_pool[m_scratch_pool_id] = 0;
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos
#endif /* defined(KOKKOS_ENABLE_CUDA) */
#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */
