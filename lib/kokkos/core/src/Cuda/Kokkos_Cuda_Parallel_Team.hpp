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
#include <Cuda/Kokkos_Cuda_Team.hpp>
#include <Kokkos_MinMax.hpp>
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
            get_cuda_func_attributes(space().cuda_device());
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
                              TeamPolicyInternal, FunctorType, void>;
    using closure_type = Impl::ParallelReduce<
        CombinedFunctorReducer<FunctorType,
                               typename functor_analysis_type::Reducer>,
        TeamPolicy<Properties...>, Kokkos::Cuda>;
    return internal_team_size_max<closure_type>(f);
  }

  template <class FunctorType, class ReducerType>
  inline int team_size_max(const FunctorType& f, const ReducerType& /*r*/,
                           const ParallelReduceTag&) const {
    using closure_type =
        Impl::ParallelReduce<CombinedFunctorReducer<FunctorType, ReducerType>,
                             TeamPolicy<Properties...>, Kokkos::Cuda>;
    return internal_team_size_max<closure_type>(f);
  }

  template <class FunctorType>
  int team_size_recommended(const FunctorType& f, const ParallelForTag&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...>>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes(space().cuda_device());
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
                              TeamPolicyInternal, FunctorType, void>;
    using closure_type = Impl::ParallelReduce<
        CombinedFunctorReducer<FunctorType,
                               typename functor_analysis_type::Reducer>,
        TeamPolicy<Properties...>, Kokkos::Cuda>;
    return internal_team_size_recommended<closure_type>(f);
  }

  template <class FunctorType, class ReducerType>
  int team_size_recommended(const FunctorType& f, const ReducerType&,
                            const ParallelReduceTag&) const {
    using closure_type =
        Impl::ParallelReduce<CombinedFunctorReducer<FunctorType, ReducerType>,
                             TeamPolicy<Properties...>, Kokkos::Cuda>;
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
    // Cuda Teams use (team_size + 2)*sizeof(double) shared memory for team
    // reductions. They also use one int64_t in static shared memory for a
    // shared ID. Furthermore, they use additional scratch memory in some
    // reduction scenarios, which depend on the size of the value_type and is
    // NOT captured here.
    constexpr size_t max_possible_team_size = 1024;
    constexpr size_t max_reserved_shared_mem_per_team =
        (max_possible_team_size + 2) * sizeof(double) + sizeof(int64_t);
    // arbitrarily setting level 1 scratch limit to 20MB, for a
    // Volta V100 that would give us about 3.2GB for 2 teams per SM
    constexpr size_t max_l1_scratch_size = 20 * 1024 * 1024;

    size_t max_shmem = Cuda().cuda_device_prop().sharedMemPerBlock;
    return (level == 0 ? max_shmem - max_reserved_shared_mem_per_team
                       : max_l1_scratch_size);
  }

  //----------------------------------------

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
    const int maxGridSizeX = m_space.cuda_device_prop().maxGridSize[0];
    if (league_size_ >= maxGridSizeX)
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
                              FunctorType, void>;

    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes(space().cuda_device());
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
                                                 int32_t* scratch_locks,
                                                 size_t num_scratch_locks) {
  int64_t threadid = 0;
  __shared__ int64_t base_thread_id;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int64_t const wraparound_len = Kokkos::max(
        int64_t(1),
        Kokkos::min(int64_t(league_size),
                    int64_t(num_scratch_locks) / (blockDim.x * blockDim.y)));
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
  size_t m_num_scratch_locks;

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
      threadid = cuda_get_scratch_index(m_league_size, m_scratch_locks,
                                        m_num_scratch_locks);
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
    auto internal_space_instance =
        m_policy.space().impl_internal_space_instance();
    m_team_size = m_team_size >= 0 ? m_team_size
                                   : arg_policy.team_size_recommended(
                                         arg_functor, ParallelForTag());

    m_shmem_begin = (sizeof(double) * (m_team_size + 2));
    m_shmem_size =
        (m_policy.scratch_size(0, m_team_size) +
         FunctorTeamShmemSize<FunctorType>::value(m_functor, m_team_size));
    m_scratch_size[0]   = m_policy.scratch_size(0, m_team_size);
    m_scratch_size[1]   = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks     = internal_space_instance->m_scratch_locks;
    m_num_scratch_locks = internal_space_instance->m_num_scratch_locks;

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    m_scratch_ptr[0] = nullptr;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      m_scratch_pool_id = internal_space_instance->acquire_team_scratch_space();
      m_scratch_ptr[1]  = internal_space_instance->resize_team_scratch_space(
          m_scratch_pool_id,
          static_cast<std::int64_t>(m_scratch_size[1]) *
              (std::min(
                  static_cast<std::int64_t>(Cuda().concurrency() /
                                            (m_team_size * m_vector_size)),
                  static_cast<std::int64_t>(m_league_size))));
    }

    const int maxShmemPerBlock =
        m_policy.space().cuda_device_prop().sharedMemPerBlock;
    const int shmem_size_total = m_shmem_begin + m_shmem_size;
    if (maxShmemPerBlock < shmem_size_total) {
      printf("%i %i\n", maxShmemPerBlock, shmem_size_total);
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< Cuda > insufficient shared memory"));
    }

    if (m_team_size > arg_policy.team_size_max(arg_functor, ParallelForTag())) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< Cuda > requested too large team size."));
    }
  }

  ~ParallelFor() {
    if (m_scratch_pool_id >= 0) {
      m_policy.space()
          .impl_internal_space_instance()
          ->release_team_scratch_space(m_scratch_pool_id);
    }
  }
};

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, Kokkos::Cuda> {
 public:
  using Policy      = TeamPolicy<Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

 public:
  using functor_type = FunctorType;
  // Conditionally set word_size_type to int16_t or int8_t if value_type is
  // smaller than int32_t (Kokkos::Cuda::size_type)
  // word_size_type is used to determine the word count, shared memory buffer
  // size, and global memory buffer size before the reduction is performed.
  // Within the reduction, the word count is recomputed based on word_size_type
  // and when calculating indexes into the shared/global memory buffers for
  // performing the reduction, word_size_type is used again.
  // For scalars > 4 bytes in size, indexing into shared/global memory relies
  // on the block and grid dimensions to ensure that we index at the correct
  // offset rather than at every 4 byte word; such that, when the join is
  // performed, we have the correct data that was copied over in chunks of 4
  // bytes.
  using word_size_type = std::conditional_t<
      sizeof(value_type) < sizeof(Kokkos::Cuda::size_type),
      std::conditional_t<sizeof(value_type) == 2, int16_t, int8_t>,
      Kokkos::Cuda::size_type>;
  using size_type    = Cuda::size_type;
  using reducer_type = ReducerType;

  static constexpr bool UseShflReduction =
      (true && (ReducerType::static_value_size() != 0));

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

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  word_size_type* m_scratch_space;
  // m_scratch_flags must be of type Cuda::size_type due to use of atomics
  // for tracking metadata in Kokkos_Cuda_ReduceScan.hpp
  Cuda::size_type* m_scratch_flags;
  word_size_type* m_unified_space;
  size_type m_team_begin;
  size_type m_shmem_begin;
  size_type m_shmem_size;
  void* m_scratch_ptr[2];
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  size_t m_num_scratch_locks;
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;

  template <class TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const Member& member, reference_type update) const {
    m_functor_reducer.get_functor()(member, update);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const Member& member, reference_type update) const {
    m_functor_reducer.get_functor()(TagType(), member, update);
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  __device__ inline void operator()() const {
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = cuda_get_scratch_index(m_league_size, m_scratch_locks,
                                        m_num_scratch_locks);
    }

    using ReductionTag = std::conditional_t<UseShflReduction, ShflReductionTag,
                                            SHMEMReductionTag>;
    run(ReductionTag{}, threadid);
    if (m_scratch_size[1] > 0) {
      cuda_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  __device__ inline void run(SHMEMReductionTag&, const int& threadid) const {
    const integral_nonzero_constant<word_size_type,
                                    ReducerType::static_value_size() /
                                        sizeof(word_size_type)>
        word_count(m_functor_reducer.get_reducer().value_size() /
                   sizeof(word_size_type));

    reference_type value = m_functor_reducer.get_reducer().init(
        kokkos_impl_cuda_shared_memory<word_size_type>() +
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
          m_functor_reducer.get_reducer(), blockIdx.x, gridDim.x,
          kokkos_impl_cuda_shared_memory<word_size_type>(), m_scratch_space,
          m_scratch_flags);

    if (do_final_reduction) {
      // This is the final block with the final result at the final threads'
      // location

      word_size_type* const shared =
          kokkos_impl_cuda_shared_memory<word_size_type>() +
          (blockDim.y - 1) * word_count.value;
      size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<word_size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        m_functor_reducer.get_reducer().final(
            reinterpret_cast<value_type*>(shared));
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      } else {
        // In the above call to final(), shared might have been updated by a
        // single thread within a warp without synchronization. Synchronize
        // threads within warp to avoid potential race condition.
        __syncwarp(0xffffffff);
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  __device__ inline void run(ShflReductionTag, const int& threadid) const {
    value_type value;
    m_functor_reducer.get_reducer().init(&value);

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
    m_functor_reducer.get_reducer().init(&init);

    if (int_league_size == 0) {
      m_functor_reducer.get_reducer().final(&value);
      *result = value;
    } else if (Impl::cuda_inter_block_reduction(
                   value, init, m_functor_reducer.get_reducer(),
                   reinterpret_cast<pointer_type>(m_scratch_space), result,
                   m_scratch_flags, blockDim.y)) {
      const unsigned id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        m_functor_reducer.get_reducer().final(&value);
        *result = value;
      }
    }
  }

  inline void execute() {
    const bool is_empty_range  = m_league_size == 0 || m_team_size == 0;
    const bool need_device_set = ReducerType::has_init_member_function() ||
                                 ReducerType::has_final_member_function() ||
                                 !m_result_ptr_host_accessible ||
                                 Policy::is_graph_kernel::value ||
                                 !std::is_same<ReducerType, InvalidType>::value;
    if (!is_empty_range || need_device_set) {
      const int block_count = std::max(
          1u, UseShflReduction ? std::min(m_league_size, size_type(1024 * 32))
                               : std::min(int(m_league_size), m_team_size));

      m_scratch_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_space(
              m_policy.space(),
              m_functor_reducer.get_reducer().value_size() * block_count));
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type));
      m_unified_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_unified(
              m_policy.space(), m_functor_reducer.get_reducer().value_size()));

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
            const int count = m_functor_reducer.get_reducer().value_count();
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = m_functor_reducer.get_reducer().value_size();
            DeepCopy<HostSpace, CudaSpace, Cuda>(m_policy.space(), m_result_ptr,
                                                 m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        m_functor_reducer.get_reducer().init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
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
    auto internal_space_instance =
        m_policy.space().impl_internal_space_instance();
    m_team_size = m_team_size >= 0 ? m_team_size
                                   : arg_policy.team_size_recommended(
                                         arg_functor_reducer.get_functor(),
                                         arg_functor_reducer.get_reducer(),
                                         ParallelReduceTag());

    m_team_begin =
        UseShflReduction
            ? 0
            : cuda_single_inter_block_reduce_scan_shmem<false, WorkTag,
                                                        value_type>(
                  arg_functor_reducer.get_functor(), m_team_size);
    m_shmem_begin = sizeof(double) * (m_team_size + 2);
    m_shmem_size  = m_policy.scratch_size(0, m_team_size) +
                   FunctorTeamShmemSize<FunctorType>::value(
                       arg_functor_reducer.get_functor(), m_team_size);
    m_scratch_size[0]   = m_shmem_size;
    m_scratch_size[1]   = m_policy.scratch_size(1, m_team_size);
    m_scratch_locks     = internal_space_instance->m_scratch_locks;
    m_num_scratch_locks = internal_space_instance->m_num_scratch_locks;
    if (m_team_size <= 0) {
      m_scratch_ptr[1] = nullptr;
    } else {
      m_scratch_pool_id = internal_space_instance->acquire_team_scratch_space();
      m_scratch_ptr[1]  = internal_space_instance->resize_team_scratch_space(
          m_scratch_pool_id,
          static_cast<std::int64_t>(m_scratch_size[1]) *
              (std::min(
                  static_cast<std::int64_t>(Cuda().concurrency() /
                                            (m_team_size * m_vector_size)),
                  static_cast<std::int64_t>(m_league_size))));
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

    const int maxShmemPerBlock =
        m_policy.space().cuda_device_prop().sharedMemPerBlock;
    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

    if (!Kokkos::Impl::is_integral_power_of_two(m_team_size) &&
        !UseShflReduction) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > bad team size"));
    }

    if (maxShmemPerBlock < shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too much "
                      "L0 scratch memory"));
    }

    if (int(m_team_size) >
        arg_policy.team_size_max(m_functor_reducer.get_functor(),
                                 m_functor_reducer.get_reducer(),
                                 ParallelReduceTag())) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< Cuda > requested too "
                      "large team size."));
    }
  }

  ~ParallelReduce() {
    if (m_scratch_pool_id >= 0) {
      m_policy.space()
          .impl_internal_space_instance()
          ->release_team_scratch_space(m_scratch_pool_id);
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos
#endif /* defined(KOKKOS_ENABLE_CUDA) */
#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */
