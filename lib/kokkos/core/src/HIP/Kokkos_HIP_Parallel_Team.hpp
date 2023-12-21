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

#ifndef KOKKO_HIP_PARALLEL_TEAM_HPP
#define KOKKO_HIP_PARALLEL_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_Team.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <Kokkos_MinMaxClamp.hpp>

namespace Kokkos {
namespace Impl {

template <typename... Properties>
class TeamPolicyInternal<HIP, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  template <typename ExecSpace, typename... OtherProperties>
  friend class TeamPolicyInternal;

 private:
  typename traits::execution_space m_space;
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];
  int m_chunk_size;
  bool m_tune_team_size;
  bool m_tune_vector_length;

 public:
  using execution_space = HIP;

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
    m_tune_team_size         = p.m_tune_team_size;
    m_tune_vector_length     = p.m_tune_vector_length;
  }

  template <typename FunctorType>
  int team_size_max(FunctorType const& f, ParallelForTag const&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...>>;

    return internal_team_size_common<BlockType::Max, closure_type, void>(f);
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
        TeamPolicy<Properties...>, Kokkos::HIP>;
    return internal_team_size_common<
        BlockType::Max, closure_type,
        typename functor_analysis_type::value_type>(f);
  }

  template <typename FunctorType, typename ReducerType>
  inline int team_size_max(const FunctorType& f, const ReducerType&,
                           const ParallelReduceTag&) const {
    using closure_type =
        Impl::ParallelReduce<CombinedFunctorReducer<FunctorType, ReducerType>,
                             TeamPolicy<Properties...>, Kokkos::HIP>;
    return internal_team_size_common<BlockType::Max, closure_type,
                                     typename ReducerType::value_type>(f);
  }

  template <typename FunctorType>
  int team_size_recommended(FunctorType const& f, ParallelForTag const&) const {
    using closure_type =
        Impl::ParallelFor<FunctorType, TeamPolicy<Properties...>>;

    return internal_team_size_common<BlockType::Preferred, closure_type, void>(
        f);
  }

  template <typename FunctorType>
  inline int team_size_recommended(FunctorType const& f,
                                   ParallelReduceTag const&) const {
    using functor_analysis_type =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                              TeamPolicyInternal, FunctorType, void>;
    using closure_type = Impl::ParallelReduce<
        CombinedFunctorReducer<FunctorType,
                               typename functor_analysis_type::Reducer>,
        TeamPolicy<Properties...>, Kokkos::HIP>;
    return internal_team_size_common<
        BlockType::Preferred, closure_type,
        typename functor_analysis_type::value_type>(f);
  }

  template <typename FunctorType, typename ReducerType>
  int team_size_recommended(FunctorType const& f, ReducerType const&,
                            ParallelReduceTag const&) const {
    using closure_type =
        Impl::ParallelReduce<CombinedFunctorReducer<FunctorType, ReducerType>,
                             TeamPolicy<Properties...>, Kokkos::HIP>;
    return internal_team_size_common<BlockType::Preferred, closure_type,
                                     typename ReducerType::value_type>(f);
  }

  inline bool impl_auto_vector_length() const { return m_tune_vector_length; }
  inline bool impl_auto_team_size() const { return m_tune_team_size; }
  static int vector_length_max() { return HIPTraits::WarpSize; }

  static int verify_requested_vector_length(int requested_vector_length) {
    int test_vector_length =
        std::min(requested_vector_length, vector_length_max());

    // Allow only power-of-two vector_length
    if (!(is_integral_power_of_two(test_vector_length))) {
      int test_pow2           = 1;
      constexpr int warp_size = HIPTraits::WarpSize;
      while (test_pow2 < warp_size) {
        test_pow2 <<= 1;
        if (test_pow2 > test_vector_length) {
          break;
        }
      }
      test_vector_length = test_pow2 >> 1;
    }

    return test_vector_length;
  }

  inline static int scratch_size_max(int level) {
    // HIP Teams use (team_size + 2)*sizeof(double) shared memory for team
    // reductions. They also use one int64_t in static shared memory for a
    // shared ID. Furthermore, they use additional scratch memory in some
    // reduction scenarios, which depend on the size of the value_type and is
    // NOT captured here
    constexpr size_t max_possible_team_size = 1024;
    constexpr size_t max_reserved_shared_mem_per_team =
        (max_possible_team_size + 2) * sizeof(double) + sizeof(int64_t);
    // arbitrarily setting level 1 scratch limit to 20MB, for a
    // MI250 that would give us about 4.4GB for 2 teams per CU
    constexpr size_t max_l1_scratch_size = 20 * 1024 * 1024;

    size_t max_shmem = HIP().hip_device_prop().sharedMemPerBlock;
    return (level == 0 ? max_shmem - max_reserved_shared_mem_per_team
                       : max_l1_scratch_size);
  }

  inline void impl_set_vector_length(size_t size) { m_vector_length = size; }
  inline void impl_set_team_size(size_t size) { m_team_size = size; }
  int impl_vector_length() const { return m_vector_length; }

  int team_size() const { return m_team_size; }

  int league_size() const { return m_league_size; }

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

  typename traits::execution_space space() const { return m_space; }

  TeamPolicyInternal()
      : m_space(typename traits::execution_space()),
        m_league_size(0),
        m_team_size(-1),
        m_vector_length(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(HIPTraits::WarpSize),
        m_tune_team_size(false),
        m_tune_vector_length(false) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request, int vector_length_request = 1)
      : m_space(space_),
        m_league_size(league_size_),
        m_team_size(team_size_request),
        m_vector_length(
            (vector_length_request > 0)
                ? verify_requested_vector_length(vector_length_request)
                : (verify_requested_vector_length(1))),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(HIPTraits::WarpSize),
        m_tune_team_size(bool(team_size_request <= 0)),
        m_tune_vector_length(bool(vector_length_request <= 0)) {
    // Make sure league size is permissible
    if (league_size_ >= static_cast<int>(hip_internal_maximum_grid_count()[0]))
      Impl::throw_runtime_exception(
          "Requested too large league_size for TeamPolicy on HIP execution "
          "space.");

    // Make sure total block size is permissible
    if (m_team_size * m_vector_length > HIPTraits::MaxThreadsPerBlock) {
      Impl::throw_runtime_exception(
          std::string("Kokkos::TeamPolicy< HIP > the team size is too large. "
                      "Team size x vector length must be smaller than 1024."));
    }
  }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : TeamPolicyInternal(space_, league_size_, -1, vector_length_request) {}
  // FLAG
  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */
                     )
      : TeamPolicyInternal(space_, league_size_, team_size_request, -1)

  {}

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(const execution_space space_, int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(space_, league_size_, -1, -1)

  {}

  TeamPolicyInternal(int league_size_, int team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, vector_length_request) {}

  TeamPolicyInternal(int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(), league_size_, -1,
                           vector_length_request) {}

  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(int league_size_, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(typename traits::execution_space(), league_size_,
                           team_size_request, -1)

  {}

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(int league_size_,
                     const Kokkos::AUTO_t& /* team_size_request */,
                     const Kokkos::AUTO_t& /* vector_length_request */

                     )
      : TeamPolicyInternal(typename traits::execution_space(), league_size_, -1,
                           -1) {}

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
  template <BlockType BlockSize, class ClosureType, class ValueType,
            class FunctorType>
  int internal_team_size_common(FunctorType const& f) const {
    const unsigned shmem_block = team_scratch_size(0) + 2 * sizeof(double);
    unsigned shmem_thread      = thread_scratch_size(0) + sizeof(double);
    using Tag = typename PatternTagFromImplSpecialization<ClosureType>::type;
    if constexpr (std::is_same_v<Tag, ParallelReduceTag>) {
      using Interface =
          typename Impl::DeduceFunctorPatternInterface<ClosureType>::type;
      using Analysis =
          Impl::FunctorAnalysis<Interface, typename ClosureType::Policy,
                                FunctorType, ValueType>;
      shmem_thread +=
          ((Analysis::StaticValueSize != 0) ? 0 : Analysis::value_size(f));
    }
    const int vector_length = impl_vector_length();

    const auto functor = [&f, shmem_block, shmem_thread, vector_length](
                             const hipFuncAttributes& attr, int block_size) {
      int functor_shmem =
          ::Kokkos::Impl::FunctorTeamShmemSize<FunctorType>::value(
              f, block_size / vector_length);
      return shmem_block + shmem_thread * (block_size / vector_length) +
             functor_shmem + attr.sharedSizeBytes;
    };
    int block_size;
    if constexpr (BlockSize == BlockType::Max) {
      block_size = hip_get_max_team_blocksize<ClosureType,
                                              typename traits::launch_bounds>(
          space().impl_internal_space_instance(), functor);
    } else {
      block_size =
          hip_get_preferred_team_blocksize<ClosureType,
                                           typename traits::launch_bounds>(
              space().impl_internal_space_instance(), functor);
    }

    if (block_size == 0) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor/Reduce< HIP > could not find a valid "
          "team size."));
    }
    if constexpr (std::is_same_v<Tag, ParallelForTag>) {
      return block_size / impl_vector_length();
    } else {
      // Currently we require Power-of-2 team size for reductions.
      int p2 = 1;
      while (p2 <= block_size) p2 *= 2;
      p2 /= 2;
      return p2 / impl_vector_length();
    }
  }
};

__device__ inline int64_t hip_get_scratch_index(HIP::size_type league_size,
                                                int32_t* scratch_locks,
                                                size_t num_scratch_locks) {
  int64_t threadid = 0;
  __shared__ int64_t base_thread_id;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int64_t const wraparound_len =
        Kokkos::min(int64_t(league_size),
                    int64_t(num_scratch_locks) / (blockDim.x * blockDim.y));
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

__device__ inline void hip_release_scratch_index(int32_t* scratch_locks,
                                                 int64_t threadid) {
  __syncthreads();
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    scratch_locks[threadid] = 0;
  }
}

template <typename FunctorType, typename... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>, HIP> {
 public:
  using Policy       = TeamPolicy<Properties...>;
  using functor_type = FunctorType;
  using size_type    = HIP::size_type;

 private:
  using member_type   = typename Policy::member_type;
  using work_tag      = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  // Algorithmic constraints: blockDim.y is a power of two AND
  // blockDim.y  == blockDim.z == 1 shared memory utilization:
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
  size_t m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  size_t m_num_scratch_locks;

  template <typename TagType>
  __device__ inline std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const member_type& member) const {
    m_functor(member);
  }

  template <typename TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const member_type& member) const {
    m_functor(TagType(), member);
  }

 public:
  ParallelFor()                   = delete;
  ParallelFor(ParallelFor const&) = default;
  ParallelFor& operator=(ParallelFor const&) = delete;

  __device__ inline void operator()() const {
    // Iterate this block through the league
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = hip_get_scratch_index(m_league_size, m_scratch_locks,
                                       m_num_scratch_locks);
    }

    int const int_league_size = static_cast<int>(m_league_size);
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<work_tag>(typename Policy::member_type(
          kokkos_impl_hip_shared_memory<void>(), m_shmem_begin, m_shmem_size,
          static_cast<void*>(static_cast<char*>(m_scratch_ptr[1]) +
                             ptrdiff_t(threadid / (blockDim.x * blockDim.y)) *
                                 m_scratch_size[1]),
          m_scratch_size[1], league_rank, m_league_size));
    }
    if (m_scratch_size[1] > 0) {
      hip_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  inline void execute() const {
    int64_t const shmem_size_total = m_shmem_begin + m_shmem_size;
    dim3 const grid(static_cast<int>(m_league_size), 1, 1);
    dim3 const block(static_cast<int>(m_vector_size),
                     static_cast<int>(m_team_size), 1);

    using closure_type =
        ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>, HIP>;
    Impl::hip_parallel_launch<closure_type, launch_bounds>(
        *this, grid, block, shmem_size_total,
        m_policy.space().impl_internal_space_instance(),
        true);  // copy to device and execute
  }

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
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
                  static_cast<std::int64_t>(HIP().concurrency() /
                                            (m_team_size * m_vector_size)),
                  static_cast<std::int64_t>(m_league_size))));
    }

    int const shmem_size_total = m_shmem_begin + m_shmem_size;
    if (internal_space_instance->m_maxShmemPerBlock < shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > insufficient shared memory"));
    }

    size_t max_size = arg_policy.team_size_max(arg_functor, ParallelForTag());
    if (static_cast<int>(m_team_size) > static_cast<int>(max_size)) {
      Kokkos::Impl::throw_runtime_exception(std::string(
          "Kokkos::Impl::ParallelFor< HIP > requested too large team size."));
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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, HIP> {
 public:
  using Policy      = TeamPolicyInternal<HIP, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using member_type   = typename Policy::member_type;
  using work_tag      = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

 public:
  using functor_type = FunctorType;
  using size_type    = HIP::size_type;

  // static int constexpr UseShflReduction = false;
  // FIXME_HIP This should be disabled unconditionally for best performance, but
  // it currently causes tests to fail.
  static constexpr int UseShflReduction =
      (ReducerType::static_value_size() != 0);

 private:
  struct ShflReductionTag {};
  struct SHMEMReductionTag {};

  // Algorithmic constraints: blockDim.y is a power of two AND
  // blockDim.y == blockDim.z == 1 shared memory utilization:
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
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
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
      member_type const& member, reference_type update) const {
    m_functor_reducer.get_functor()(member, update);
  }

  template <class TagType>
  __device__ inline std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      member_type const& member, reference_type update) const {
    m_functor_reducer.get_functor()(TagType(), member, update);
  }

  __device__ inline void iterate_through_league(int const threadid,
                                                reference_type value) const {
    int const int_league_size = static_cast<int>(m_league_size);
    for (int league_rank = blockIdx.x; league_rank < int_league_size;
         league_rank += gridDim.x) {
      this->template exec_team<work_tag>(
          member_type(
              kokkos_impl_hip_shared_memory<char>() + m_team_begin,
              m_shmem_begin, m_shmem_size,
              reinterpret_cast<void*>(
                  reinterpret_cast<char*>(m_scratch_ptr[1]) +
                  static_cast<ptrdiff_t>(threadid / (blockDim.x * blockDim.y)) *
                      m_scratch_size[1]),
              m_scratch_size[1], league_rank, m_league_size),
          value);
    }
  }

  int compute_block_count() const {
    constexpr auto light_weight =
        Kokkos::Experimental::WorkItemProperty::HintLightWeight;
    constexpr typename Policy::work_item_property property;
    // Numbers were tuned on MI210 using dot product and yAx benchmarks
    constexpr int block_max =
        (property & light_weight) == light_weight ? 2097152 : 65536;
    constexpr int preferred_block_min = 1024;
    int block_count                   = m_league_size;
    if (block_count < preferred_block_min) {
      // keep blocks as is, already low parallelism
    } else if (block_count >= block_max) {
      block_count = block_max;

    } else {
      int nwork = m_league_size * m_team_size;
      int items_per_thread =
          (nwork + block_count * m_team_size - 1) / (block_count * m_team_size);
      if (items_per_thread < 4) {
        int ratio = std::min(
            (block_count + preferred_block_min - 1) / preferred_block_min,
            (4 + items_per_thread - 1) / items_per_thread);
        block_count /= ratio;
      }
    }

    return block_count;
  }

 public:
  __device__ inline void operator()() const {
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = hip_get_scratch_index(m_league_size, m_scratch_locks,
                                       m_num_scratch_locks);
    }

    using ReductionTag = std::conditional_t<UseShflReduction, ShflReductionTag,
                                            SHMEMReductionTag>;
    run(ReductionTag{}, threadid);

    if (m_scratch_size[1] > 0) {
      hip_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  __device__ inline void run(SHMEMReductionTag, int const threadid) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    integral_nonzero_constant<size_type, ReducerType::static_value_size() /
                                             sizeof(size_type)> const
        word_count(reducer.value_size() / sizeof(size_type));

    reference_type value =
        reducer.init(kokkos_impl_hip_shared_memory<size_type>() +
                     threadIdx.y * word_count.value);
    // Iterate this block through the league
    iterate_through_league(threadid, value);

    // Reduce with final value at blockDim.y - 1 location.
    bool do_final_reduce = (m_league_size == 0);
    if (!do_final_reduce)
      do_final_reduce =
          hip_single_inter_block_reduce_scan<false, FunctorType, work_tag>(
              reducer, blockIdx.x, gridDim.x,
              kokkos_impl_hip_shared_memory<size_type>(), m_scratch_space,
              m_scratch_flags);
    if (do_final_reduce) {
      // This is the final block with the final result at the final threads'
      // location

      size_type* const shared = kokkos_impl_hip_shared_memory<size_type>() +
                                (blockDim.y - 1) * word_count.value;
      size_type* const global = m_result_ptr_device_accessible
                                    ? reinterpret_cast<size_type*>(m_result_ptr)
                                    : m_scratch_space;

      if (threadIdx.y == 0) {
        reducer.final(reinterpret_cast<value_type*>(shared));
      }

      if (HIPTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  __device__ inline void run(ShflReductionTag, int const threadid) const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    value_type value;
    reducer.init(&value);

    // Iterate this block through the league
    iterate_through_league(threadid, value);

    pointer_type const result =
        m_result_ptr_device_accessible
            ? m_result_ptr
            : reinterpret_cast<pointer_type>(m_scratch_space);

    value_type init;
    reducer.init(&init);
    if (m_league_size == 0) {
      reducer.final(&value);
      *result = value;
    } else if (Impl::hip_inter_block_shuffle_reduction(
                   value, init, reducer, m_scratch_space, result,
                   m_scratch_flags, blockDim.y)) {
      unsigned int const id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        reducer.final(&value);
        *result = value;
      }
    }
  }

  inline void execute() {
    const ReducerType& reducer = m_functor_reducer.get_reducer();

    const bool is_empty_range  = m_league_size == 0 || m_team_size == 0;
    const bool need_device_set = ReducerType::has_init_member_function() ||
                                 ReducerType::has_final_member_function() ||
                                 !m_result_ptr_host_accessible ||
                                 Policy::is_graph_kernel::value ||
                                 !std::is_same<ReducerType, InvalidType>::value;
    if (!is_empty_range || need_device_set) {
      int const block_count = compute_block_count();

      m_scratch_space = hip_internal_scratch_space(
          m_policy.space(), reducer.value_size() * block_count);
      m_scratch_flags =
          hip_internal_scratch_flags(m_policy.space(), sizeof(size_type));

      dim3 block(m_vector_size, m_team_size, 1);
      dim3 grid(block_count, 1, 1);
      if (is_empty_range) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }
      const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

      Impl::hip_parallel_launch<ParallelReduce, launch_bounds>(
          *this, grid, block, shmem_size_total,
          m_policy.space().impl_internal_space_instance(),
          true);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().impl_internal_space_instance()->fence();

        if (m_result_ptr) {
          const int size = reducer.value_size();
          DeepCopy<HostSpace, HIPSpace>(m_result_ptr, m_scratch_space, size);
        }
      }
    } else {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(CombinedFunctorReducerType const& arg_functor_reducer,
                 Policy const& arg_policy, ViewType const& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<HIPSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_host_accessible(
            MemorySpaceAccess<Kokkos::HostSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
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
            : hip_single_inter_block_reduce_scan_shmem<false, work_tag,
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
                  static_cast<std::int64_t>(HIP().concurrency() /
                                            (m_team_size * m_vector_size)),
                  static_cast<std::int64_t>(m_league_size))));
    }

    // The global parallel_reduce does not support vector_length other than 1 at
    // the moment
    if ((arg_policy.impl_vector_length() > 1) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a vector length of "
          "greater than 1 is not currently supported for HIP for dynamic "
          "sized reduction types.");

    if ((m_team_size < HIPTraits::WarpSize) && !UseShflReduction)
      Impl::throw_runtime_exception(
          "Kokkos::parallel_reduce with a TeamPolicy using a team_size smaller "
          "than 64 is not currently supported with HIP for dynamic sized "
          "reduction types.");

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.

    const int shmem_size_total = m_team_begin + m_shmem_begin + m_shmem_size;

    if (!Kokkos::Impl::is_integral_power_of_two(m_team_size) &&
        !UseShflReduction) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > bad team size"));
    }

    if (internal_space_instance->m_maxShmemPerBlock < shmem_size_total) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > requested too much "
                      "L0 scratch memory"));
    }

    size_t max_size = arg_policy.team_size_max(
        arg_functor_reducer.get_functor(), arg_functor_reducer.get_reducer(),
        ParallelReduceTag());
    if (static_cast<int>(m_team_size) > static_cast<int>(max_size)) {
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelReduce< HIP > requested too "
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

#endif

#endif
