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

#ifndef KOKKOS_HIP_TEAM_POLICY_INTERNAL_HPP
#define KOKKOS_HIP_TEAM_POLICY_INTERNAL_HPP

#include <Kokkos_MinMax.hpp>

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

}  // namespace Impl
}  // namespace Kokkos

#endif
