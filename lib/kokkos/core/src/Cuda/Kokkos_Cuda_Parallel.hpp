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

#ifndef KOKKOS_CUDA_PARALLEL_HPP
#define KOKKOS_CUDA_PARALLEL_HPP

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

#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/KokkosExp_IterateTileGPU.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

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
  enum { MAX_WARP = 8 };

  typename traits::execution_space m_space;
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  int m_team_scratch_size[2];
  int m_thread_scratch_size[2];
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
  inline int scratch_size(int level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }
  inline int team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }
  inline int thread_scratch_size(int level) const {
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
    using functor_value_traits =
        Impl::FunctorValueTraits<FunctorType, typename traits::work_tag>;

    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type, typename traits::launch_bounds>::
            get_cuda_func_attributes();
    const int block_size = std::forward<BlockSizeCallable>(block_size_callable)(
        space().impl_internal_space_instance(), attr, f,
        (size_t)impl_vector_length(),
        (size_t)team_scratch_size(0) + 2 * sizeof(double),
        (size_t)thread_scratch_size(0) + sizeof(double) +
            ((functor_value_traits::StaticValueSize != 0)
                 ? 0
                 : functor_value_traits::value_size(f)));
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

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

  ParallelFor()        = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;

  template <class TagType>
  inline __device__
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member i) const {
    m_functor(i);
  }

  template <class TagType>
  inline __device__
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member i) const {
    m_functor(TagType(), i);
  }

 public:
  using functor_type = FunctorType;

  Policy const& get_policy() const { return m_policy; }

  inline __device__ void operator()() const {
    const Member work_stride = blockDim.y * gridDim.x;
    const Member work_end    = m_policy.end();

    for (Member iwork =
             m_policy.begin() + threadIdx.y + blockDim.y * blockIdx.x;
         iwork < work_end;
         iwork = iwork < work_end - work_stride ? iwork + work_stride
                                                : work_end) {
      this->template exec_range<WorkTag>(iwork);
    }
  }

  inline void execute() const {
    const typename Policy::index_type nwork = m_policy.end() - m_policy.begin();

    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelFor,
                           LaunchBounds>::get_cuda_func_attributes();
    const int block_size =
        Kokkos::Impl::cuda_get_opt_block_size<FunctorType, LaunchBounds>(
            m_policy.space().impl_internal_space_instance(), attr, m_functor, 1,
            0, 0);
    KOKKOS_ASSERT(block_size > 0);
    dim3 block(1, block_size, 1);
    dim3 grid(
        std::min(
            typename Policy::index_type((nwork + block.y - 1) / block.y),
            typename Policy::index_type(cuda_internal_maximum_grid_count()[0])),
        1, 1);
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (Kokkos::Impl::CudaInternal::cuda_use_serial_execution()) {
      block = dim3(1, 1, 1);
      grid  = dim3(1, 1, 1);
    }
#endif

    CudaParallelLaunch<ParallelFor, LaunchBounds>(
        *this, grid, block, 0, m_policy.space().impl_internal_space_instance(),
        false);
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

// MDRangePolicy impl
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>, Kokkos::Cuda> {
 public:
  using Policy       = Kokkos::MDRangePolicy<Traits...>;
  using functor_type = FunctorType;

 private:
  using RP               = Policy;
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using LaunchBounds     = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_rp;

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& pol, const Functor&) {
    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelFor,
                           LaunchBounds>::get_cuda_func_attributes();
    auto const& prop = pol.space().cuda_device_prop();
    // Limits due to registers/SM, MDRange doesn't have
    // shared memory constraints
    int const regs_per_sm        = prop.regsPerMultiprocessor;
    int const regs_per_thread    = attr.numRegs;
    int const max_threads_per_sm = regs_per_sm / regs_per_thread;
    return std::min(
        max_threads_per_sm,
        static_cast<int>(Kokkos::Impl::CudaTraits::MaxHierarchicalParallelism));
  }
  Policy const& get_policy() const { return m_rp; }
  inline __device__ void operator()() const {
    Kokkos::Impl::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                    typename Policy::work_tag>(m_rp, m_functor)
        .exec_range();
  }

  inline void execute() const {
    using namespace std;

    if (m_rp.m_num_tiles == 0) return;
    const auto maxblocks = cuda_internal_maximum_grid_count();
    if (RP::rank == 2) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], 1);
      KOKKOS_ASSERT(block.x > 0);
      KOKKOS_ASSERT(block.y > 0);
      const dim3 grid(
          std::min<array_index_type>(
              (m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1) / block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1) / block.y,
              maxblocks[1]),
          1);
      CudaParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0, m_rp.space().impl_internal_space_instance(),
          false);
    } else if (RP::rank == 3) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], m_rp.m_tile[2]);
      KOKKOS_ASSERT(block.x > 0);
      KOKKOS_ASSERT(block.y > 0);
      KOKKOS_ASSERT(block.z > 0);
      const dim3 grid(
          std::min<array_index_type>(
              (m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1) / block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1) / block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_rp.m_upper[2] - m_rp.m_lower[2] + block.z - 1) / block.z,
              maxblocks[2]));
      CudaParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0, m_rp.space().impl_internal_space_instance(),
          false);
    } else if (RP::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1], m_rp.m_tile[2],
                       m_rp.m_tile[3]);
      KOKKOS_ASSERT(block.y > 0);
      KOKKOS_ASSERT(block.z > 0);
      const dim3 grid(
          std::min<array_index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1],
                                     maxblocks[0]),
          std::min<array_index_type>(
              (m_rp.m_upper[2] - m_rp.m_lower[2] + block.y - 1) / block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_rp.m_upper[3] - m_rp.m_lower[3] + block.z - 1) / block.z,
              maxblocks[2]));
      CudaParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0, m_rp.space().impl_internal_space_instance(),
          false);
    } else if (RP::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3], m_rp.m_tile[4]);
      KOKKOS_ASSERT(block.z > 0);
      const dim3 grid(
          std::min<array_index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1],
                                     maxblocks[0]),
          std::min<array_index_type>(m_rp.m_tile_end[2] * m_rp.m_tile_end[3],
                                     maxblocks[1]),
          std::min<array_index_type>(
              (m_rp.m_upper[4] - m_rp.m_lower[4] + block.z - 1) / block.z,
              maxblocks[2]));
      CudaParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0, m_rp.space().impl_internal_space_instance(),
          false);
    } else if (RP::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4,id5 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3],
                       m_rp.m_tile[4] * m_rp.m_tile[5]);
      const dim3 grid(
          std::min<array_index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1],
                                     maxblocks[0]),
          std::min<array_index_type>(m_rp.m_tile_end[2] * m_rp.m_tile_end[3],
                                     maxblocks[1]),
          std::min<array_index_type>(m_rp.m_tile_end[4] * m_rp.m_tile_end[5],
                                     maxblocks[2]));
      CudaParallelLaunch<ParallelFor, LaunchBounds>(
          *this, grid, block, 0, m_rp.space().impl_internal_space_instance(),
          false);
    } else {
      Kokkos::abort("Kokkos::MDRange Error: Exceeded rank bounds with Cuda\n");
    }

  }  // end execute

  //  inline
  ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_rp(arg_policy) {}
};

__device__ inline int64_t cuda_get_scratch_index(Cuda::size_type league_size,
                                                 int32_t* scratch_locks) {
  int64_t threadid = 0;
  __shared__ int64_t base_thread_id;
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int64_t const wraparound_len = Kokkos::Experimental::min(
        int64_t(league_size),
        (int64_t(Kokkos::Impl::g_device_cuda_lock_arrays.n)) /
            (blockDim.x * blockDim.y));
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
  int m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;

  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_team(const Member& member) const {
    m_functor(member);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_team(const Member& member) const {
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
        m_policy.space().impl_internal_space_instance(),
        true);  // copy to device and execute
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

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using WorkRange    = typename Policy::WorkRange;
  using WorkTag      = typename Policy::work_tag;
  using Member       = typename Policy::member_type;
  using LaunchBounds = typename Policy::launch_bounds;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using value_type     = typename ValueTraits::value_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
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
  using word_size_type = typename std::conditional<
      sizeof(value_type) < sizeof(Kokkos::Cuda::size_type),
      typename std::conditional<sizeof(value_type) == 2, int16_t, int8_t>::type,
      Kokkos::Cuda::size_type>::type;
  using index_type   = typename Policy::index_type;
  using reducer_type = ReducerType;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  const bool m_result_ptr_host_accessible;
  word_size_type* m_scratch_space;
  // m_scratch_flags must be of type Cuda::size_type due to use of atomics
  // for tracking metadata in Kokkos_Cuda_ReduceScan.hpp
  Cuda::size_type* m_scratch_flags;
  word_size_type* m_unified_space;

  // Shall we use the shfl based reduction or not (only use it for static sized
  // types of more than 128bit)
  enum {
    UseShflReduction = false
  };  //((sizeof(value_type)>2*sizeof(double)) && ValueTraits::StaticValueSize)
      //};
      // Some crutch to do function overloading
 private:
  using DummyShflReductionType  = double;
  using DummySHMEMReductionType = int;

 public:
  Policy const& get_policy() const { return m_policy; }

  // Make the exec_range calls call to Reduce::DeviceIterateTile
  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update) const {
    m_functor(i, update);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update) const {
    m_functor(TagType(), i, update);
  }

  __device__ inline void operator()() const {
    /*    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType,
      DummySHMEMReductionType>::select(1,1.0) );
      }

      __device__ inline
      void run(const DummySHMEMReductionType& ) const
      {*/
    const integral_nonzero_constant<
        word_size_type, ValueTraits::StaticValueSize / sizeof(word_size_type)>
        word_count(ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(word_size_type));

    {
      reference_type value =
          ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                          kokkos_impl_cuda_shared_memory<word_size_type>() +
                              threadIdx.y * word_count.value);

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmatically
      // equivalent.

      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
           iwork < iwork_end; iwork += blockDim.y) {
        this->template exec_range<WorkTag>(iwork, value);
      }
    }

    // Doing code duplication here to fix issue #3428
    // Suspect optimizer bug??
    // Reduce with final value at blockDim.y - 1 location.
    // Shortcut for length zero reduction
    if (m_policy.begin() == m_policy.end()) {
      // This is the final block with the final result at the final threads'
      // location

      word_size_type* const shared =
          kokkos_impl_cuda_shared_memory<word_size_type>() +
          (blockDim.y - 1) * word_count.value;
      word_size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<word_size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), shared);
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
      // return ;
    }

    if (m_policy.begin() != m_policy.end()) {
      {
        if (cuda_single_inter_block_reduce_scan<false, ReducerTypeFwd,
                                                WorkTagFwd>(
                ReducerConditional::select(m_functor, m_reducer), blockIdx.x,
                gridDim.x, kokkos_impl_cuda_shared_memory<word_size_type>(),
                m_scratch_space, m_scratch_flags)) {
          // This is the final block with the final result at the final threads'
          // location

          word_size_type* const shared =
              kokkos_impl_cuda_shared_memory<word_size_type>() +
              (blockDim.y - 1) * word_count.value;
          word_size_type* const global =
              m_result_ptr_device_accessible
                  ? reinterpret_cast<word_size_type*>(m_result_ptr)
                  : (m_unified_space ? m_unified_space : m_scratch_space);

          if (threadIdx.y == 0) {
            Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
                ReducerConditional::select(m_functor, m_reducer), shared);
          }

          if (CudaTraits::WarpSize < word_count.value) {
            __syncthreads();
          }

          for (unsigned i = threadIdx.y; i < word_count.value;
               i += blockDim.y) {
            global[i] = shared[i];
          }
        }
      }
    }
  }
  /*  __device__ inline
     void run(const DummyShflReductionType&) const
     {
       value_type value;
       ValueInit::init( ReducerConditional::select(m_functor , m_reducer) ,
     &value);
       // Number of blocks is bounded so that the reduction can be limited to
     two passes.
       // Each thread block is given an approximately equal amount of work to
     perform.
       // Accumulate the values for this block.
       // The accumulation ordering does not match the final pass, but is
     arithmatically equivalent.

       const WorkRange range( m_policy , blockIdx.x , gridDim.x );

       for ( Member iwork = range.begin() + threadIdx.y , iwork_end =
     range.end() ; iwork < iwork_end ; iwork += blockDim.y ) { this-> template
     exec_range< WorkTag >( iwork , value );
       }

       pointer_type const result = (pointer_type) (m_unified_space ?
     m_unified_space : m_scratch_space) ;

       int max_active_thread = range.end()-range.begin() < blockDim.y ?
     range.end() - range.begin():blockDim.y;

       max_active_thread = (max_active_thread ==
     0)?blockDim.y:max_active_thread;

      value_type init;
      ValueInit::init( ReducerConditional::select(m_functor , m_reducer) ,
     &init);
       if(Impl::cuda_inter_block_reduction<ReducerTypeFwd,ValueJoin,WorkTagFwd>
              (value,init,ValueJoin(ReducerConditional::select(m_functor ,
     m_reducer)),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
         const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
         if(id==0) {
           Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final(
     ReducerConditional::select(m_functor , m_reducer) , (void*) &value );
           *result = value;
         }
       }
     }*/

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    unsigned n = CudaTraits::WarpSize * 8;
    int shmem_size =
        cuda_single_inter_block_reduce_scan_shmem<false, FunctorType, WorkTag>(
            f, n);
    using closure_type = Impl::ParallelReduce<FunctorType, Policy, ReducerType>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type,
                           LaunchBounds>::get_cuda_func_attributes();
    while (
        (n &&
         (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
          shmem_size)) ||
        (n >
         static_cast<unsigned>(
             Kokkos::Impl::cuda_get_max_block_size<FunctorType, LaunchBounds>(
                 m_policy.space().impl_internal_space_instance(), attr, f, 1,
                 shmem_size, 0)))) {
      n >>= 1;
      shmem_size = cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                             WorkTag>(f, n);
    }
    return n;
  }

  inline void execute() {
    const index_type nwork     = m_policy.end() - m_policy.begin();
    const bool need_device_set = ReduceFunctorHasInit<FunctorType>::value ||
                                 ReduceFunctorHasFinal<FunctorType>::value ||
                                 !m_result_ptr_host_accessible ||
#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
                                 Policy::is_graph_kernel::value ||
#endif
                                 !std::is_same<ReducerType, InvalidType>::value;
    if ((nwork > 0) || need_device_set) {
      const int block_size = local_block_size(m_functor);

      KOKKOS_ASSERT(block_size > 0);

      // TODO: down casting these uses more space than required?
      m_scratch_space = (word_size_type*)cuda_internal_scratch_space(
          m_policy.space(), ValueTraits::value_size(ReducerConditional::select(
                                m_functor, m_reducer)) *
                                block_size /* block_size == max block_count */);

      // Intentionally do not downcast to word_size_type since we use Cuda
      // atomics in Kokkos_Cuda_ReduceScan.hpp
      m_scratch_flags = cuda_internal_scratch_flags(m_policy.space(),
                                                    sizeof(Cuda::size_type));
      m_unified_space =
          reinterpret_cast<word_size_type*>(cuda_internal_scratch_unified(
              m_policy.space(),
              ValueTraits::value_size(
                  ReducerConditional::select(m_functor, m_reducer))));

      // REQUIRED ( 1 , N , 1 )
      dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      dim3 grid(std::min(int(block.y), int((nwork + block.y - 1) / block.y)), 1,
                1);

      // TODO @graph We need to effectively insert this in to the graph
      const int shmem =
          UseShflReduction
              ? 0
              : cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                          WorkTag>(m_functor,
                                                                   block.y);

      if ((nwork == 0)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
          || Kokkos::Impl::CudaInternal::cuda_use_serial_execution()
#endif
      ) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      }

      CudaParallelLaunch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().fence(
            "Kokkos::Impl::ParallelReduce<Cuda, RangePolicy>::execute: Result "
            "Not Device Accessible");

        if (m_result_ptr) {
          if (m_unified_space) {
            const int count = ValueTraits::value_count(
                ReducerConditional::select(m_functor, m_reducer));
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = ValueTraits::value_size(
                ReducerConditional::select(m_functor, m_reducer));
            DeepCopy<HostSpace, CudaSpace>(m_result_ptr, m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<ViewType>::value,
                                         void*>::type = nullptr)
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
        m_unified_space(nullptr) {
    check_reduced_view_shmem_size<WorkTag>(m_policy, m_functor);
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
        m_unified_space(nullptr) {
    check_reduced_view_shmem_size<WorkTag>(m_policy, m_functor);
  }
};

// MDRangePolicy impl
template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Cuda> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;

  using WorkTag      = typename Policy::work_tag;
  using Member       = typename Policy::member_type;
  using LaunchBounds = typename Policy::launch_bounds;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using value_type     = typename ValueTraits::value_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Cuda::size_type;
  using reducer_type   = ReducerType;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;  // used for workrange and nwork
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type* m_unified_space;

  using DeviceIteratePattern = typename Kokkos::Impl::Reduce::DeviceIterateTile<
      Policy::rank, Policy, FunctorType, typename Policy::work_tag,
      reference_type>;

  // Shall we use the shfl based reduction or not (only use it for static sized
  // types of more than 128bit
  static constexpr bool UseShflReduction = false;
  //((sizeof(value_type)>2*sizeof(double)) && ValueTraits::StaticValueSize)
  // Some crutch to do function overloading
 private:
  using DummyShflReductionType  = double;
  using DummySHMEMReductionType = int;

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& pol, const Functor&) {
    cudaFuncAttributes attr =
        CudaParallelLaunch<ParallelReduce,
                           LaunchBounds>::get_cuda_func_attributes();
    auto const& prop = pol.space().cuda_device_prop();
    // Limits due do registers/SM
    int const regs_per_sm        = prop.regsPerMultiprocessor;
    int const regs_per_thread    = attr.numRegs;
    int const max_threads_per_sm = regs_per_sm / regs_per_thread;
    return std::min(
        max_threads_per_sm,
        static_cast<int>(Kokkos::Impl::CudaTraits::MaxHierarchicalParallelism));
  }
  Policy const& get_policy() const { return m_policy; }
  inline __device__ void exec_range(reference_type update) const {
    Kokkos::Impl::Reduce::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                            typename Policy::work_tag,
                                            reference_type>(m_policy, m_functor,
                                                            update)
        .exec_range();
  }

  inline __device__ void operator()() const {
    /*    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType,
      DummySHMEMReductionType>::select(1,1.0) );
      }

      __device__ inline
      void run(const DummySHMEMReductionType& ) const
      {*/
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    {
      reference_type value =
          ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                          kokkos_impl_cuda_shared_memory<size_type>() +
                              threadIdx.y * word_count.value);

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmatically
      // equivalent.

      this->exec_range(value);
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Problem: non power-of-two blockDim
    if (cuda_single_inter_block_reduce_scan<false, ReducerTypeFwd, WorkTagFwd>(
            ReducerConditional::select(m_functor, m_reducer), blockIdx.x,
            gridDim.x, kokkos_impl_cuda_shared_memory<size_type>(),
            m_scratch_space, m_scratch_flags)) {
      // This is the final block with the final result at the final threads'
      // location
      size_type* const shared = kokkos_impl_cuda_shared_memory<size_type>() +
                                (blockDim.y - 1) * word_count.value;
      size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), shared);
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }
  }

  /*  __device__ inline
     void run(const DummyShflReductionType&) const
     {

       value_type value;
       ValueInit::init( ReducerConditional::select(m_functor , m_reducer) ,
     &value);
       // Number of blocks is bounded so that the reduction can be limited to
     two passes.
       // Each thread block is given an approximately equal amount of work to
     perform.
       // Accumulate the values for this block.
       // The accumulation ordering does not match the final pass, but is
     arithmatically equivalent.

       const Member work_part =
         ( ( m_policy.m_num_tiles + ( gridDim.x - 1 ) ) / gridDim.x ); //portion
     of tiles handled by each block

       this-> exec_range( value );

       pointer_type const result = (pointer_type) (m_unified_space ?
     m_unified_space : m_scratch_space) ;

       int max_active_thread = work_part < blockDim.y ? work_part:blockDim.y;
       max_active_thread = (max_active_thread ==
     0)?blockDim.y:max_active_thread;

       value_type init;
       ValueInit::init( ReducerConditional::select(m_functor , m_reducer) ,
     &init);
       if(Impl::cuda_inter_block_reduction<ReducerTypeFwd,ValueJoin,WorkTagFwd>
           (value,init,ValueJoin(ReducerConditional::select(m_functor ,
     m_reducer)),m_scratch_space,result,m_scratch_flags,max_active_thread)) {
         const unsigned id = threadIdx.y*blockDim.x + threadIdx.x;
         if(id==0) {
           Kokkos::Impl::FunctorFinal< ReducerTypeFwd , WorkTagFwd >::final(
     ReducerConditional::select(m_functor , m_reducer) , (void*) &value );
           *result = value;
         }
       }
     }
  */
  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    unsigned n = CudaTraits::WarpSize * 8;
    int shmem_size =
        cuda_single_inter_block_reduce_scan_shmem<false, FunctorType, WorkTag>(
            f, n);
    using closure_type = Impl::ParallelReduce<FunctorType, Policy, ReducerType>;
    cudaFuncAttributes attr =
        CudaParallelLaunch<closure_type,
                           LaunchBounds>::get_cuda_func_attributes();
    while (
        (n &&
         (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
          shmem_size)) ||
        (n >
         static_cast<unsigned>(
             Kokkos::Impl::cuda_get_max_block_size<FunctorType, LaunchBounds>(
                 m_policy.space().impl_internal_space_instance(), attr, f, 1,
                 shmem_size, 0)))) {
      n >>= 1;
      shmem_size = cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                             WorkTag>(f, n);
    }
    return n;
  }

  inline void execute() {
    const auto nwork = m_policy.m_num_tiles;
    if (nwork) {
      int block_size = m_policy.m_prod_tile_dims;
      // CONSTRAINT: Algorithm requires block_size >= product of tile dimensions
      // Nearest power of two
      int exponent_pow_two    = std::ceil(std::log2(block_size));
      block_size              = std::pow(2, exponent_pow_two);
      int suggested_blocksize = local_block_size(m_functor);

      block_size = (block_size > suggested_blocksize)
                       ? block_size
                       : suggested_blocksize;  // Note: block_size must be less
                                               // than or equal to 512

      m_scratch_space = cuda_internal_scratch_space(
          m_policy.space(), ValueTraits::value_size(ReducerConditional::select(
                                m_functor, m_reducer)) *
                                block_size /* block_size == max block_count */);
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type));
      m_unified_space = cuda_internal_scratch_unified(
          m_policy.space(), ValueTraits::value_size(ReducerConditional::select(
                                m_functor, m_reducer)));

      // REQUIRED ( 1 , N , 1 )
      const dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      const dim3 grid(std::min(int(block.y), int(nwork)), 1, 1);

      // TODO @graph We need to effectively insert this in to the graph
      const int shmem =
          UseShflReduction
              ? 0
              : cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                          WorkTag>(m_functor,
                                                                   block.y);

      CudaParallelLaunch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().fence(
            "Kokkos::Impl::ParallelReduce<Cuda, MDRangePolicy>::execute: "
            "Result Not Device Accessible");

        if (m_result_ptr) {
          if (m_unified_space) {
            const int count = ValueTraits::value_count(
                ReducerConditional::select(m_functor, m_reducer));
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = ValueTraits::value_size(
                ReducerConditional::select(m_functor, m_reducer));
            DeepCopy<HostSpace, CudaSpace>(m_result_ptr, m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<ViewType>::value,
                                         void*>::type = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::CudaSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_unified_space(nullptr) {
    check_reduced_view_shmem_size<WorkTag>(m_policy, m_functor);
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
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_unified_space(nullptr) {
    check_reduced_view_shmem_size<WorkTag>(m_policy, m_functor);
  }
};

//----------------------------------------------------------------------------

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

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;
  using value_type     = typename ValueTraits::value_type;

 public:
  using functor_type = FunctorType;
  using size_type    = Cuda::size_type;
  using reducer_type = ReducerType;

  enum : bool {
    UseShflReduction = (true && (ValueTraits::StaticValueSize != 0))
  };

 private:
  using DummyShflReductionType  = double;
  using DummySHMEMReductionType = int;

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
  int m_scratch_size[2];
  int m_scratch_pool_id = -1;
  int32_t* m_scratch_locks;
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;

  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_team(const Member& member, reference_type update) const {
    m_functor(member, update);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_team(const Member& member, reference_type update) const {
    m_functor(TagType(), member, update);
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  __device__ inline void operator()() const {
    int64_t threadid = 0;
    if (m_scratch_size[1] > 0) {
      threadid = cuda_get_scratch_index(m_league_size, m_scratch_locks);
    }

    run(Kokkos::Impl::if_c<UseShflReduction, DummyShflReductionType,
                           DummySHMEMReductionType>::select(1, 1.0),
        threadid);
    if (m_scratch_size[1] > 0) {
      cuda_release_scratch_index(m_scratch_locks, threadid);
    }
  }

  __device__ inline void run(const DummySHMEMReductionType&,
                             const int& threadid) const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    reference_type value =
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        kokkos_impl_cuda_shared_memory<size_type>() +
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
    // Doing code duplication here to fix issue #3428
    // Suspect optimizer bug??
    if (m_league_size == 0) {
      // This is the final block with the final result at the final threads'
      // location

      size_type* const shared = kokkos_impl_cuda_shared_memory<size_type>() +
                                (blockDim.y - 1) * word_count.value;
      size_type* const global =
          m_result_ptr_device_accessible
              ? reinterpret_cast<size_type*>(m_result_ptr)
              : (m_unified_space ? m_unified_space : m_scratch_space);

      if (threadIdx.y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), shared);
      }

      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }

      for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
        global[i] = shared[i];
      }
    }

    if (m_league_size != 0) {
      if (cuda_single_inter_block_reduce_scan<false, FunctorType, WorkTag>(
              ReducerConditional::select(m_functor, m_reducer), blockIdx.x,
              gridDim.x, kokkos_impl_cuda_shared_memory<size_type>(),
              m_scratch_space, m_scratch_flags)) {
        // This is the final block with the final result at the final threads'
        // location

        size_type* const shared = kokkos_impl_cuda_shared_memory<size_type>() +
                                  (blockDim.y - 1) * word_count.value;
        size_type* const global =
            m_result_ptr_device_accessible
                ? reinterpret_cast<size_type*>(m_result_ptr)
                : (m_unified_space ? m_unified_space : m_scratch_space);

        if (threadIdx.y == 0) {
          Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
              ReducerConditional::select(m_functor, m_reducer), shared);
        }

        if (CudaTraits::WarpSize < word_count.value) {
          __syncthreads();
        }

        for (unsigned i = threadIdx.y; i < word_count.value; i += blockDim.y) {
          global[i] = shared[i];
        }
      }
    }
  }

  __device__ inline void run(const DummyShflReductionType&,
                             const int& threadid) const {
    value_type value;
    ValueInit::init(ReducerConditional::select(m_functor, m_reducer), &value);

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
    ValueInit::init(ReducerConditional::select(m_functor, m_reducer), &init);

    if (int_league_size == 0) {
      Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
          ReducerConditional::select(m_functor, m_reducer), (void*)&value);
      *result = value;
    } else if (
        Impl::cuda_inter_block_reduction<FunctorType, ValueJoin, WorkTag>(
            value, init,
            ValueJoin(ReducerConditional::select(m_functor, m_reducer)),
            m_scratch_space, result, m_scratch_flags, blockDim.y)
        // This breaks a test
        //   Kokkos::Impl::CudaReductionsFunctor<FunctorType,WorkTag,false,true>::scalar_inter_block_reduction(ReducerConditional::select(m_functor
        //   , m_reducer) , blockIdx.x , gridDim.x ,
        //              kokkos_impl_cuda_shared_memory<size_type>() ,
        //              m_scratch_space , m_scratch_flags)
    ) {
      const unsigned id = threadIdx.y * blockDim.x + threadIdx.x;
      if (id == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), (void*)&value);
        *result = value;
      }
    }
  }

  inline void execute() {
    const bool is_empty_range  = m_league_size == 0 || m_team_size == 0;
    const bool need_device_set = ReduceFunctorHasInit<FunctorType>::value ||
                                 ReduceFunctorHasFinal<FunctorType>::value ||
                                 !m_result_ptr_host_accessible ||
#ifdef KOKKOS_CUDA_ENABLE_GRAPHS
                                 Policy::is_graph_kernel::value ||
#endif
                                 !std::is_same<ReducerType, InvalidType>::value;
    if (!is_empty_range || need_device_set) {
      const int block_count =
          UseShflReduction ? std::min(m_league_size, size_type(1024 * 32))
                           : std::min(int(m_league_size), m_team_size);

      m_scratch_space = cuda_internal_scratch_space(
          m_policy.space(), ValueTraits::value_size(ReducerConditional::select(
                                m_functor, m_reducer)) *
                                block_count);
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type));
      m_unified_space = cuda_internal_scratch_unified(
          m_policy.space(), ValueTraits::value_size(ReducerConditional::select(
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
          m_policy.space().impl_internal_space_instance(),
          true);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        m_policy.space().fence(
            "Kokkos::Impl::ParallelReduce<Cuda, TeamPolicy>::execute: Result "
            "Not Device Accessible");

        if (m_result_ptr) {
          if (m_unified_space) {
            const int count = ValueTraits::value_count(
                ReducerConditional::select(m_functor, m_reducer));
            for (int i = 0; i < count; ++i) {
              m_result_ptr[i] = pointer_type(m_unified_space)[i];
            }
          } else {
            const int size = ValueTraits::value_size(
                ReducerConditional::select(m_functor, m_reducer));
            DeepCopy<HostSpace, CudaSpace>(m_result_ptr, m_scratch_space, size);
          }
        }
      }
    } else {
      if (m_result_ptr) {
        // TODO @graph We need to effectively insert this in to the graph
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<ViewType>::value,
                                         void*>::type = nullptr)
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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, WorkTag>;
  using ValueInit   = Kokkos::Impl::FunctorValueInit<FunctorType, WorkTag>;
  using ValueOps    = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Cuda::size_type;

 private:
  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type m_final;
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  bool m_run_serial;
#endif

  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(i, update, final_result);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(TagType(), i, update, final_result);
  }

  //----------------------------------------

  __device__ inline void initial() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    size_type* const shared_value =
        kokkos_impl_cuda_shared_memory<size_type>() +
        word_count.value * threadIdx.y;

    ValueInit::init(m_functor, shared_value);

    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmatically equivalent.

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
         iwork < iwork_end; iwork += blockDim.y) {
      this->template exec_range<WorkTag>(
          iwork, ValueOps::reference(shared_value), false);
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups'
    // totals. Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i <
    // gridDim.x
    cuda_single_inter_block_reduce_scan<true, FunctorType, WorkTag>(
        m_functor, blockIdx.x, gridDim.x,
        kokkos_impl_cuda_shared_memory<size_type>(), m_scratch_space,
        m_scratch_flags);
  }

  //----------------------------------------

  __device__ inline void final() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] ,
    // value[2] , ... }
    size_type* const shared_data = kokkos_impl_cuda_shared_memory<size_type>();
    size_type* const shared_prefix =
        shared_data + word_count.value * threadIdx.y;
    size_type* const shared_accum =
        shared_data + word_count.value * (blockDim.y + 1);

    // Starting value for this thread block is the previous block's total.
    if (blockIdx.x) {
      size_type* const block_total =
          m_scratch_space + word_count.value * (blockIdx.x - 1);
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_accum[i] = block_total[i];
      }
    } else if (0 == threadIdx.y) {
      ValueInit::init(m_functor, shared_accum);
    }

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (typename Policy::member_type iwork_base = range.begin();
         iwork_base < range.end(); iwork_base += blockDim.y) {
      unsigned MASK                            = __activemask();
      const typename Policy::member_type iwork = iwork_base + threadIdx.y;

      __syncthreads();  // Don't overwrite previous iteration values until they
                        // are used

      ValueInit::init(m_functor, shared_prefix + word_count.value);

      // Copy previous block's accumulation total into thread[0] prefix and
      // inclusive scan value of this block
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i];
      }
      __syncwarp(MASK);
      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }  // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix + word_count.value),
            false);
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true, FunctorType, WorkTag>(
          m_functor,
          typename ValueTraits::pointer_type(shared_data + word_count.value));

      {
        size_type* const block_total =
            shared_data + word_count.value * blockDim.y;
        for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
          shared_accum[i] = block_total[i];
        }
      }

      // Call functor with exclusive scan value
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix), true);
      }
    }
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  //----------------------------------------

  __device__ inline void operator()() const {
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (m_run_serial) {
      typename ValueTraits::value_type value;
      ValueInit::init(m_functor, (void*)&value);
      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (typename Policy::member_type iwork_base = range.begin();
           iwork_base < range.end(); iwork_base++) {
        this->template exec_range<WorkTag>(iwork_base, value, true);
      }
    } else {
#endif
      if (!m_final) {
        initial();
      } else {
        final();
      }
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    }
#endif
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512
    // (16 warps) gridDim.x <= blockDim.y * blockDim.y
    //
    // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit
    // testing

    unsigned n = CudaTraits::WarpSize * 4;
    while (n &&
           unsigned(m_policy.space()
                        .impl_internal_space_instance()
                        ->m_maxShmemPerBlock) <
               cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                         WorkTag>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void execute() {
    const auto nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      enum { GridMaxComputeCapability_2x = 0x0ffff };

      const int block_size = local_block_size(m_functor);
      KOKKOS_ASSERT(block_size > 0);

      const int grid_max =
          (block_size * block_size) < GridMaxComputeCapability_2x
              ? (block_size * block_size)
              : GridMaxComputeCapability_2x;

      // At most 'max_grid' blocks:
      const int max_grid =
          std::min(int(grid_max), int((nwork + block_size - 1) / block_size));

      // How much work per block:
      const int work_per_block = (nwork + max_grid - 1) / max_grid;

      // How many block are really needed for this much work:
      const int grid_x = (nwork + work_per_block - 1) / work_per_block;

      m_scratch_space = cuda_internal_scratch_space(
          m_policy.space(), ValueTraits::value_size(m_functor) * grid_x);
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type) * 1);

      dim3 grid(grid_x, 1, 1);
      dim3 block(1, block_size, 1);  // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = ValueTraits::value_size(m_functor) * (block_size + 2);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      } else {
#endif
        m_final = false;
        CudaParallelLaunch<ParallelScan, LaunchBounds>(
            *this, grid, block, shmem,
            m_policy.space().impl_internal_space_instance(),
            false);  // copy to device and execute
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      }
#endif
      m_final = true;
      CudaParallelLaunch<ParallelScan, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute
    }
  }

  ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_final(false)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        ,
        m_run_serial(Kokkos::Impl::CudaInternal::cuda_use_serial_execution())
#endif
  {
  }
};

//----------------------------------------------------------------------------
template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Cuda> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, WorkTag>;
  using ValueInit   = Kokkos::Impl::FunctorValueInit<FunctorType, WorkTag>;
  using ValueOps    = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Cuda::size_type;

 private:
  // Algorithmic constraints:
  //  (a) blockDim.y is a power of two
  //  (b) blockDim.y == blockDim.z == 1
  //  (c) gridDim.x  <= blockDim.y * blockDim.y
  //  (d) gridDim.y  == gridDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;
  size_type* m_scratch_space;
  size_type* m_scratch_flags;
  size_type m_final;
  ReturnType& m_returnvalue;
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
  bool m_run_serial;
#endif

  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(i, update, final_result);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(TagType(), i, update, final_result);
  }

  //----------------------------------------

  __device__ inline void initial() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    size_type* const shared_value =
        kokkos_impl_cuda_shared_memory<size_type>() +
        word_count.value * threadIdx.y;

    ValueInit::init(m_functor, shared_value);

    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmatically equivalent.

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (Member iwork = range.begin() + threadIdx.y, iwork_end = range.end();
         iwork < iwork_end; iwork += blockDim.y) {
      this->template exec_range<WorkTag>(
          iwork, ValueOps::reference(shared_value), false);
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups'
    // totals. Blocks' scan values are written to 'blockIdx.x' location.
    // Block-groups' scan values are at: i = ( j * blockDim.y - 1 ) for i <
    // gridDim.x
    cuda_single_inter_block_reduce_scan<true, FunctorType, WorkTag>(
        m_functor, blockIdx.x, gridDim.x,
        kokkos_impl_cuda_shared_memory<size_type>(), m_scratch_space,
        m_scratch_flags);
  }

  //----------------------------------------

  __device__ inline void final() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] ,
    // value[2] , ... }
    size_type* const shared_data = kokkos_impl_cuda_shared_memory<size_type>();
    size_type* const shared_prefix =
        shared_data + word_count.value * threadIdx.y;
    size_type* const shared_accum =
        shared_data + word_count.value * (blockDim.y + 1);

    // Starting value for this thread block is the previous block's total.
    if (blockIdx.x) {
      size_type* const block_total =
          m_scratch_space + word_count.value * (blockIdx.x - 1);
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_accum[i] = block_total[i];
      }
    } else if (0 == threadIdx.y) {
      ValueInit::init(m_functor, shared_accum);
    }

    const WorkRange range(m_policy, blockIdx.x, gridDim.x);

    for (typename Policy::member_type iwork_base = range.begin();
         iwork_base < range.end(); iwork_base += blockDim.y) {
      unsigned MASK = __activemask();

      const typename Policy::member_type iwork = iwork_base + threadIdx.y;

      __syncthreads();  // Don't overwrite previous iteration values until they
                        // are used

      ValueInit::init(m_functor, shared_prefix + word_count.value);

      // Copy previous block's accumulation total into thread[0] prefix and
      // inclusive scan value of this block
      for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i];
      }

      __syncwarp(MASK);
      if (CudaTraits::WarpSize < word_count.value) {
        __syncthreads();
      }  // Protect against large scan values.

      // Call functor to accumulate inclusive scan value for this work item
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix + word_count.value),
            false);
      }

      // Scan block values into locations shared_data[1..blockDim.y]
      cuda_intra_block_reduce_scan<true, FunctorType, WorkTag>(
          m_functor,
          typename ValueTraits::pointer_type(shared_data + word_count.value));

      {
        size_type* const block_total =
            shared_data + word_count.value * blockDim.y;
        for (unsigned i = threadIdx.y; i < word_count.value; ++i) {
          shared_accum[i] = block_total[i];
        }
      }

      // Call functor with exclusive scan value
      if (iwork < range.end()) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix), true);
      }
    }
  }

 public:
  Policy const& get_policy() const { return m_policy; }

  //----------------------------------------

  __device__ inline void operator()() const {
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    if (m_run_serial) {
      typename ValueTraits::value_type value;
      ValueInit::init(m_functor, (void*)&value);
      const WorkRange range(m_policy, blockIdx.x, gridDim.x);

      for (typename Policy::member_type iwork_base = range.begin();
           iwork_base < range.end(); iwork_base++) {
        this->template exec_range<WorkTag>(iwork_base, value, true);
      }
      *((typename ValueTraits::value_type*)m_scratch_space) = value;
    } else {
#endif
      if (!m_final) {
        initial();
      } else {
        final();
      }
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    }
#endif
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // blockDim.y must be power of two = 128 (4 warps) or 256 (8 warps) or 512
    // (16 warps) gridDim.x <= blockDim.y * blockDim.y
    //
    // 4 warps was 10% faster than 8 warps and 20% faster than 16 warps in unit
    // testing

    unsigned n = CudaTraits::WarpSize * 4;
    while (n &&
           unsigned(m_policy.space()
                        .impl_internal_space_instance()
                        ->m_maxShmemPerBlock) <
               cuda_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                         WorkTag>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void execute() {
    const auto nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      enum { GridMaxComputeCapability_2x = 0x0ffff };

      const int block_size = local_block_size(m_functor);
      KOKKOS_ASSERT(block_size > 0);

      const int grid_max =
          (block_size * block_size) < GridMaxComputeCapability_2x
              ? (block_size * block_size)
              : GridMaxComputeCapability_2x;

      // At most 'max_grid' blocks:
      const int max_grid =
          std::min(int(grid_max), int((nwork + block_size - 1) / block_size));

      // How much work per block:
      const int work_per_block = (nwork + max_grid - 1) / max_grid;

      // How many block are really needed for this much work:
      const int grid_x = (nwork + work_per_block - 1) / work_per_block;

      m_scratch_space = cuda_internal_scratch_space(
          m_policy.space(), ValueTraits::value_size(m_functor) * grid_x);
      m_scratch_flags =
          cuda_internal_scratch_flags(m_policy.space(), sizeof(size_type) * 1);

      dim3 grid(grid_x, 1, 1);
      dim3 block(1, block_size, 1);  // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = ValueTraits::value_size(m_functor) * (block_size + 2);

#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial) {
        block = dim3(1, 1, 1);
        grid  = dim3(1, 1, 1);
      } else {
#endif

        m_final = false;
        CudaParallelLaunch<ParallelScanWithTotal, LaunchBounds>(
            *this, grid, block, shmem,
            m_policy.space().impl_internal_space_instance(),
            false);  // copy to device and execute
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      }
#endif
      m_final = true;
      CudaParallelLaunch<ParallelScanWithTotal, LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      const int size = ValueTraits::value_size(m_functor);
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
      if (m_run_serial)
        DeepCopy<HostSpace, CudaSpace>(&m_returnvalue, m_scratch_space, size);
      else
#endif
        DeepCopy<HostSpace, CudaSpace>(
            &m_returnvalue, m_scratch_space + (grid_x - 1) * size / sizeof(int),
            size);
    }
  }

  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const Policy& arg_policy, ReturnType& arg_returnvalue)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_scratch_space(nullptr),
        m_scratch_flags(nullptr),
        m_final(false),
        m_returnvalue(arg_returnvalue)
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        ,
        m_run_serial(Kokkos::Impl::CudaInternal::cuda_use_serial_execution())
#endif
  {
  }
};

}  // namespace Impl

}  // namespace Kokkos

#endif /* defined(KOKKOS_ENABLE_CUDA) */
#endif /* #ifndef KOKKOS_CUDA_PARALLEL_HPP */
