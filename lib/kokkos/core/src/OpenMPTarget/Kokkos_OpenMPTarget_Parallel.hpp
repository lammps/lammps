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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Exec.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

#define KOKKOS_IMPL_LOCK_FREE_HIERARCHICAL

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  inline void execute() const { execute_impl<WorkTag>(); }
  /*
    template <class TagType>
    inline typename std::enable_if<std::is_same<TagType, void>::value>::type
    execute_impl() const {
      OpenMPTargetExec::verify_is_process(
          "Kokkos::Experimental::OpenMPTarget parallel_for");
      OpenMPTargetExec::verify_initialized(
          "Kokkos::Experimental::OpenMPTarget parallel_for");
      const typename Policy::member_type begin = m_policy.begin();
      const typename Policy::member_type end   = m_policy.end();

  #pragma omp target teams distribute parallel for map(to: this->m_functor)
      for (int i = begin; i < end; i++) m_functor(i);
    }
  */
  template <class TagType>
  inline void execute_impl() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = m_policy.begin();
    const auto end   = m_policy.end();

    if (end <= begin) return;

    FunctorType a_functor(m_functor);

    if constexpr (std::is_same<TagType, void>::value) {
#pragma omp target teams distribute parallel for map(to : a_functor)
      for (auto i = begin; i < end; i++) a_functor(i);
    } else {
#pragma omp target teams distribute parallel for map(to : a_functor)
      for (auto i = begin; i < end; i++) a_functor(TagType(), i);
    }
  }

  inline ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class PolicyType, class ReducerType,
          class PointerType, class ValueType, bool FunctorHasJoin,
          bool UseReducerType>
struct ParallelReduceSpecialize {
  static inline void execute(const FunctorType& /*f*/, const PolicyType& /*p*/,
                             PointerType /*result_ptr*/) {
    std::stringstream error_message;
    error_message << "Error: Invalid Specialization " << FunctorHasJoin << ' '
                  << UseReducerType << '\n';
    // FIXME_OPENMPTARGET
    OpenMPTarget_abort(error_message.str().c_str());
  }
};

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, Kokkos::RangePolicy<PolicyArgs...>,
                                ReducerType, PointerType, ValueType, false,
                                false> {
  using PolicyType = Kokkos::RangePolicy<PolicyArgs...>;
  template <class TagType>
  inline static void execute_impl(const FunctorType& f, const PolicyType& p,
                                  PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = p.begin();
    const auto end   = p.end();

    if (end <= begin) return;

    ValueType result = ValueType();
    if constexpr (std::is_same<TagType, void>::value) {
#pragma omp target teams distribute parallel for num_teams(512) \
                map(to:f) map(tofrom:result) reduction(+: result)
      for (auto i = begin; i < end; i++) f(i, result);
    } else {
#pragma omp target teams distribute parallel for num_teams(512) \
                map(to:f) map(tofrom:result) reduction(+: result)
      for (auto i = begin; i < end; i++) f(TagType(), i, result);
    }

    *result_ptr = result;
  }

  inline static void execute(const FunctorType& f, const PolicyType& p,
                             PointerType ptr) {
    execute_impl<typename PolicyType::work_tag>(f, p, ptr);
  }
};

template <class FunctorType, class PolicyType, class ReducerType,
          class PointerType, class ValueType>
struct ParallelReduceSpecialize<FunctorType, PolicyType, ReducerType,
                                PointerType, ValueType, false, true> {
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  template <class TagType>
  inline static void execute_impl(const FunctorType& f, const PolicyType& p,
                                  PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const typename PolicyType::member_type begin = p.begin();
    const typename PolicyType::member_type end   = p.end();

    if (end <= begin) return;

    ValueType result = ValueType();
    OpenMPTargetReducerWrapper<ReducerType>::init(result);

    if constexpr (std::is_same<TagType, void>::value) {
#pragma omp target teams distribute parallel for num_teams(512) map(to   \
                                                                    : f) \
    reduction(custom                                                     \
              : result)
      for (auto i = begin; i < end; i++) f(i, result);
      *result_ptr = result;
    } else {
#pragma omp target teams distribute parallel for num_teams(512) map(to   \
                                                                    : f) \
    reduction(custom                                                     \
              : result)
      for (auto i = begin; i < end; i++) f(TagType(), i, result);
      *result_ptr = result;
    }
  }

  inline static void execute(const FunctorType& f, const PolicyType& p,
                             PointerType ptr) {
    execute_impl<typename PolicyType::work_tag>(f, p, ptr);
  }
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, WorkTag,
                         void>;

  // Static Assert WorkTag void if ReducerType not InvalidType

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

  enum { HasJoin = ReduceFunctorHasJoin<FunctorType>::value };
  enum { UseReducer = is_reducer_type<ReducerType>::value };

  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;

  using ParReduceSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy, ReducerType, pointer_type,
                               typename ValueTraits::value_type, HasJoin,
                               UseReducer>;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;

 public:
  inline void execute() const {
    ParReduceSpecialize::execute(m_functor, m_policy, m_result_ptr);
  }

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& arg_functor, Policy arg_policy,
      const ViewType& arg_result_view,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {}

  inline ParallelReduce(const FunctorType& arg_functor, Policy arg_policy,
                        const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {}
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::OpenMPTarget> {
 protected:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;
  using idx_type  = typename Policy::index_type;

  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, WorkTag>;
  using ValueInit   = Kokkos::Impl::FunctorValueInit<FunctorType, WorkTag>;
  using ValueJoin   = Kokkos::Impl::FunctorValueJoin<FunctorType, WorkTag>;
  using ValueOps    = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

  using value_type     = typename ValueTraits::value_type;
  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline typename std::enable_if<std::is_same<TagType, void>::value>::type
  call_with_tag(const FunctorType& f, const idx_type& idx, value_type& val,
                const bool& is_final) const {
    f(idx, val, is_final);
  }
  template <class TagType>
  inline typename std::enable_if<!std::is_same<TagType, void>::value>::type
  call_with_tag(const FunctorType& f, const idx_type& idx, value_type& val,
                const bool& is_final) const {
    f(WorkTag(), idx, val, is_final);
  }

 public:
  inline void impl_execute(
      Kokkos::View<value_type**, Kokkos::LayoutRight,
                   Kokkos::Experimental::OpenMPTargetSpace>
          element_values,
      Kokkos::View<value_type*, Kokkos::Experimental::OpenMPTargetSpace>
          chunk_values,
      Kokkos::View<int64_t, Kokkos::Experimental::OpenMPTargetSpace> count)
      const {
    const idx_type N          = m_policy.end() - m_policy.begin();
    const idx_type chunk_size = 128;
    const idx_type n_chunks   = (N + chunk_size - 1) / chunk_size;
    idx_type nteams           = n_chunks > 512 ? 512 : n_chunks;
    idx_type team_size        = 128;

    FunctorType a_functor(m_functor);
#pragma omp target teams distribute map(to                             \
                                        : a_functor) num_teams(nteams) \
    thread_limit(team_size)
    for (idx_type team_id = 0; team_id < n_chunks; team_id++) {
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;

#pragma omp for
        for (idx_type i = 0; i < chunk_size; i++) {
          const idx_type idx = local_offset + i;
          value_type val;
          ValueInit::init(a_functor, &val);
          if (idx < N) call_with_tag<WorkTag>(a_functor, idx, val, false);
          element_values(team_id, i) = val;
        }
#pragma omp barrier
        if (omp_get_thread_num() == 0) {
          value_type sum;
          ValueInit::init(a_functor, &sum);
          for (idx_type i = 0; i < chunk_size; i++) {
            ValueJoin::join(a_functor, &sum, &element_values(team_id, i));
            element_values(team_id, i) = sum;
          }
          chunk_values(team_id) = sum;
        }
#pragma omp barrier
        if (omp_get_thread_num() == 0) {
          if (Kokkos::atomic_fetch_add(&count(), 1) == n_chunks - 1) {
            value_type sum;
            ValueInit::init(a_functor, &sum);
            for (idx_type i = 0; i < n_chunks; i++) {
              ValueJoin::join(a_functor, &sum, &chunk_values(i));
              chunk_values(i) = sum;
            }
          }
        }
      }
    }

#pragma omp target teams distribute map(to                             \
                                        : a_functor) num_teams(nteams) \
    thread_limit(team_size)
    for (idx_type team_id = 0; team_id < n_chunks; team_id++) {
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;
        value_type offset_value;
        if (team_id > 0)
          offset_value = chunk_values(team_id - 1);
        else
          ValueInit::init(a_functor, &offset_value);

#pragma omp for
        for (idx_type i = 0; i < chunk_size; i++) {
          const idx_type idx = local_offset + i;
          value_type local_offset_value;
          if (i > 0) {
            local_offset_value = element_values(team_id, i - 1);
            ValueJoin::join(a_functor, &local_offset_value, &offset_value);
          } else
            local_offset_value = offset_value;
          if (idx < N)
            call_with_tag<WorkTag>(a_functor, idx, local_offset_value, true);
        }
      }
    }
  }

  inline void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const idx_type N          = m_policy.end() - m_policy.begin();
    const idx_type chunk_size = 128;
    const idx_type n_chunks   = (N + chunk_size - 1) / chunk_size;

    // This could be scratch memory per team
    Kokkos::View<value_type**, Kokkos::LayoutRight,
                 Kokkos::Experimental::OpenMPTargetSpace>
        element_values("element_values", n_chunks, chunk_size);
    Kokkos::View<value_type*, Kokkos::Experimental::OpenMPTargetSpace>
        chunk_values("chunk_values", n_chunks);
    Kokkos::View<int64_t, Kokkos::Experimental::OpenMPTargetSpace> count(
        "Count");

    impl_execute(element_values, chunk_values, count);
  }

  //----------------------------------------

  inline ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}

  //----------------------------------------
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::OpenMPTarget>
    : public ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                          Kokkos::Experimental::OpenMPTarget> {
  using base_t     = ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                              Kokkos::Experimental::OpenMPTarget>;
  using value_type = typename base_t::value_type;
  value_type& m_returnvalue;

 public:
  inline void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const int64_t N        = base_t::m_policy.end() - base_t::m_policy.begin();
    const int chunk_size   = 128;
    const int64_t n_chunks = (N + chunk_size - 1) / chunk_size;

    if (N > 0) {
      // This could be scratch memory per team
      Kokkos::View<value_type**, Kokkos::LayoutRight,
                   Kokkos::Experimental::OpenMPTargetSpace>
          element_values("element_values", n_chunks, chunk_size);
      Kokkos::View<value_type*, Kokkos::Experimental::OpenMPTargetSpace>
          chunk_values("chunk_values", n_chunks);
      Kokkos::View<int64_t, Kokkos::Experimental::OpenMPTargetSpace> count(
          "Count");

      base_t::impl_execute(element_values, chunk_values, count);

      const int size = base_t::ValueTraits::value_size(base_t::m_functor);
      DeepCopy<HostSpace, Kokkos::Experimental::OpenMPTargetSpace>(
          &m_returnvalue, chunk_values.data() + (n_chunks - 1), size);
    } else {
      m_returnvalue = 0;
    }
  }

  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const typename base_t::Policy& arg_policy,
                        ReturnType& arg_returnvalue)
      : base_t(arg_functor, arg_policy), m_returnvalue(arg_returnvalue) {}
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::OpenMPTarget,
                                       Properties...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const int m_shmem_size;

 public:
  inline void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    execute_impl<WorkTag>();
  }

 private:
  template <class TagType>
  inline void execute_impl() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto league_size   = m_policy.league_size();
    const auto team_size     = m_policy.team_size();
    const auto vector_length = m_policy.impl_vector_length();

    const size_t shmem_size_L0 = m_policy.scratch_size(0, team_size);
    const size_t shmem_size_L1 = m_policy.scratch_size(1, team_size);
    OpenMPTargetExec::resize_scratch(team_size, shmem_size_L0, shmem_size_L1);

    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();
    FunctorType a_functor(m_functor);

    // FIXME_OPENMPTARGET - If the team_size is not a multiple of 32, the
    // scratch implementation does not work in the Release or RelWithDebugInfo
    // mode but works in the Debug mode.

    // Maximum active teams possible.
    int max_active_teams = OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size;
    // nteams should not exceed the maximum in-flight teams possible.
    const auto nteams =
        league_size < max_active_teams ? league_size : max_active_teams;

#ifdef KOKKOS_IMPL_LOCK_FREE_HIERARCHICAL
// Performing our own scheduling of teams to avoid separation of code between
// teams-distribute and parallel. Gave a 2x performance boost in test cases with
// the clang compiler. atomic_compare_exchange can be avoided since the standard
// guarantees that the number of teams specified in the `num_teams` clause is
// always less than or equal to the maximum concurrently running teams.
#pragma omp target teams num_teams(nteams) thread_limit(team_size) \
    map(to                                                         \
        : a_functor) is_device_ptr(scratch_ptr)
#pragma omp parallel
    {
      const int blockIdx = omp_get_team_num();
      const int gridDim  = omp_get_num_teams();

      // Iterate through the number of teams until league_size and assign the
      // league_id accordingly
      // Guarantee that the compilers respect the `num_teams` clause
      if (gridDim <= nteams) {
        for (int league_id = blockIdx; league_id < league_size;
             league_id += gridDim) {
          typename Policy::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              blockIdx, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_same<TagType, void>::value)
            m_functor(team);
          else
            m_functor(TagType(), team);
        }
      } else
        Kokkos::abort("`num_teams` clause was not respected.\n");
    }

#else
// Saving the older implementation that uses `atomic_compare_exchange` to
// calculate the shared memory block index and `distribute` clause to distribute
// teams.
#pragma omp target teams distribute map(to                   \
                                        : a_functor)         \
    is_device_ptr(scratch_ptr, lock_array) num_teams(nteams) \
        thread_limit(team_size)
    for (int i = 0; i < league_size; i++) {
      int shmem_block_index = -1, lock_team = 99999, iter = -1;
      iter = (omp_get_team_num() % max_active_teams);

      // Loop as long as a shmem_block_index is not found.
      while (shmem_block_index == -1) {
        // Try and acquire a lock on the index.
        lock_team = atomic_compare_exchange(&lock_array[iter], 0, 1);

        // If lock is acquired assign it to the block index.
        // lock_team = 0, implies atomic_compare_exchange is successfull.
        if (lock_team == 0)
          shmem_block_index = iter;
        else
          iter = ++iter % max_active_teams;
      }

#pragma omp parallel num_threads(team_size)
      {
        typename Policy::member_type team(
            i, league_size, team_size, vector_length, scratch_ptr,
            shmem_block_index, shmem_size_L0, shmem_size_L1);
        m_functor(team);
      }

      // Free the locked block and increment the number of available free
      // blocks.
      lock_team = atomic_compare_exchange(&lock_array[shmem_block_index], 1, 0);
    }
#endif
  }

 public:
  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {}
};

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, TeamPolicyInternal<PolicyArgs...>,
                                ReducerType, PointerType, ValueType, false,
                                false> {
  using PolicyType = TeamPolicyInternal<PolicyArgs...>;

  template <class TagType>
  inline static void execute_impl(const FunctorType& f, const PolicyType& p,
                                  PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");

    const int league_size   = p.league_size();
    const int team_size     = p.team_size();
    const int vector_length = p.impl_vector_length();

    const size_t shmem_size_L0 = p.scratch_size(0, team_size);
    const size_t shmem_size_L1 = p.scratch_size(1, team_size);
    OpenMPTargetExec::resize_scratch(PolicyType::member_type::TEAM_REDUCE_SIZE,
                                     shmem_size_L0, shmem_size_L1);
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

    ValueType result = ValueType();

    // Maximum active teams possible.
    int max_active_teams = OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size;
    const auto nteams =
        league_size < max_active_teams ? league_size : max_active_teams;

#ifdef KOKKOS_IMPL_LOCK_FREE_HIERARCHICAL
#pragma omp target teams num_teams(nteams) thread_limit(team_size) map(to   \
                                                                       : f) \
    is_device_ptr(scratch_ptr) reduction(+: result)
#pragma omp parallel reduction(+ : result)
    {
      const int blockIdx = omp_get_team_num();
      const int gridDim  = omp_get_num_teams();

      // Guarantee that the compilers respect the `num_teams` clause
      if (gridDim <= nteams) {
        for (int league_id = blockIdx; league_id < league_size;
             league_id += gridDim) {
          typename PolicyType::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              blockIdx, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_same<TagType, void>::value)
            f(team, result);
          else
            f(TagType(), team, result);
        }
      } else
        Kokkos::abort("`num_teams` clause was not respected.\n");
    }

    *result_ptr = result;
#else
// Saving the older implementation that uses `atomic_compare_exchange` to
// calculate the shared memory block index and `distribute` clause to distribute
// teams.
#pragma omp target teams distribute num_teams(nteams) thread_limit(team_size) \
         map(to:f) map(tofrom:result) reduction(+: result) \
    is_device_ptr(scratch_ptr, lock_array)
    for (int i = 0; i < league_size; i++) {
      ValueType inner_result = ValueType();
      int shmem_block_index = -1, lock_team = 99999, iter = -1;
      iter = (omp_get_team_num() % max_active_teams);

      // Loop as long as a shmem_block_index is not found.
      while (shmem_block_index == -1) {
        // Try and acquire a lock on the index.
        lock_team = atomic_compare_exchange(&lock_array[iter], 0, 1);

        // If lock is acquired assign it to the block index.
        // lock_team = 0, implies atomic_compare_exchange is successfull.
        if (lock_team == 0)
          shmem_block_index = iter;
        else
          iter = ++iter % max_active_teams;
      }
#pragma omp parallel num_threads(team_size) reduction(+ : inner_result)
      {
        typename PolicyType::member_type team(
            i, league_size, team_size, vector_length, scratch_ptr,
            shmem_block_index, shmem_size_L0, shmem_size_L1);
        f(team, inner_result);
      }
      result = inner_result;

      // Free the locked block and increment the number of available free
      // blocks.
      lock_team = atomic_compare_exchange(&lock_array[shmem_block_index], 1, 0);
    }

    *result_ptr = result;
#endif
  }

  inline static void execute(const FunctorType& f, const PolicyType& p,
                             PointerType ptr) {
    execute_impl<typename PolicyType::work_tag>(f, p, ptr);
  }
};

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, TeamPolicyInternal<PolicyArgs...>,
                                ReducerType, PointerType, ValueType, false,
                                true> {
  using PolicyType = TeamPolicyInternal<PolicyArgs...>;
  template <class TagType>
  inline static void execute_impl(const FunctorType& f, const PolicyType& p,
                                  PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");

#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))
    const int league_size      = p.league_size();
    const int team_size        = p.team_size();
    const int vector_length    = p.impl_vector_length();
    const size_t shmem_size_L0 = p.scratch_size(0, team_size);
    const size_t shmem_size_L1 = p.scratch_size(1, team_size);
    OpenMPTargetExec::resize_scratch(team_size, shmem_size_L0, shmem_size_L1);
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

    ValueType result = ValueType();

    // Maximum active teams possible.
    int max_active_teams = OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size;
    const auto nteams =
        league_size < max_active_teams ? league_size : max_active_teams;

#pragma omp target teams num_teams(nteams) thread_limit(team_size) map(to   \
                                                                       : f) \
    is_device_ptr(scratch_ptr) reduction(custom                             \
                                         : result)
#pragma omp parallel reduction(custom : result)
    {
      const int blockIdx = omp_get_team_num();
      const int gridDim  = omp_get_num_teams();

      // Guarantee that the compilers respect the `num_teams` clause
      if (gridDim <= nteams) {
        for (int league_id = blockIdx; league_id < league_size;
             league_id += gridDim) {
          typename PolicyType::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              blockIdx, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_same<TagType, void>::value)
            f(team, result);
          else
            f(TagType(), team, result);
        }
      } else
        Kokkos::abort("`num_teams` clause was not respected.\n");
    }

    *result_ptr = result;
  }

  inline static void execute(const FunctorType& f, const PolicyType& p,
                             PointerType ptr) {
    execute_impl<typename PolicyType::work_tag>(f, p, ptr);
  }
};

template <class FunctorType, class ReducerType, class... Properties>
class ParallelReduce<FunctorType, Kokkos::TeamPolicy<Properties...>,
                     ReducerType, Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::OpenMPTarget,
                                       Properties...>;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, WorkTag,
                         void>;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;
  using value_type     = typename ValueTraits::value_type;

  enum { HasJoin = ReduceFunctorHasJoin<FunctorType>::value };
  enum { UseReducer = is_reducer_type<ReducerType>::value };

  using ParForSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy, ReducerType, pointer_type,
                               typename ValueTraits::value_type, HasJoin,
                               UseReducer>;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const int m_shmem_size;

 public:
  inline void execute() const {
    ParForSpecialize::execute(m_functor, m_policy, m_result_ptr);
  }

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& arg_functor, const Policy& arg_policy,
      const ViewType& arg_result,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {}

  inline ParallelReduce(const FunctorType& arg_functor, Policy arg_policy,
                        const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {}
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <typename iType>
struct TeamThreadRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const iType start;
  const iType end;
  const OpenMPTargetExecTeamMember& team;

  inline TeamThreadRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, iType count)
      : start(0), end(count), team(thread_) {}
  inline TeamThreadRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, iType begin_, iType end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  inline ThreadVectorRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, index_type count)
      : start(0), end(count), team(thread_) {}
  inline ThreadVectorRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, index_type begin_,
      index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  inline TeamVectorRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, index_type count)
      : start(0), end(count), team(thread_) {}
  inline TeamVectorRangeBoundariesStruct(
      const OpenMPTargetExecTeamMember& thread_, index_type begin_,
      index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

}  // namespace Impl

}  // namespace Kokkos
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef KOKKOS_IMPL_LOCK_FREE_HIERARCHICAL
#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */
