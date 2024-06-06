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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_COMMON_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_COMMON_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Reducer.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Macros.hpp>

namespace Kokkos {
namespace Impl {

// This class has the memcpy routine that is commonly used by ParallelReduce
// over RangePolicy and TeamPolicy.
template <class PointerType>
struct ParallelReduceCopy {
  // Copy the result back to device if the view is on the device.
  static void memcpy_result(PointerType dest, PointerType src, size_t size,
                            bool ptr_on_device) {
    if (ptr_on_device) {
      if (0 < size) {
        KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(dest, src, size, 0, 0,
                                                     omp_get_default_device(),
                                                     omp_get_initial_device()));
      }

    } else {
      *dest = *src;
    }
  }
};

// template <class FunctorType, class PolicyType, class ReducerType,
// class PointerType, class ValueType>
template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class PolicyType>
struct ParallelReduceSpecialize {
  inline static void execute(const FunctorType& /*f*/, const PolicyType& /*p*/,
                             PointerType /*result_ptr*/) {
    constexpr int FunctorHasJoin =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                              FunctorType,
                              ValueType>::Reducer::has_join_member_function();
    constexpr int UseReducerType = is_reducer_v<ReducerType>;

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
                                ReducerType, PointerType, ValueType> {
  using PolicyType = Kokkos::RangePolicy<PolicyArgs...>;
  using TagType    = typename PolicyType::work_tag;
  using ReducerTypeFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                         PolicyType, ReducerTypeFwd, ValueType>;
  using ReferenceType = typename Analysis::reference_type;

  using ParReduceCopy = ParallelReduceCopy<PointerType>;

  static void execute_reducer(const FunctorType& f, const PolicyType& p,
                              PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:reducer");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:reducer");
    const auto begin = p.begin();
    const auto end   = p.end();

    ValueType result;
    OpenMPTargetReducerWrapper<ReducerType>::init(result);

    // Initialize and copy back the result even if it is a zero length
    // reduction.
    if (end <= begin) {
      ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                   ptr_on_device);
      return;
    }

#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp target teams distribute parallel for map(to                    \
                                                     : f) reduction(custom \
                                                                    : result)
    for (auto i = begin; i < end; ++i) {
      if constexpr (std::is_void_v<TagType>) {
        f(i, result);
      } else {
        f(TagType(), i, result);
      }
    }

    ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                 ptr_on_device);
  }

  template <class TagType, int NumReductions>
  static void execute_array(const FunctorType& f, const PolicyType& p,
                            PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:array_reduction");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:array_reduction");
    const auto begin = p.begin();
    const auto end   = p.end();

    // Enter the loop if the reduction is on a scalar type.
    if constexpr (NumReductions == 1) {
      ValueType result = ValueType();

      // Initialize and copy back the result even if it is a zero length
      // reduction.
      if (end <= begin) {
        ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                     ptr_on_device);
        return;
      }

      // Case where reduction is on a native data type.
      if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp target teams distribute parallel for \
         map(to:f) reduction(+: result)
        for (auto i = begin; i < end; ++i)

          if constexpr (std::is_void_v<TagType>) {
            f(i, result);
          } else {
            f(TagType(), i, result);
          }
      } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)
#pragma omp target teams distribute parallel for map(to                    \
                                                     : f) reduction(custom \
                                                                    : result)
        for (auto i = begin; i < end; ++i)

          if constexpr (std::is_void_v<TagType>) {
            f(i, result);
          } else {
            f(TagType(), i, result);
          }
      }

      ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                   ptr_on_device);
    } else {
      ValueType result[NumReductions] = {};

      // Initialize and copy back the result even if it is a zero length
      // reduction.
      if (end <= begin) {
        ParReduceCopy::memcpy_result(result_ptr, result,
                                     NumReductions * sizeof(ValueType),
                                     ptr_on_device);
        return;
      }
#pragma omp target teams distribute parallel for map(to:f) reduction(+:result[:NumReductions])
      for (auto i = begin; i < end; ++i) {
        if constexpr (std::is_void_v<TagType>) {
          f(i, result);
        } else {
          f(TagType(), i, result);
        }
      }

      ParReduceCopy::memcpy_result(
          result_ptr, result, NumReductions * sizeof(ValueType), ptr_on_device);
    }
  }

  static void execute_init_join(const FunctorType& f, const PolicyType& p,
                                PointerType ptr, const bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:init_join");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget RangePolicy "
        "parallel_reduce:init_join");
    const auto begin = p.begin();
    const auto end   = p.end();

    using FunctorAnalysis =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                              FunctorType, ValueType>;

    // Initialize the result pointer.

    const auto size = end - begin;

    // FIXME_OPENMPTARGET: The team size and MAX_ACTIVE_THREADS are currently
    // based on NVIDIA-V100 and should be modifid to be based on the
    // architecture in the future.
    const int max_team_threads = 32;
    const int max_teams =
        OpenMPTargetExec::MAX_ACTIVE_THREADS / max_team_threads;
    // Number of elements in the reduction
    const auto value_count = FunctorAnalysis::value_count(f);

    // Allocate scratch per active thread. Achieved by setting the first
    // parameter of `resize_scratch=1`.
    OpenMPTargetExec::resize_scratch(1, 0, value_count * sizeof(ValueType),
                                     std::numeric_limits<int64_t>::max());
    ValueType* scratch_ptr =
        static_cast<ValueType*>(OpenMPTargetExec::get_scratch_ptr());

    typename FunctorAnalysis::Reducer final_reducer(f);

    if (end <= begin) {
#pragma omp target map(to : final_reducer) is_device_ptr(scratch_ptr)
      {
        // If there is no work to be done, copy back the initialized values and
        // exit.
        final_reducer.init(scratch_ptr);
        final_reducer.final(scratch_ptr);
      }
      if (0 < value_count) {
        if (!ptr_on_device)
          KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
              ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
              omp_get_initial_device(), omp_get_default_device()));
        else
          KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
              ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
              omp_get_default_device(), omp_get_default_device()));
      }

      return;
    }

#pragma omp target teams num_teams(max_teams) thread_limit(max_team_threads) \
    map(to                                                                   \
        : final_reducer) is_device_ptr(scratch_ptr)
    {
#pragma omp parallel
      {
        const int team_num    = omp_get_team_num();
        const int num_teams   = omp_get_num_teams();
        const auto chunk_size = size / num_teams;
        const auto team_begin = begin + team_num * chunk_size;
        const auto team_end =
            (team_num == num_teams - 1) ? end : (team_begin + chunk_size);
        ValueType* team_scratch =
            scratch_ptr + team_num * max_team_threads * value_count;
        ReferenceType result = final_reducer.init(
            &team_scratch[omp_get_thread_num() * value_count]);

        // Accumulate partial results in thread specific storage.
#pragma omp for simd
        for (auto i = team_begin; i < team_end; ++i) {
          if constexpr (std::is_void_v<TagType>) {
            f(i, result);
          } else {
            f(TagType(), i, result);
          }
        }

        // Reduce all paritial results within a team.
        const int team_size      = max_team_threads;
        int tree_neighbor_offset = 1;
        do {
#pragma omp for simd
          for (int i = 0; i < team_size - tree_neighbor_offset;
               i += 2 * tree_neighbor_offset) {
            const int neighbor = i + tree_neighbor_offset;
            final_reducer.join(&team_scratch[i * value_count],
                               &team_scratch[neighbor * value_count]);
          }
          tree_neighbor_offset *= 2;
        } while (tree_neighbor_offset < team_size);
      }  // end parallel
    }    // end target

    int tree_neighbor_offset = 1;
    do {
#pragma omp target teams distribute parallel for simd map(to   \
                                                          : f) \
    is_device_ptr(scratch_ptr)
      for (int i = 0; i < max_teams - tree_neighbor_offset;
           i += 2 * tree_neighbor_offset) {
        ValueType* team_scratch = scratch_ptr;
        const int team_offset   = max_team_threads * value_count;
        final_reducer.join(
            &team_scratch[i * team_offset],
            &team_scratch[(i + tree_neighbor_offset) * team_offset]);

        // If `final` is provided by the functor.
        // Do the final only once at the end.
        if (tree_neighbor_offset * 2 >= max_teams && omp_get_team_num() == 0 &&
            omp_get_thread_num() == 0) {
          final_reducer.final(scratch_ptr);
        }
      }
      tree_neighbor_offset *= 2;
    } while (tree_neighbor_offset < max_teams);

    // If the result view is on the host, copy back the values via memcpy.
    if (0 < value_count) {
      if (!ptr_on_device)
        KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_initial_device(), omp_get_default_device()));
      else
        KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_default_device(), omp_get_default_device()));
    }
  }
};

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, TeamPolicyInternal<PolicyArgs...>,
                                ReducerType, PointerType, ValueType> {
  using PolicyType = TeamPolicyInternal<PolicyArgs...>;
  using TagType    = typename PolicyType::work_tag;
  using ReducerTypeFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                         PolicyType, ReducerTypeFwd, ValueType>;

  using ReferenceType = typename Analysis::reference_type;

  using ParReduceCopy = ParallelReduceCopy<PointerType>;

  static void execute_reducer(const FunctorType& f, const PolicyType& p,
                              PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:reducer");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:reducer");

    const int league_size   = p.league_size();
    const int team_size     = p.team_size();
    const int vector_length = p.impl_vector_length();

    const size_t shmem_size_L0 = p.scratch_size(0, team_size);
    const size_t shmem_size_L1 = p.scratch_size(1, team_size);
    OpenMPTargetExec::resize_scratch(PolicyType::member_type::TEAM_REDUCE_SIZE,
                                     shmem_size_L0, shmem_size_L1, league_size);
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

    ValueType result = ValueType();

    // Maximum active teams possible.
    // FIXME_OPENMPTARGET: Cray compiler did not yet implement
    // omp_get_max_teams.
#if !defined(KOKKOS_COMPILER_CRAY_LLVM)
    int max_active_teams = omp_get_max_teams();
#else
    int max_active_teams =
        std::min(OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size, league_size);
#endif

    // If the league size is <=0, do not launch the kernel.
    if (max_active_teams <= 0) return;

#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#if !defined(KOKKOS_IMPL_OPENMPTARGET_HIERARCHICAL_INTEL_GPU)
    KOKKOS_IMPL_OMPTARGET_PRAGMA(
        teams num_teams(max_active_teams) thread_limit(team_size)
            firstprivate(f) is_device_ptr(scratch_ptr) reduction(custom
                                                                 : result)
                KOKKOS_IMPL_OMPX_DYN_CGROUP_MEM(shmem_size_L0))
#pragma omp parallel reduction(custom : result)
    {
      if (omp_get_num_teams() > max_active_teams)
        Kokkos::abort("`omp_set_num_teams` call was not respected.\n");

      const int blockIdx = omp_get_team_num();
      const int gridDim  = omp_get_num_teams();

      // Guarantee that the compilers respect the `num_teams` clause
      for (int league_id = blockIdx; league_id < league_size;
           league_id += gridDim) {
        typename PolicyType::member_type team(
            league_id, league_size, team_size, vector_length, scratch_ptr,
            blockIdx, shmem_size_L0, shmem_size_L1);
        if constexpr (std::is_void_v<TagType>)
          f(team, result);
        else
          f(TagType(), team, result);
      }
    }
#else
#pragma omp target teams distribute firstprivate(f) is_device_ptr(scratch_ptr) \
    num_teams(max_active_teams) thread_limit(team_size) reduction(custom       \
                                                                  : result)
    for (int i = 0; i < league_size; i++) {
#pragma omp parallel reduction(custom : result)
      {
        if (omp_get_num_teams() > max_active_teams)
          Kokkos::abort("`omp_set_num_teams` call was not respected.\n");

        typename PolicyType::member_type team(i, league_size, team_size,
                                              vector_length, scratch_ptr, i,
                                              shmem_size_L0, shmem_size_L1);
        if constexpr (std::is_void_v<TagType>)
          f(team, result);
        else
          f(TagType(), team, result);
      }
    }
#endif

    // Copy results back to device if `parallel_reduce` is on a device view.
    ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                 ptr_on_device);
  }

  template <int NumReductions>
  static void execute_array(const FunctorType& f, const PolicyType& p,
                            PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:array_reduction");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:array_reduction");

    const int league_size   = p.league_size();
    const int team_size     = p.team_size();
    const int vector_length = p.impl_vector_length();

    const size_t shmem_size_L0 = p.scratch_size(0, team_size);
    const size_t shmem_size_L1 = p.scratch_size(1, team_size);
    OpenMPTargetExec::resize_scratch(PolicyType::member_type::TEAM_REDUCE_SIZE,
                                     shmem_size_L0, shmem_size_L1, league_size);
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

    // Maximum active teams possible.
    // FIXME_OPENMPTARGET: Cray compiler did not yet implement
    // omp_get_max_teams.
#if !defined(KOKKOS_COMPILER_CRAY_LLVM)
    int max_active_teams = omp_get_max_teams();
#else
    int max_active_teams =
        std::min(OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size, league_size);
#endif

    // If the league size is <=0, do not launch the kernel.
    if (max_active_teams <= 0) return;

    // Case where the number of reduction items is 1.
    if constexpr (NumReductions == 1) {
      ValueType result = ValueType();

      // Case where reduction is on a native data type.
      if constexpr (std::is_arithmetic<ValueType>::value) {
        // Use scratch memory extensions to request dynamic shared memory for
        // the right compiler/architecture combination.
        KOKKOS_IMPL_OMPTARGET_PRAGMA(teams num_teams(max_active_teams) thread_limit(team_size) map(to: f) \
    is_device_ptr(scratch_ptr) reduction(+: result)               \
        KOKKOS_IMPL_OMPX_DYN_CGROUP_MEM(shmem_size_L0))
#pragma omp parallel reduction(+ : result)
        {
          if (omp_get_num_teams() > max_active_teams)
            Kokkos::abort("`omp_set_num_teams` call was not respected.\n");

          const int blockIdx = omp_get_team_num();
          const int gridDim  = omp_get_num_teams();

          // Guarantee that the compilers respect the `num_teams` clause
          for (int league_id = blockIdx; league_id < league_size;
               league_id += gridDim) {
            typename PolicyType::member_type team(
                league_id, league_size, team_size, vector_length, scratch_ptr,
                blockIdx, shmem_size_L0, shmem_size_L1);
            if constexpr (std::is_void_v<TagType>)
              f(team, result);
            else
              f(TagType(), team, result);
          }
        }
      } else {
        // Case where the reduction is on a non-native data type.
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)
#pragma omp target teams num_teams(max_active_teams) thread_limit(team_size) \
    map(to                                                                   \
        : f) is_device_ptr(scratch_ptr) reduction(custom                     \
                                                  : result)
#pragma omp parallel reduction(custom : result)
        {
          if (omp_get_num_teams() > max_active_teams)
            Kokkos::abort("`omp_set_num_teams` call was not respected.\n");

          const int blockIdx = omp_get_team_num();
          const int gridDim  = omp_get_num_teams();

          // Guarantee that the compilers respect the `num_teams` clause
          for (int league_id = blockIdx; league_id < league_size;
               league_id += gridDim) {
            typename PolicyType::member_type team(
                league_id, league_size, team_size, vector_length, scratch_ptr,
                blockIdx, shmem_size_L0, shmem_size_L1);
            if constexpr (std::is_void_v<TagType>)
              f(team, result);
            else
              f(TagType(), team, result);
          }
        }
      }

      // Copy results back to device if `parallel_reduce` is on a device view.
      ParReduceCopy::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                   ptr_on_device);
    } else {
      ValueType result[NumReductions] = {};
      // Case where the reduction is on an array.
#pragma omp target teams num_teams(max_active_teams) thread_limit(team_size) map(to   \
                                                                       : f) \
    is_device_ptr(scratch_ptr) reduction(+ : result[:NumReductions])
#pragma omp parallel reduction(+ : result[:NumReductions])
      {
        if (omp_get_num_teams() > max_active_teams)
          Kokkos::abort("`omp_set_num_teams` call was not respected.\n");

        const int blockIdx = omp_get_team_num();
        const int gridDim  = omp_get_num_teams();

        // Guarantee that the compilers respect the `num_teams` clause
        for (int league_id = blockIdx; league_id < league_size;
             league_id += gridDim) {
          typename PolicyType::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              blockIdx, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_void_v<TagType>)
            f(team, result);
          else
            f(TagType(), team, result);
        }
      }

      // Copy results back to device if `parallel_reduce` is on a device view.
      ParReduceCopy::memcpy_result(
          result_ptr, result, NumReductions * sizeof(ValueType), ptr_on_device);
    }
  }

  // FIXME_OPENMPTARGET : This routine is a copy from `parallel_reduce` over
  // RangePolicy. Need a new implementation.
  static void execute_init_join(const FunctorType& f, const PolicyType& p,
                                PointerType ptr, const bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:init_join ");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget TeamPolicy "
        "parallel_reduce:init_join");
    using FunctorAnalysis =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                              FunctorType, ValueType>;

    const int league_size   = p.league_size();
    const int team_size     = p.team_size();
    const int vector_length = p.impl_vector_length();

    auto begin = 0;
    auto end   = league_size * team_size + team_size * vector_length;

    const size_t shmem_size_L0 = p.scratch_size(0, team_size);
    const size_t shmem_size_L1 = p.scratch_size(1, team_size);

    // FIXME_OPENMPTARGET: This would oversubscribe scratch memory since we are
    // already using the available scratch memory to create temporaries for each
    // thread.
    if ((shmem_size_L0 + shmem_size_L1) > 0) {
      Kokkos::abort(
          "OpenMPTarget: Scratch memory is not supported in `parallel_reduce` "
          "over functors with init/join.");
    }

    const auto nteams = league_size;

    // Number of elements in the reduction
    const auto value_count = FunctorAnalysis::value_count(f);

    // Allocate scratch per active thread.
    OpenMPTargetExec::resize_scratch(1, 0, value_count * sizeof(ValueType),
                                     league_size);
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();
    typename FunctorAnalysis::Reducer final_reducer(f);

    if (end <= begin) {
// If there is no work to be done, copy back the initialized values and
// exit.
#pragma omp target map(to : final_reducer) is_device_ptr(scratch_ptr)
      {
        final_reducer.init(scratch_ptr);
        final_reducer.final(scratch_ptr);
      }

      if (0 < value_count) {
        if (!ptr_on_device)
          KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
              ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
              omp_get_initial_device(), omp_get_default_device()));
        else
          KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
              ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
              omp_get_default_device(), omp_get_default_device()));
      }

      return;
    }
    // Use scratch memory extensions to request dynamic shared memory for the
    // right compiler/architecture combination.
    KOKKOS_IMPL_OMPTARGET_PRAGMA(
        teams num_teams(nteams) thread_limit(team_size) map(to
                                                            : f)
            is_device_ptr(scratch_ptr)
                KOKKOS_IMPL_OMPX_DYN_CGROUP_MEM(shmem_size_L0)) {
#pragma omp parallel
      {
        const int team_num      = omp_get_team_num();
        const int num_teams     = omp_get_num_teams();
        ValueType* team_scratch = static_cast<ValueType*>(scratch_ptr) +
                                  team_num * team_size * value_count;
        ReferenceType result = final_reducer.init(&team_scratch[0]);

        for (int league_id = team_num; league_id < league_size;
             league_id += num_teams) {
          typename PolicyType::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              team_num, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_void_v<TagType>) {
            f(team, result);
          } else {
            f(TagType(), team, result);
          }
        }
      }  // end parallel
    }    // end target

    int tree_neighbor_offset = 1;
    do {
#pragma omp target teams distribute parallel for simd firstprivate( \
    final_reducer) is_device_ptr(scratch_ptr)
      for (int i = 0; i < nteams - tree_neighbor_offset;
           i += 2 * tree_neighbor_offset) {
        ValueType* team_scratch = static_cast<ValueType*>(scratch_ptr);
        const int team_offset   = team_size * value_count;
        final_reducer.join(
            &team_scratch[i * team_offset],
            &team_scratch[(i + tree_neighbor_offset) * team_offset]);

        // If `final` is provided by the functor.
        // Do the final only once at the end.
        if (tree_neighbor_offset * 2 >= nteams && omp_get_team_num() == 0 &&
            omp_get_thread_num() == 0) {
          final_reducer.final(scratch_ptr);
        }
      }
      tree_neighbor_offset *= 2;
    } while (tree_neighbor_offset < nteams);

    // If the result view is on the host, copy back the values via memcpy.
    if (0 < value_count) {
      if (!ptr_on_device)
        KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_initial_device(), omp_get_default_device()));
      else
        KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_default_device(), omp_get_default_device()));
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
