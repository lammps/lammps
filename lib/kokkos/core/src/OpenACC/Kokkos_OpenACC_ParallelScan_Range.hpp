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

#ifndef KOKKOS_OPENACC_PARALLEL_SCAN_RANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_SCAN_RANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <Kokkos_Parallel.hpp>

namespace Kokkos::Impl {

template <class Functor, class GivenValueType, class... Traits>
class ParallelScanOpenACCBase {
 protected:
  using Policy = Kokkos::RangePolicy<Traits...>;
  using Analysis =
      Kokkos::Impl::FunctorAnalysis<Kokkos::Impl::FunctorPatternInterface::SCAN,
                                    Policy, Functor, GivenValueType>;
  using PointerType = typename Analysis::pointer_type;
  using ValueType   = typename Analysis::value_type;
  using MemberType  = typename Policy::member_type;
  using IndexType   = typename Policy::index_type;
  Functor m_functor;
  Policy m_policy;
  ValueType* m_result_ptr;
  bool m_result_ptr_device_accessible;
  static constexpr MemberType default_scan_chunk_size = 128;

 public:
  ParallelScanOpenACCBase(Functor const& arg_functor, Policy const& arg_policy,
                          ValueType* arg_result_ptr,
                          bool arg_result_ptr_device_accessible)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_result_ptr(arg_result_ptr),
        m_result_ptr_device_accessible(arg_result_ptr_device_accessible) {}

  // This function implements the parallel scan alogithm based on the parallel
  // prefix sum algorithm proposed by Hillis and Steele (doi:10.1145/7902.7903),
  // which offers a shorter span and more parallelism but may not be
  // work-efficient.
  void OpenACCParallelScanRangePolicy(const IndexType begin,
                                      const IndexType end, IndexType chunk_size,
                                      const int async_arg) const {
    if (chunk_size > 1) {
      if (!Impl::is_integral_power_of_two(chunk_size))
        Kokkos::abort(
            "RangePolicy blocking granularity must be power of two to be used "
            "with OpenACC parallel_scan()");
    } else {
      chunk_size = default_scan_chunk_size;
    }
    const Kokkos::Experimental::Impl::FunctorAdapter<
        Functor, Policy, Kokkos::Experimental::Impl::RoutineClause::seq>
        functor(m_functor);
    const IndexType N        = end - begin;
    const IndexType n_chunks = (N + chunk_size - 1) / chunk_size;
    Kokkos::View<ValueType*, Kokkos::Experimental::OpenACCSpace> chunk_values(
        "Kokkos::OpenACCParallelScan::chunk_values", n_chunks);
    Kokkos::View<ValueType*, Kokkos::Experimental::OpenACCSpace> offset_values(
        "Kokkos::OpenACCParallelScan::offset_values", n_chunks);
    Kokkos::View<ValueType, Kokkos::Experimental::OpenACCSpace> m_result_total(
        "Kokkos::OpenACCParallelScan::m_result_total");
    std::unique_ptr<ValueType[]> element_values_owner(
        new ValueType[2 * chunk_size]);
    ValueType* element_values = element_values_owner.get();
    typename Analysis::Reducer final_reducer(m_functor);

#pragma acc enter data copyin(functor, final_reducer) \
    copyin(chunk_values, offset_values) async(async_arg)

#pragma acc parallel loop gang vector_length(chunk_size) private( \
    element_values [0:2 * chunk_size])                            \
    present(functor, chunk_values, final_reducer) async(async_arg)
    for (IndexType team_id = 0; team_id < n_chunks; ++team_id) {
      IndexType current_step = 0;
      IndexType next_step    = 1;
      IndexType temp;
#pragma acc loop vector
      for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
        const IndexType local_offset = team_id * chunk_size;
        const IndexType idx          = local_offset + thread_id;
        ValueType update;
        final_reducer.init(&update);
        if ((idx > 0) && (idx < N)) functor(idx - 1, update, false);
        element_values[thread_id] = update;
      }
      for (IndexType step_size = 1; step_size < chunk_size; step_size *= 2) {
#pragma acc loop vector
        for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
          if (thread_id < step_size) {
            element_values[next_step * chunk_size + thread_id] =
                element_values[current_step * chunk_size + thread_id];
          } else {
            ValueType localValue =
                element_values[current_step * chunk_size + thread_id];
            final_reducer.join(&localValue,
                               &element_values[current_step * chunk_size +
                                               thread_id - step_size]);
            element_values[next_step * chunk_size + thread_id] = localValue;
          }
        }
        temp         = current_step;
        current_step = next_step;
        next_step    = temp;
      }
      chunk_values(team_id) =
          element_values[current_step * chunk_size + chunk_size - 1];
    }

    ValueType tempValue;
#pragma acc serial loop present(chunk_values, offset_values, final_reducer) \
    async(async_arg)
    for (IndexType team_id = 0; team_id < n_chunks; ++team_id) {
      if (team_id == 0) {
        final_reducer.init(&offset_values(0));
        final_reducer.init(&tempValue);
      } else {
        final_reducer.join(&tempValue, &chunk_values(team_id - 1));
        offset_values(team_id) = tempValue;
      }
    }

#pragma acc parallel loop gang vector_length(chunk_size) private(         \
    element_values [0:2 * chunk_size])                                    \
    present(functor, offset_values, final_reducer) copyin(m_result_total) \
        async(async_arg)
    for (IndexType team_id = 0; team_id < n_chunks; ++team_id) {
      IndexType current_step = 0;
      IndexType next_step    = 1;
      IndexType temp;
#pragma acc loop vector
      for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
        const IndexType local_offset = team_id * chunk_size;
        const IndexType idx          = local_offset + thread_id;
        ValueType update;
        final_reducer.init(&update);
        if (thread_id == 0) {
          final_reducer.join(&update, &offset_values(team_id));
        }
        if ((idx > 0) && (idx < N)) functor(idx - 1, update, false);
        element_values[thread_id] = update;
      }
      for (IndexType step_size = 1; step_size < chunk_size; step_size *= 2) {
#pragma acc loop vector
        for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
          if (thread_id < step_size) {
            element_values[next_step * chunk_size + thread_id] =
                element_values[current_step * chunk_size + thread_id];
          } else {
            ValueType localValue =
                element_values[current_step * chunk_size + thread_id];
            final_reducer.join(&localValue,
                               &element_values[current_step * chunk_size +
                                               thread_id - step_size]);
            element_values[next_step * chunk_size + thread_id] = localValue;
          }
        }
        temp         = current_step;
        current_step = next_step;
        next_step    = temp;
      }
#pragma acc loop vector
      for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
        const IndexType local_offset = team_id * chunk_size;
        const IndexType idx          = local_offset + thread_id;
        ValueType update =
            element_values[current_step * chunk_size + thread_id];
        if (idx < N) functor(idx, update, true);
        if (idx == N - 1) {
          if (m_result_ptr_device_accessible) {
            *m_result_ptr = update;
          } else {
            m_result_total() = update;
          }
        }
      }
    }
    if (!m_result_ptr_device_accessible && m_result_ptr != nullptr) {
      DeepCopy<HostSpace, Kokkos::Experimental::OpenACCSpace,
               Kokkos::Experimental::OpenACC>(m_policy.space(), m_result_ptr,
                                              m_result_total.data(),
                                              sizeof(ValueType));
    }

#pragma acc exit data delete (functor, chunk_values, offset_values, \
                              final_reducer)async(async_arg)
    acc_wait(async_arg);
  }

  void execute() const {
    const IndexType begin = m_policy.begin();
    const IndexType end   = m_policy.end();
    IndexType chunk_size  = m_policy.chunk_size();

    if (end <= begin) {
      if (!m_result_ptr_device_accessible && m_result_ptr != nullptr) {
        *m_result_ptr = 0;
      }
      return;
    }

    int const async_arg = m_policy.space().acc_async_queue();

    OpenACCParallelScanRangePolicy(begin, end, chunk_size, async_arg);
  }
};

}  // namespace Kokkos::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class Functor, class... Traits>
class Kokkos::Impl::ParallelScan<Functor, Kokkos::RangePolicy<Traits...>,
                                 Kokkos::Experimental::OpenACC>
    : public ParallelScanOpenACCBase<Functor, void, Traits...> {
  using base_t    = ParallelScanOpenACCBase<Functor, void, Traits...>;
  using IndexType = typename base_t::IndexType;

 public:
  void execute() const {
    const IndexType begin = base_t::m_policy.begin();
    const IndexType end   = base_t::m_policy.end();
    IndexType chunk_size  = base_t::m_policy.chunk_size();

    int const async_arg = base_t::m_policy.space().acc_async_queue();

    OpenACCParallelScanRangePolicy(begin, end, chunk_size, async_arg);
  }

  ParallelScan(const Functor& arg_functor,
               const typename base_t::Policy& arg_policy)
      : base_t(arg_functor, arg_policy, nullptr, false) {}
};

template <class FunctorType, class ReturnType, class... Traits>
class Kokkos::Impl::ParallelScanWithTotal<
    FunctorType, Kokkos::RangePolicy<Traits...>, ReturnType,
    Kokkos::Experimental::OpenACC>
    : public ParallelScanOpenACCBase<FunctorType, ReturnType, Traits...> {
  using base_t    = ParallelScanOpenACCBase<FunctorType, ReturnType, Traits...>;
  using IndexType = typename base_t::IndexType;

 public:
  void execute() const {
    const IndexType begin = base_t::m_policy.begin();
    const IndexType end   = base_t::m_policy.end();
    IndexType chunk_size  = base_t::m_policy.chunk_size();

    if (end <= begin) {
      if (!base_t::m_result_ptr_device_accessible &&
          base_t::m_result_ptr != nullptr) {
        *base_t::m_result_ptr = 0;
      }
      return;
    }

    int const async_arg = base_t::m_policy.space().acc_async_queue();

    OpenACCParallelScanRangePolicy(begin, end, chunk_size, async_arg);
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const typename base_t::Policy& arg_policy,
                        const ViewType& arg_result_view)
      : base_t(arg_functor, arg_policy, arg_result_view.data(),
               MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                                 typename ViewType::memory_space>::accessible) {
  }
};

#endif
