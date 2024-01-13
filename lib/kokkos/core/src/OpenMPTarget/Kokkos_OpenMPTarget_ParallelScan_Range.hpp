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

#ifndef KOKKOS_OPENMPTARGET_PARALLELSCAN_RANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLELSCAN_RANGE_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::OpenMPTarget> {
 protected:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using WorkTag  = typename Policy::work_tag;
  using Member   = typename Policy::member_type;
  using idx_type = typename Policy::index_type;

  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         Policy, FunctorType, void>;

  using value_type     = typename Analysis::value_type;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;

  value_type* m_result_ptr;
  const bool m_result_ptr_device_accessible;

  template <class TagType>
  std::enable_if_t<std::is_void<TagType>::value> call_with_tag(
      const FunctorType& f, const idx_type& idx, value_type& val,
      const bool& is_final) const {
    f(idx, val, is_final);
  }
  template <class TagType>
  std::enable_if_t<!std::is_void<TagType>::value> call_with_tag(
      const FunctorType& f, const idx_type& idx, value_type& val,
      const bool& is_final) const {
    f(WorkTag(), idx, val, is_final);
  }

 public:
  void impl_execute(
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

    auto a_functor_reducer = m_functor_reducer;
#pragma omp target teams distribute map(to \
                                        : a_functor_reducer) num_teams(nteams)
    for (idx_type team_id = 0; team_id < n_chunks; ++team_id) {
      const typename Analysis::Reducer& reducer =
          a_functor_reducer.get_reducer();
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;

#pragma omp for
        for (idx_type i = 0; i < chunk_size; ++i) {
          const idx_type idx = local_offset + i;
          value_type val;
          reducer.init(&val);
          if (idx < N)
            call_with_tag<WorkTag>(a_functor_reducer.get_functor(), idx, val,
                                   false);
          element_values(team_id, i) = val;
        }
#pragma omp barrier
        if (omp_get_thread_num() == 0) {
          value_type sum;
          reducer.init(&sum);
          for (idx_type i = 0; i < chunk_size; ++i) {
            reducer.join(&sum, &element_values(team_id, i));
            element_values(team_id, i) = sum;
          }
          chunk_values(team_id) = sum;
        }
#pragma omp barrier
        if (omp_get_thread_num() == 0) {
          if (Kokkos::atomic_fetch_add(&count(), 1) == n_chunks - 1) {
            value_type sum;
            reducer.init(&sum);
            for (idx_type i = 0; i < n_chunks; ++i) {
              reducer.join(&sum, &chunk_values(i));
              chunk_values(i) = sum;
            }
          }
        }
      }
    }

#pragma omp target teams distribute map(to                                     \
                                        : a_functor_reducer) num_teams(nteams) \
    thread_limit(team_size)
    for (idx_type team_id = 0; team_id < n_chunks; ++team_id) {
      const typename Analysis::Reducer& reducer =
          a_functor_reducer.get_reducer();
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;
        value_type offset_value;
        if (team_id > 0)
          offset_value = chunk_values(team_id - 1);
        else
          reducer.init(&offset_value);

#pragma omp for
        for (idx_type i = 0; i < chunk_size; ++i) {
          const idx_type idx = local_offset + i;
          value_type local_offset_value;
          if (i > 0) {
            local_offset_value = element_values(team_id, i - 1);
            // FIXME_OPENMPTARGET We seem to access memory illegaly on AMD GPUs
#if defined(KOKKOS_ARCH_AMD_GPU) && !defined(KOKKOS_ARCH_AMD_GFX1030) && \
    !defined(KOKKOS_ARCH_AMD_GFX1100)
            if constexpr (Analysis::Reducer::has_join_member_function()) {
              if constexpr (std::is_void_v<WorkTag>)
                a_functor_reducer.get_functor().join(local_offset_value,
                                                     offset_value);
              else
                a_functor_reducer.get_functor().join(
                    WorkTag{}, local_offset_value, offset_value);
            } else
              local_offset_value += offset_value;
#else
            reducer.join(&local_offset_value, &offset_value);
#endif
          } else
            local_offset_value = offset_value;
          if (idx < N)
            call_with_tag<WorkTag>(a_functor_reducer.get_functor(), idx,
                                   local_offset_value, true);
          if (idx == N - 1 && m_result_ptr_device_accessible)
            *m_result_ptr = local_offset_value;
        }
      }
    }
  }

  void execute() const {
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

  ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy,
               pointer_type arg_result_ptr           = nullptr,
               bool arg_result_ptr_device_accessible = false)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy),
        m_result_ptr(arg_result_ptr),
        m_result_ptr_device_accessible(arg_result_ptr_device_accessible) {}

  //----------------------------------------
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::OpenMPTarget>
    : public ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                          Kokkos::Experimental::OpenMPTarget> {
  using base_t     = ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                              Kokkos::Experimental::OpenMPTarget>;
  using value_type = typename base_t::value_type;

 public:
  void execute() const {
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

      if (!base_t::m_result_ptr_device_accessible) {
        const int size = base_t::m_functor_reducer.get_reducer().value_size();
        DeepCopy<HostSpace, Kokkos::Experimental::OpenMPTargetSpace>(
            base_t::m_result_ptr, chunk_values.data() + (n_chunks - 1), size);
      }
    } else if (!base_t::m_result_ptr_device_accessible) {
      base_t::m_functor_reducer.get_reducer().init(base_t::m_result_ptr);
    }
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const typename base_t::Policy& arg_policy,
                        const ViewType& arg_result_view)
      : base_t(arg_functor, arg_policy, arg_result_view.data(),
               MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                                 typename ViewType::memory_space>::accessible) {
  }
};
}  // namespace Impl
}  // namespace Kokkos

#endif
