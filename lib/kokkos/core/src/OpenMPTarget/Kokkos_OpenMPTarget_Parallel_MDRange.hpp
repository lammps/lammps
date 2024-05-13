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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_MDRANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_MDRANGE_HPP

#include <omp.h>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel_Common.hpp>

// WORKAROUND OPENMPTARGET: sometimes tile sizes don't make it correctly,
// this was tracked down to a bug in clang with regards of mapping structs
// with arrays of long in it. Arrays of int might be fine though ...
#define KOKKOS_IMPL_MDRANGE_USE_NO_TILES  // undef EOF

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy  = Kokkos::MDRangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;
  using Index   = typename Policy::index_type;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  inline void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    FunctorType functor(m_functor);
    Policy policy = m_policy;

#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    typename Policy::point_type unused;

    execute_tile<Policy::rank>(unused, functor, policy);
#else
    const int64_t begin = 0;
    const int64_t end   = m_policy.m_num_tiles;

#pragma omp target teams distribute map(to : functor) num_teams(end - begin)
    {
      for (ptrdiff_t tile_idx = begin; tile_idx < end; ++tile_idx) {

#pragma omp parallel
        {
          typename Policy::point_type offset;
          if (Policy::outer_direction == Policy::Left) {
            for (int i = 0; i < Policy::rank; ++i) {
              offset[i] = (tile_idx % policy.m_tile_end[i]) * policy.m_tile[i] +
                          policy.m_lower[i];
              tile_idx /= policy.m_tile_end[i];
            }
          } else {
            for (int i = Policy::rank - 1; i >= 0; --i) {
              offset[i] = (tile_idx % policy.m_tile_end[i]) * policy.m_tile[i] +
                          policy.m_lower[i];
              tile_idx /= policy.m_tile_end[i];
            }
          }
          execute_tile<Policy::rank>(offset, functor, policy);
        }
      }
    }
#endif
  }

  template <int Rank>
  inline std::enable_if_t<Rank == 2> execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    (void)offset;
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];

#pragma omp target teams distribute parallel for collapse(2) map(to : functor)
    for (auto i0 = begin_0; i0 < end_0; ++i0) {
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        if constexpr (std::is_void<typename Policy::work_tag>::value)
          functor(i0, i1);
        else
          functor(typename Policy::work_tag(), i0, i1);
      }
    }
#else
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

#pragma omp for collapse(2)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; ++i0)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; ++i1) {
        if constexpr (std::is_void<typename Policy::work_tag>::value)
          functor(i0, i1);
        else
          functor(typename Policy::work_tag(), i0, i1);
      }
#endif
  }

  template <int Rank>
  inline std::enable_if_t<Rank == 3> execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    (void)offset;
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];

#pragma omp target teams distribute parallel for collapse(3) map(to : functor)
    for (auto i0 = begin_0; i0 < end_0; ++i0) {
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          if constexpr (std::is_void<typename Policy::work_tag>::value)
            functor(i0, i1, i2);
          else
            functor(typename Policy::work_tag(), i0, i1, i2);
        }
      }
    }
#else
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

    const ptrdiff_t begin_2 = offset[2];
    ptrdiff_t end_2         = begin_2 + policy.m_tile[2];
    end_2 = end_2 < policy.m_upper[2] ? end_2 : policy.m_upper[2];

#pragma omp for collapse(3)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; ++i0)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; ++i1)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; ++i2) {
          if constexpr (std::is_void<typename Policy::work_tag>::value)
            functor(i0, i1, i2);
          else
            functor(typename Policy::work_tag(), i0, i1, i2);
        }
#endif
  }

  template <int Rank>
  inline std::enable_if_t<Rank == 4> execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    (void)offset;
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];
    const Index begin_3 = policy.m_lower[3];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];

#pragma omp target teams distribute parallel for collapse(4) map(to : functor)
    for (auto i0 = begin_0; i0 < end_0; ++i0) {
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          for (auto i3 = begin_3; i3 < end_3; ++i3) {
            if constexpr (std::is_void<typename Policy::work_tag>::value)
              functor(i0, i1, i2, i3);
            else
              functor(typename Policy::work_tag(), i0, i1, i2, i3);
          }
        }
      }
    }
#else
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

    const ptrdiff_t begin_2 = offset[2];
    ptrdiff_t end_2         = begin_2 + policy.m_tile[2];
    end_2 = end_2 < policy.m_upper[2] ? end_2 : policy.m_upper[2];

    const ptrdiff_t begin_3 = offset[3];
    ptrdiff_t end_3         = begin_3 + policy.m_tile[3];
    end_3 = end_3 < policy.m_upper[3] ? end_3 : policy.m_upper[3];

#pragma omp for collapse(4)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; ++i0)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; ++i1)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; ++i2)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; ++i3) {
            if constexpr (std::is_void<typename Policy::work_tag>::value)
              functor(i0, i1, i2, i3);
            else
              functor(typename Policy::work_tag(), i0, i1, i2, i3);
          }
#endif
  }

  template <int Rank>
  inline std::enable_if_t<Rank == 5> execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    (void)offset;
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];
    const Index begin_3 = policy.m_lower[3];
    const Index begin_4 = policy.m_lower[4];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];
    const Index end_4 = policy.m_upper[4];

#pragma omp target teams distribute parallel for collapse(5) map(to : functor)
    for (auto i0 = begin_0; i0 < end_0; ++i0) {
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          for (auto i3 = begin_3; i3 < end_3; ++i3) {
            for (auto i4 = begin_4; i4 < end_4; ++i4) {
              if constexpr (std::is_same<typename Policy::work_tag,
                                         void>::value)
                functor(i0, i1, i2, i3, i4);
              else
                functor(typename Policy::work_tag(), i0, i1, i2, i3, i4);
            }
          }
        }
      }
    }
#else
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

    const ptrdiff_t begin_2 = offset[2];
    ptrdiff_t end_2         = begin_2 + policy.m_tile[2];
    end_2 = end_2 < policy.m_upper[2] ? end_2 : policy.m_upper[2];

    const ptrdiff_t begin_3 = offset[3];
    ptrdiff_t end_3         = begin_3 + policy.m_tile[3];
    end_3 = end_3 < policy.m_upper[3] ? end_3 : policy.m_upper[3];

    const ptrdiff_t begin_4 = offset[4];
    ptrdiff_t end_4         = begin_4 + policy.m_tile[4];
    end_4 = end_4 < policy.m_upper[4] ? end_4 : policy.m_upper[4];

#pragma omp for collapse(5)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; ++i0)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; ++i1)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; ++i2)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; ++i3)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; ++i4) {
              if constexpr (std::is_same<typename Policy::work_tag,
                                         void>::value)
                functor(i0, i1, i2, i3, i4);
              else
                functor(typename Policy::work_tag(), i0, i1, i2, i3, i4);
            }
#endif
  }

  template <int Rank>
  inline std::enable_if_t<Rank == 6> execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
#ifdef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
    (void)offset;
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];
    const Index begin_3 = policy.m_lower[3];
    const Index begin_4 = policy.m_lower[4];
    const Index begin_5 = policy.m_lower[5];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];
    const Index end_4 = policy.m_upper[4];
    const Index end_5 = policy.m_upper[5];

#pragma omp target teams distribute parallel for collapse(6) map(to : functor)
    for (auto i0 = begin_0; i0 < end_0; ++i0) {
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          for (auto i3 = begin_3; i3 < end_3; ++i3) {
            for (auto i4 = begin_4; i4 < end_4; ++i4) {
              for (auto i5 = begin_5; i5 < end_5; ++i5) {
                {
                  if constexpr (std::is_same<typename Policy::work_tag,
                                             void>::value)
                    functor(i0, i1, i2, i3, i4, i5);
                  else
                    functor(typename Policy::work_tag(), i0, i1, i2, i3, i4,
                            i5);
                }
              }
            }
          }
        }
      }
    }
#else
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

    const ptrdiff_t begin_2 = offset[2];
    ptrdiff_t end_2         = begin_2 + policy.m_tile[2];
    end_2 = end_2 < policy.m_upper[2] ? end_2 : policy.m_upper[2];

    const ptrdiff_t begin_3 = offset[3];
    ptrdiff_t end_3         = begin_3 + policy.m_tile[3];
    end_3 = end_3 < policy.m_upper[3] ? end_3 : policy.m_upper[3];

    const ptrdiff_t begin_4 = offset[4];
    ptrdiff_t end_4         = begin_4 + policy.m_tile[4];
    end_4 = end_4 < policy.m_upper[4] ? end_4 : policy.m_upper[4];

    const ptrdiff_t begin_5 = offset[5];
    ptrdiff_t end_5         = begin_5 + policy.m_tile[5];
    end_5 = end_5 < policy.m_upper[5] ? end_5 : policy.m_upper[5];

#pragma omp for collapse(6)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; ++i0)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; ++i1)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; ++i2)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; ++i3)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; ++i4)
              for (ptrdiff_t i5 = begin_5; i5 < end_5; ++i5) {
                if constexpr (std::is_same<typename Policy::work_tag,
                                           void>::value)
                  functor(i0, i1, i2, i3, i4, i5);
                else
                  functor(typename Policy::work_tag(), i0, i1, i2, i3, i4, i5);
              }
#endif
  }

  inline ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
  // TODO DZP: based on a conversation with Christian, we're using 256 as a
  // heuristic here. We need something better once we can query these kinds of
  // properties
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    return 256;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy      = Kokkos::MDRangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;
  using Index   = typename Policy::index_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  static constexpr bool UseReducer =
      !std::is_same_v<FunctorType, typename ReducerType::functor_type>;

  const pointer_type m_result_ptr;
  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;

  using ParReduceCopy = ParallelReduceCopy<pointer_type>;

  bool m_result_ptr_on_device;

  // Only let one ParallelReduce instance at a time use the scratch memory.
  // The constructor acquires the mutex which is released in the destructor.
  std::scoped_lock<std::mutex> m_scratch_memory_lock;

 public:
  inline void execute() const {
    execute_tile<Policy::rank, typename ReducerType::value_type>(
        m_functor_reducer.get_functor(), m_policy, m_result_ptr);
  }

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                        Policy arg_policy, const ViewType& arg_result_view)
      : m_result_ptr(arg_result_view.data()),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible),
        m_scratch_memory_lock(OpenMPTargetExec::m_mutex_scratch_ptr) {}

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 2> execute_tile(const FunctorType& functor,
                                                  const Policy& policy,
                                                  pointer_type ptr) const {
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];

    ValueType result = ValueType();

    // FIXME_OPENMPTARGET: Unable to separate directives and their companion
    // loops which leads to code duplication for different reduction types.
    if constexpr (UseReducer) {
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp target teams distribute parallel for collapse(2) map(to         \
                                                                 : functor) \
    reduction(custom                                                        \
              : result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          if constexpr (std::is_void<typename Policy::work_tag>::value)
            functor(i0, i1, result);
          else
            functor(typename Policy::work_tag(), i0, i1, result);
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(2) map(to : functor) \
reduction(+:result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          if constexpr (std::is_void<typename Policy::work_tag>::value)
            functor(i0, i1, result);
          else
            functor(typename Policy::work_tag(), i0, i1, result);
        }
      }
    }

    ParReduceCopy::memcpy_result(ptr, &result, sizeof(ValueType),
                                 m_result_ptr_on_device);
  }

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 3> execute_tile(const FunctorType& functor,
                                                  const Policy& policy,
                                                  pointer_type ptr) const {
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];

    ValueType result = ValueType();

    // FIXME_OPENMPTARGET: Unable to separate directives and their companion
    // loops which leads to code duplication for different reduction types.
    if constexpr (UseReducer) {
#pragma omp declare reduction(                                                 \
    custom:ValueType                                                           \
    : OpenMPTargetReducerWrapper <typename ReducerType::functor_type>::join(   \
        omp_out, omp_in))                                                      \
    initializer(                                                               \
        OpenMPTargetReducerWrapper <typename ReducerType::functor_type>::init( \
            omp_priv))

#pragma omp target teams distribute parallel for collapse(3) map(to         \
                                                                 : functor) \
    reduction(custom                                                        \
              : result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            if constexpr (std::is_void<typename Policy::work_tag>::value)
              functor(i0, i1, i2, result);
            else
              functor(typename Policy::work_tag(), i0, i1, i2, result);
          }
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(3) map(to : functor) \
reduction(+:result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            if constexpr (std::is_void<typename Policy::work_tag>::value)
              functor(i0, i1, i2, result);
            else
              functor(typename Policy::work_tag(), i0, i1, i2, result);
          }
        }
      }
    }

    ParReduceCopy::memcpy_result(ptr, &result, sizeof(ValueType),
                                 m_result_ptr_on_device);
  }

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 4> execute_tile(const FunctorType& functor,
                                                  const Policy& policy,
                                                  pointer_type ptr) const {
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[3];
    const Index begin_3 = policy.m_lower[2];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];

    ValueType result = ValueType();

    // FIXME_OPENMPTARGET: Unable to separate directives and their companion
    // loops which leads to code duplication for different reduction types.
    if constexpr (UseReducer) {
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp target teams distribute parallel for collapse(4) map(to         \
                                                                 : functor) \
    reduction(custom                                                        \
              : result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              if constexpr (std::is_same<typename Policy::work_tag,
                                         void>::value)
                functor(i0, i1, i2, i3, result);
              else
                functor(typename Policy::work_tag(), i0, i1, i2, i3, result);
            }
          }
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(4) map(to : functor) \
reduction(+:result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              if constexpr (std::is_same<typename Policy::work_tag,
                                         void>::value)
                functor(i0, i1, i2, i3, result);
              else
                functor(typename Policy::work_tag(), i0, i1, i2, i3, result);
            }
          }
        }
      }
    }

    ParReduceCopy::memcpy_result(ptr, &result, sizeof(ValueType),
                                 m_result_ptr_on_device);
  }

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 5> execute_tile(const FunctorType& functor,
                                                  const Policy& policy,
                                                  pointer_type ptr) const {
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];
    const Index begin_3 = policy.m_lower[3];
    const Index begin_4 = policy.m_lower[4];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];
    const Index end_4 = policy.m_upper[4];

    ValueType result = ValueType();

    // FIXME_OPENMPTARGET: Unable to separate directives and their companion
    // loops which leads to code duplication for different reduction types.
    if constexpr (UseReducer) {
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp target teams distribute parallel for collapse(5) map(to         \
                                                                 : functor) \
    reduction(custom                                                        \
              : result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              for (auto i4 = begin_4; i4 < end_4; ++i4) {
                if constexpr (std::is_same<typename Policy::work_tag,
                                           void>::value)
                  functor(i0, i1, i2, i3, i4, result);
                else
                  functor(typename Policy::work_tag(), i0, i1, i2, i3, i4,
                          result);
              }
            }
          }
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(5) map(to : functor) \
reduction(+:result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              for (auto i4 = begin_4; i4 < end_4; ++i4) {
                if constexpr (std::is_same<typename Policy::work_tag,
                                           void>::value)
                  functor(i0, i1, i2, i3, i4, result);
                else
                  functor(typename Policy::work_tag(), i0, i1, i2, i3, i4,
                          result);
              }
            }
          }
        }
      }
    }

    ParReduceCopy::memcpy_result(ptr, &result, sizeof(ValueType),
                                 m_result_ptr_on_device);
  }

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 6> execute_tile(const FunctorType& functor,
                                                  const Policy& policy,
                                                  pointer_type ptr) const {
    const Index begin_0 = policy.m_lower[0];
    const Index begin_1 = policy.m_lower[1];
    const Index begin_2 = policy.m_lower[2];
    const Index begin_3 = policy.m_lower[3];
    const Index begin_4 = policy.m_lower[4];
    const Index begin_5 = policy.m_lower[5];

    const Index end_0 = policy.m_upper[0];
    const Index end_1 = policy.m_upper[1];
    const Index end_2 = policy.m_upper[2];
    const Index end_3 = policy.m_upper[3];
    const Index end_4 = policy.m_upper[4];
    const Index end_5 = policy.m_upper[5];

    ValueType result = ValueType();

    // FIXME_OPENMPTARGET: Unable to separate directives and their companion
    // loops which leads to code duplication for different reduction types.
    if constexpr (UseReducer) {
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp target teams distribute parallel for collapse(6) map(to         \
                                                                 : functor) \
    reduction(custom                                                        \
              : result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              for (auto i4 = begin_4; i4 < end_4; ++i4) {
                for (auto i5 = begin_5; i5 < end_5; ++i5) {
                  if constexpr (std::is_same<typename Policy::work_tag,
                                             void>::value)
                    functor(i0, i1, i2, i3, i4, i5, result);
                  else
                    functor(typename Policy::work_tag(), i0, i1, i2, i3, i4, i5,
                            result);
                }
              }
            }
          }
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(6) map(to : functor) \
reduction(+:result)
      for (auto i0 = begin_0; i0 < end_0; ++i0) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i3 = begin_3; i3 < end_3; ++i3) {
              for (auto i4 = begin_4; i4 < end_4; ++i4) {
                for (auto i5 = begin_5; i5 < end_5; ++i5) {
                  if constexpr (std::is_same<typename Policy::work_tag,
                                             void>::value)
                    functor(i0, i1, i2, i3, i4, i5, result);
                  else
                    functor(typename Policy::work_tag(), i0, i1, i2, i3, i4, i5,
                            result);
                }
              }
            }
          }
        }
      }
    }

    ParReduceCopy::memcpy_result(ptr, &result, sizeof(ValueType),
                                 m_result_ptr_on_device);
  }

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    return 256;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#undef KOKKOS_IMPL_MDRANGE_USE_NO_TILES
#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */
