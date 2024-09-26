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

#ifndef KOKKOS_OPENMPTARGET_PARALLELREDUCE_MDRANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLELREDUCE_MDRANGE_HPP

#include <omp.h>
#include <Kokkos_Parallel.hpp>
#include "Kokkos_OpenMPTarget_MDRangePolicy.hpp"
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel_Common.hpp>

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

 public:
  inline void execute() const {
    // Only let one ParallelReduce instance at a time use the scratch memory.
    std::scoped_lock<std::mutex> scratch_memory_lock(
        OpenMPTargetExec::m_mutex_scratch_ptr);
    execute_tile<Policy::rank, typename ReducerType::value_type>(
        m_functor_reducer.get_functor(), m_policy, m_result_ptr,
        std::integral_constant<Iterate, Policy::inner_direction>());
  }

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                        Policy arg_policy, const ViewType& arg_result_view)
      : m_result_ptr(arg_result_view.data()),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible) {}

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 2> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateLeft) const {
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
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i0 = begin_0; i0 < end_0; ++i0) {
          if constexpr (std::is_void<typename Policy::work_tag>::value)
            functor(i0, i1, result);
          else
            functor(typename Policy::work_tag(), i0, i1, result);
        }
      }
    } else {
#pragma omp target teams distribute parallel for collapse(2) map(to : functor) \
reduction(+:result)
      for (auto i1 = begin_1; i1 < end_1; ++i1) {
        for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
  inline std::enable_if_t<Rank == 3> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateLeft) const {
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
      for (auto i2 = begin_2; i2 < end_2; ++i2) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
      for (auto i2 = begin_2; i2 < end_2; ++i2) {
        for (auto i1 = begin_1; i1 < end_1; ++i1) {
          for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
  inline std::enable_if_t<Rank == 4> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateLeft) const {
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
      for (auto i3 = begin_3; i3 < end_3; ++i3) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          for (auto i1 = begin_1; i1 < end_1; ++i1) {
            for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
      for (auto i3 = begin_3; i3 < end_3; ++i3) {
        for (auto i2 = begin_2; i2 < end_2; ++i2) {
          for (auto i1 = begin_1; i1 < end_1; ++i1) {
            for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
  inline std::enable_if_t<Rank == 5> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateLeft) const {
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
      for (auto i4 = begin_4; i4 < end_4; ++i4) {
        for (auto i3 = begin_3; i3 < end_3; ++i3) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i1 = begin_1; i1 < end_1; ++i1) {
              for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
      for (auto i4 = begin_4; i4 < end_4; ++i4) {
        for (auto i3 = begin_3; i3 < end_3; ++i3) {
          for (auto i2 = begin_2; i2 < end_2; ++i2) {
            for (auto i1 = begin_1; i1 < end_1; ++i1) {
              for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
  inline std::enable_if_t<Rank == 6> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateLeft) const {
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
      for (auto i5 = begin_5; i5 < end_5; ++i5) {
        for (auto i4 = begin_4; i4 < end_4; ++i4) {
          for (auto i3 = begin_3; i3 < end_3; ++i3) {
            for (auto i2 = begin_2; i2 < end_2; ++i2) {
              for (auto i1 = begin_1; i1 < end_1; ++i1) {
                for (auto i0 = begin_0; i0 < end_0; ++i0) {
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
      for (auto i5 = begin_5; i5 < end_5; ++i5) {
        for (auto i4 = begin_4; i4 < end_4; ++i4) {
          for (auto i3 = begin_3; i3 < end_3; ++i3) {
            for (auto i2 = begin_2; i2 < end_2; ++i2) {
              for (auto i1 = begin_1; i1 < end_1; ++i1) {
                for (auto i0 = begin_0; i0 < end_0; ++i0) {
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

  template <int Rank, class ValueType>
  inline std::enable_if_t<Rank == 2> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateRight) const {
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
  inline std::enable_if_t<Rank == 3> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateRight) const {
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
  inline std::enable_if_t<Rank == 4> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateRight) const {
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
  inline std::enable_if_t<Rank == 5> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateRight) const {
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
  inline std::enable_if_t<Rank == 6> execute_tile(
      const FunctorType& functor, const Policy& policy, pointer_type ptr,
      OpenMPTargetIterateRight) const {
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
#endif /* KOKKOS_OPENMPTARGET_PARALLELREDUCE_MDRANGE_HPP */
