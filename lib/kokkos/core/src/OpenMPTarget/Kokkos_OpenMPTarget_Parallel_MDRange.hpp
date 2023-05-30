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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_MDRANGE_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_MDRANGE_HPP

#include <omp.h>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Exec.hpp>

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

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;
  using Index   = typename Policy::index_type;

  using ReducerConditional =
      std::conditional<std::is_same<InvalidType, ReducerType>::value,
                       FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                         Policy, ReducerTypeFwd>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  enum {
    HasJoin =
        Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, Policy,
                              FunctorType>::has_join_member_function
  };
  enum { UseReducer = is_reducer<ReducerType>::value };

  const pointer_type m_result_ptr;
  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;

  using ParReduceCommon = ParallelReduceCommon<pointer_type>;

  bool m_result_ptr_on_device;

 public:
  inline void execute() const {
    execute_tile<Policy::rank, typename Analysis::value_type>(
        m_functor, m_policy, m_result_ptr);
  }

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& arg_functor, Policy arg_policy,
      const ViewType& arg_result_view,
      std::enable_if_t<Kokkos::is_view<ViewType>::value &&
                           !Kokkos::is_reducer<ReducerType>::value,
                       void*> = NULL)
      : m_result_ptr(arg_result_view.data()),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible) {}

  inline ParallelReduce(const FunctorType& arg_functor, Policy arg_policy,
                        const ReducerType& reducer)
      : m_result_ptr(reducer.view().data()),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible) {}

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

    ParReduceCommon::memcpy_result(ptr, &result, sizeof(ValueType),
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
#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

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

    ParReduceCommon::memcpy_result(ptr, &result, sizeof(ValueType),
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

    ParReduceCommon::memcpy_result(ptr, &result, sizeof(ValueType),
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

    ParReduceCommon::memcpy_result(ptr, &result, sizeof(ValueType),
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

    ParReduceCommon::memcpy_result(ptr, &result, sizeof(ValueType),
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
