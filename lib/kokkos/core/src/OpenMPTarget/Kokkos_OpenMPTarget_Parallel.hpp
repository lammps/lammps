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
  void execute() const { execute_impl<WorkTag>(); }

  template <class TagType>
  void execute_impl() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = m_policy.begin();
    const auto end   = m_policy.end();

    if (end <= begin) return;

    FunctorType a_functor(m_functor);

#pragma omp target teams distribute parallel for map(to : a_functor)
    for (auto i = begin; i < end; ++i) {
      if constexpr (std::is_same<TagType, void>::value) {
        a_functor(i);
      } else {
        a_functor(TagType(), i);
      }
    }
  }

  ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// This class has the memcpy routine that is commonly used by ParallelReduce
// over RangePolicy and TeamPolicy.
template <class PointerType>
struct ParallelReduceCommon {
  // Copy the result back to device if the view is on the device.
  static void memcpy_result(PointerType dest, PointerType src, size_t size,
                            bool ptr_on_device) {
    if (ptr_on_device) {
      OMPT_SAFE_CALL(omp_target_memcpy(dest, src, size, 0, 0,
                                       omp_get_default_device(),
                                       omp_get_initial_device()));
    } else {
      *dest = *src;
    }
  }
};

template <class FunctorType, class PolicyType, class ReducerType,
          class PointerType, class ValueType>
struct ParallelReduceSpecialize {
  inline static void execute(const FunctorType& /*f*/, const PolicyType& /*p*/,
                             PointerType /*result_ptr*/) {
    constexpr int FunctorHasJoin = ReduceFunctorHasJoin<FunctorType>::value;
    constexpr int UseReducerType = is_reducer_type<ReducerType>::value;

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
      typename std::conditional<std::is_same<InvalidType, ReducerType>::value,
                                FunctorType, ReducerType>::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, TagType,
                         void>;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit     = Kokkos::Impl::FunctorValueInit<FunctorType, TagType>;
  using ValueJoin     = Kokkos::Impl::FunctorValueJoin<FunctorType, TagType>;
  using ReferenceType = typename ValueTraits::reference_type;

  using ParReduceCommon = ParallelReduceCommon<PointerType>;

  static void execute_reducer(const FunctorType& f, const PolicyType& p,
                              PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = p.begin();
    const auto end   = p.end();

    ValueType result;
    OpenMPTargetReducerWrapper<ReducerType>::init(result);

    // Initialize and copy back the result even if it is a zero length
    // reduction.
    if (end <= begin) {
      ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
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
      if constexpr (std::is_same<TagType, void>::value) {
        f(i, result);
      } else {
        f(TagType(), i, result);
      }
    }

    ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                   ptr_on_device);
  }

  template <class TagType, int NumReductions>
  static void execute_array(const FunctorType& f, const PolicyType& p,
                            PointerType result_ptr, bool ptr_on_device) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const auto begin = p.begin();
    const auto end   = p.end();

    // Enter the loop if the reduction is on a scalar type.
    if constexpr (NumReductions == 1) {
      ValueType result = ValueType();

      // Initialize and copy back the result even if it is a zero length
      // reduction.
      if (end <= begin) {
        ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                       ptr_on_device);
        return;
      }
      // Case where reduction is on a native data type.
      if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp target teams distribute parallel for \
         map(to:f) reduction(+: result)
        for (auto i = begin; i < end; ++i)

          if constexpr (std::is_same<TagType, void>::value) {
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

          if constexpr (std::is_same<TagType, void>::value) {
            f(i, result);
          } else {
            f(TagType(), i, result);
          }
      }

      ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                     ptr_on_device);
    } else {
      ValueType result[NumReductions] = {};

      // Initialize and copy back the result even if it is a zero length
      // reduction.
      if (end <= begin) {
        ParReduceCommon::memcpy_result(result_ptr, result,
                                       NumReductions * sizeof(ValueType),
                                       ptr_on_device);
        return;
      }
#pragma omp target teams distribute parallel for map(to:f) reduction(+:result[:NumReductions])
      for (auto i = begin; i < end; ++i) {
        if constexpr (std::is_same<TagType, void>::value) {
          f(i, result);
        } else {
          f(TagType(), i, result);
        }
      }

      ParReduceCommon::memcpy_result(
          result_ptr, result, NumReductions * sizeof(ValueType), ptr_on_device);
    }
  }

  static void execute_init_join(const FunctorType& f, const PolicyType& p,
                                PointerType ptr, const bool ptr_on_device) {
    const auto begin = p.begin();
    const auto end   = p.end();

    constexpr int HasInit = ReduceFunctorHasInit<FunctorType>::value;

    // Initialize the result pointer.

    const auto size = end - begin;

    // FIXME_OPENMPTARGET: The team size and MAX_ACTIVE_THREADS are currently
    // based on NVIDIA-V100 and should be modifid to be based on the
    // architecture in the future.
    const int max_team_threads = 32;
    const int max_teams =
        OpenMPTargetExec::MAX_ACTIVE_THREADS / max_team_threads;
    // Number of elements in the reduction
    const auto value_count =
        FunctorValueTraits<FunctorType, TagType>::value_count(f);

    // Allocate scratch per active thread. Achieved by setting the first
    // parameter of `resize_scratch=1`.
    OpenMPTargetExec::resize_scratch(1, 0, value_count * sizeof(ValueType));
    ValueType* scratch_ptr =
        static_cast<ValueType*>(OpenMPTargetExec::get_scratch_ptr());

#pragma omp target map(to : f) is_device_ptr(scratch_ptr)
    {
      // Enter this loop if the functor has an `init`
      if constexpr (HasInit) {
        // The `init` routine needs to be called on the device since it might
        // need device members.
        ValueInit::init(f, scratch_ptr);
        if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
          FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
      } else {
        for (int i = 0; i < value_count; ++i) {
          static_cast<ValueType*>(scratch_ptr)[i] = ValueType();
        }

        if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
          FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
      }
    }

    if (end <= begin) {
      // If there is no work to be done, copy back the initialized values and
      // exit.
      if (!ptr_on_device)
        OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_initial_device(), omp_get_default_device()));
      else
        OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_default_device(), omp_get_default_device()));

      return;
    }

#pragma omp target teams num_teams(max_teams) thread_limit(max_team_threads) \
    map(to                                                                   \
        : f) is_device_ptr(scratch_ptr)
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
        ReferenceType result = ValueInit::init(
            f, &team_scratch[omp_get_thread_num() * value_count]);

        // Accumulate partial results in thread specific storage.
#pragma omp for simd
        for (auto i = team_begin; i < team_end; ++i) {
          if constexpr (std::is_same<TagType, void>::value) {
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
            ValueJoin::join(f, &team_scratch[i * value_count],
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
        ValueJoin::join(
            f, &team_scratch[i * team_offset],
            &team_scratch[(i + tree_neighbor_offset) * team_offset]);

        // If `final` is provided by the functor.
        if constexpr (ReduceFunctorHasFinal<FunctorType>::value) {
          // Do the final only once at the end.
          if (tree_neighbor_offset * 2 >= max_teams &&
              omp_get_team_num() == 0 && omp_get_thread_num() == 0)
            FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
        }
      }
      tree_neighbor_offset *= 2;
    } while (tree_neighbor_offset < max_teams);

    // If the result view is on the host, copy back the values via memcpy.
    if (!ptr_on_device)
      OMPT_SAFE_CALL(omp_target_memcpy(
          ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
          omp_get_initial_device(), omp_get_default_device()));
    else
      OMPT_SAFE_CALL(omp_target_memcpy(
          ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
          omp_get_default_device(), omp_get_default_device()));
  }
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;

  using ReducerTypeFwd =
      typename std::conditional<std::is_same<InvalidType, ReducerType>::value,
                                FunctorType, ReducerType>::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, WorkTag,
                         void>;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;

  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;

  static constexpr int HasJoin    = ReduceFunctorHasJoin<FunctorType>::value;
  static constexpr int UseReducer = is_reducer_type<ReducerType>::value;
  static constexpr int IsArray    = std::is_pointer<reference_type>::value;

  using ParReduceSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy, ReducerType, pointer_type,
                               typename ValueTraits::value_type>;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  bool m_result_ptr_on_device;
  const int m_result_ptr_num_elems;
  using TagType = typename Policy::work_tag;

 public:
  void execute() const {
    if constexpr (HasJoin) {
      // Enter this loop if the Functor has a init-join.
      ParReduceSpecialize::execute_init_join(m_functor, m_policy, m_result_ptr,
                                             m_result_ptr_on_device);
    } else if constexpr (UseReducer) {
      // Enter this loop if the Functor is a reducer type.
      ParReduceSpecialize::execute_reducer(m_functor, m_policy, m_result_ptr,
                                           m_result_ptr_on_device);
    } else if constexpr (IsArray) {
      // Enter this loop if the reduction is on an array and the routine is
      // templated over the size of the array.
      if (m_result_ptr_num_elems <= 2) {
        ParReduceSpecialize::template execute_array<TagType, 2>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 4) {
        ParReduceSpecialize::template execute_array<TagType, 4>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 8) {
        ParReduceSpecialize::template execute_array<TagType, 8>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 16) {
        ParReduceSpecialize::template execute_array<TagType, 16>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 32) {
        ParReduceSpecialize::template execute_array<TagType, 32>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else {
        Kokkos::abort("array reduction length must be <= 32");
      }
    } else {
      // This loop handles the basic scalar reduction.
      ParReduceSpecialize::template execute_array<TagType, 1>(
          m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
    }
  }

  template <class ViewType>
  ParallelReduce(
      const FunctorType& arg_functor, Policy& arg_policy,
      const ViewType& arg_result_view,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = nullptr)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_num_elems(arg_result_view.size()) {}

  ParallelReduce(const FunctorType& arg_functor, Policy& arg_policy,
                 const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible),
        m_result_ptr_num_elems(reducer.view().size()) {}
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
  typename std::enable_if<std::is_same<TagType, void>::value>::type
  call_with_tag(const FunctorType& f, const idx_type& idx, value_type& val,
                const bool& is_final) const {
    f(idx, val, is_final);
  }
  template <class TagType>
  typename std::enable_if<!std::is_same<TagType, void>::value>::type
  call_with_tag(const FunctorType& f, const idx_type& idx, value_type& val,
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

    FunctorType a_functor(m_functor);
#pragma omp target teams distribute map(to                             \
                                        : a_functor) num_teams(nteams) \
    thread_limit(team_size)
    for (idx_type team_id = 0; team_id < n_chunks; ++team_id) {
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;

#pragma omp for
        for (idx_type i = 0; i < chunk_size; ++i) {
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
          for (idx_type i = 0; i < chunk_size; ++i) {
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
            for (idx_type i = 0; i < n_chunks; ++i) {
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
    for (idx_type team_id = 0; team_id < n_chunks; ++team_id) {
#pragma omp parallel num_threads(team_size)
      {
        const idx_type local_offset = team_id * chunk_size;
        value_type offset_value;
        if (team_id > 0)
          offset_value = chunk_values(team_id - 1);
        else
          ValueInit::init(a_functor, &offset_value);

#pragma omp for
        for (idx_type i = 0; i < chunk_size; ++i) {
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

  ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
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
  void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    execute_impl<WorkTag>();
  }

 private:
  template <class TagType>
  void execute_impl() const {
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
  }

 public:
  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {}
};

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType, TeamPolicyInternal<PolicyArgs...>,
                                ReducerType, PointerType, ValueType> {
  using PolicyType = TeamPolicyInternal<PolicyArgs...>;
  using TagType    = typename PolicyType::work_tag;
  using ReducerTypeFwd =
      typename std::conditional<std::is_same<InvalidType, ReducerType>::value,
                                FunctorType, ReducerType>::type;
  using WorkTagFwd =
      std::conditional_t<std::is_same<InvalidType, ReducerType>::value, TagType,
                         void>;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit     = Kokkos::Impl::FunctorValueInit<FunctorType, TagType>;
  using ValueJoin     = Kokkos::Impl::FunctorValueJoin<FunctorType, TagType>;
  using ReferenceType = typename ValueTraits::reference_type;

  using ParReduceCommon = ParallelReduceCommon<PointerType>;

  static void execute_reducer(const FunctorType& f, const PolicyType& p,
                              PointerType result_ptr, bool ptr_on_device) {
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

#pragma omp declare reduction(                                         \
    custom:ValueType                                                   \
    : OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

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

    // Copy results back to device if `parallel_reduce` is on a device view.
    ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                   ptr_on_device);
  }

  template <int NumReductions>
  static void execute_array(const FunctorType& f, const PolicyType& p,
                            PointerType result_ptr, bool ptr_on_device) {
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

    // Maximum active teams possible.
    int max_active_teams = OpenMPTargetExec::MAX_ACTIVE_THREADS / team_size;
    const auto nteams =
        league_size < max_active_teams ? league_size : max_active_teams;

    // Case where the number of reduction items is 1.
    if constexpr (NumReductions == 1) {
      ValueType result = ValueType();

      // Case where reduction is on a native data type.
      if constexpr (std::is_arithmetic<ValueType>::value) {
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
      } else {
        // Case where the reduction is on a non-native data type.
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)
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
      }

      // Copy results back to device if `parallel_reduce` is on a device view.
      ParReduceCommon::memcpy_result(result_ptr, &result, sizeof(ValueType),
                                     ptr_on_device);
    } else {
      ValueType result[NumReductions] = {};
      // Case where the reduction is on an array.
#pragma omp target teams num_teams(nteams) thread_limit(team_size) map(to   \
                                                                       : f) \
    is_device_ptr(scratch_ptr) reduction(+ : result[:NumReductions])
#pragma omp parallel reduction(+ : result[:NumReductions])
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

      // Copy results back to device if `parallel_reduce` is on a device view.
      ParReduceCommon::memcpy_result(
          result_ptr, result, NumReductions * sizeof(ValueType), ptr_on_device);
    }
  }

  // FIXME_OPENMPTARGET : This routine is a copy from `parallel_reduce` over
  // RangePolicy. Need a new implementation.
  static void execute_init_join(const FunctorType& f, const PolicyType& p,
                                PointerType ptr, const bool ptr_on_device) {
    constexpr int HasInit = ReduceFunctorHasInit<FunctorType>::value;

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
    const auto value_count =
        FunctorValueTraits<FunctorType, TagType>::value_count(f);

    // Allocate scratch per active thread.
    OpenMPTargetExec::resize_scratch(1, 0, value_count * sizeof(ValueType));
    void* scratch_ptr = OpenMPTargetExec::get_scratch_ptr();

    // Enter this loop if the functor has an `init`
    if constexpr (HasInit) {
      // The `init` routine needs to be called on the device since it might need
      // device members.
#pragma omp target map(to : f) is_device_ptr(scratch_ptr)
      {
        ValueInit::init(f, scratch_ptr);

        if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
          FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
      }
    } else {
#pragma omp target map(to : f) is_device_ptr(scratch_ptr)
      {
        for (int i = 0; i < value_count; ++i) {
          static_cast<ValueType*>(scratch_ptr)[i] = ValueType();
        }

        if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
          FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
      }
    }

    if (end <= begin) {
      // If there is no work to be done, copy back the initialized values and
      // exit.
      if (!ptr_on_device)
        OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_initial_device(), omp_get_default_device()));
      else
        OMPT_SAFE_CALL(omp_target_memcpy(
            ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
            omp_get_default_device(), omp_get_default_device()));

      return;
    }

#pragma omp target teams num_teams(nteams) thread_limit(team_size) map(to   \
                                                                       : f) \
    is_device_ptr(scratch_ptr)
    {
#pragma omp parallel
      {
        const int team_num      = omp_get_team_num();
        const int num_teams     = omp_get_num_teams();
        ValueType* team_scratch = static_cast<ValueType*>(scratch_ptr) +
                                  team_num * team_size * value_count;
        ReferenceType result = ValueInit::init(f, &team_scratch[0]);

        for (int league_id = team_num; league_id < league_size;
             league_id += num_teams) {
          typename PolicyType::member_type team(
              league_id, league_size, team_size, vector_length, scratch_ptr,
              team_num, shmem_size_L0, shmem_size_L1);
          if constexpr (std::is_same<TagType, void>::value) {
            f(team, result);
          } else {
            f(TagType(), team, result);
          }
        }
      }  // end parallel
    }    // end target

    int tree_neighbor_offset = 1;
    do {
#pragma omp target teams distribute parallel for simd map(to   \
                                                          : f) \
    is_device_ptr(scratch_ptr)
      for (int i = 0; i < nteams - tree_neighbor_offset;
           i += 2 * tree_neighbor_offset) {
        ValueType* team_scratch = static_cast<ValueType*>(scratch_ptr);
        const int team_offset   = team_size * value_count;
        ValueJoin::join(
            f, &team_scratch[i * team_offset],
            &team_scratch[(i + tree_neighbor_offset) * team_offset]);

        // If `final` is provided by the functor.
        if constexpr (ReduceFunctorHasFinal<FunctorType>::value) {
          // Do the final only once at the end.
          if (tree_neighbor_offset * 2 >= nteams && omp_get_team_num() == 0 &&
              omp_get_thread_num() == 0)
            FunctorFinal<FunctorType, TagType>::final(f, scratch_ptr);
        }
      }
      tree_neighbor_offset *= 2;
    } while (tree_neighbor_offset < nteams);

    // If the result view is on the host, copy back the values via memcpy.
    if (!ptr_on_device)
      OMPT_SAFE_CALL(omp_target_memcpy(
          ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
          omp_get_initial_device(), omp_get_default_device()));
    else
      OMPT_SAFE_CALL(omp_target_memcpy(
          ptr, scratch_ptr, value_count * sizeof(ValueType), 0, 0,
          omp_get_default_device(), omp_get_default_device()));
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
  using ReducerTypeFwd =
      typename std::conditional<std::is_same<InvalidType, ReducerType>::value,
                                FunctorType, ReducerType>::type;
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

  bool m_result_ptr_on_device;
  const int m_result_ptr_num_elems;

  static constexpr int HasJoin    = ReduceFunctorHasJoin<FunctorType>::value;
  static constexpr int UseReducer = is_reducer_type<ReducerType>::value;
  static constexpr int IsArray    = std::is_pointer<reference_type>::value;

  using ParReduceSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy, ReducerType, pointer_type,
                               typename ValueTraits::value_type>;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const int m_shmem_size;

 public:
  void execute() const {
    if constexpr (HasJoin) {
      ParReduceSpecialize::execute_init_join(m_functor, m_policy, m_result_ptr,
                                             m_result_ptr_on_device);
    } else if constexpr (UseReducer) {
      ParReduceSpecialize::execute_reducer(m_functor, m_policy, m_result_ptr,
                                           m_result_ptr_on_device);
    } else if constexpr (IsArray) {
      if (m_result_ptr_num_elems <= 2) {
        ParReduceSpecialize::template execute_array<2>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 4) {
        ParReduceSpecialize::template execute_array<4>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 8) {
        ParReduceSpecialize::template execute_array<8>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 16) {
        ParReduceSpecialize::template execute_array<16>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 32) {
        ParReduceSpecialize::template execute_array<32>(
            m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else {
        Kokkos::abort("array reduction length must be <= 32");
      }
    } else {
      ParReduceSpecialize::template execute_array<1>(
          m_functor, m_policy, m_result_ptr, m_result_ptr_on_device);
    }
  }

  template <class ViewType>
  ParallelReduce(
      const FunctorType& arg_functor, const Policy& arg_policy,
      const ViewType& arg_result,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = nullptr)
      : m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_num_elems(arg_result.size()),
        m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_shmem_size(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                     FunctorTeamShmemSize<FunctorType>::value(
                         arg_functor, arg_policy.team_size())) {}

  ParallelReduce(const FunctorType& arg_functor, Policy& arg_policy,
                 const ReducerType& reducer)
      : m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible),
        m_result_ptr_num_elems(reducer.view().size()),
        m_functor(arg_functor),
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

  TeamThreadRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  iType count)
      : start(0), end(count), team(thread_) {}
  TeamThreadRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  iType begin_, iType end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  ThreadVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                    index_type count)
      : start(0), end(count), team(thread_) {}
  ThreadVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                    index_type begin_, index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  TeamVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  index_type count)
      : start(0), end(count), team(thread_) {}
  TeamVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  index_type begin_, index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

}  // namespace Impl

}  // namespace Kokkos
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */
