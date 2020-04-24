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
#include <iostream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Exec.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Experimental::OpenMPTarget> {
 private:
  typedef Kokkos::MDRangePolicy<Traits...> Policy;
  typedef typename Policy::work_tag WorkTag;
  typedef typename Policy::member_type Member;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  inline void execute() const {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const int64_t begin = 0;
    const int64_t end   = m_policy.m_num_tiles;
    FunctorType functor(m_functor);
    Policy policy = m_policy;
#pragma omp target teams distribute map(to : functor) num_teams(end - begin)
    {
      for (ptrdiff_t tile_idx = begin; tile_idx < end; tile_idx++) {
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
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 1>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];
#pragma omp for
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++) {
      functor(i0);
    }
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 2>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
    const ptrdiff_t begin_0 = offset[0];
    ptrdiff_t end_0         = begin_0 + policy.m_tile[0];
    end_0 = end_0 < policy.m_upper[0] ? end_0 : policy.m_upper[0];

    const ptrdiff_t begin_1 = offset[1];
    ptrdiff_t end_1         = begin_1 + policy.m_tile[1];
    end_1 = end_1 < policy.m_upper[1] ? end_1 : policy.m_upper[1];

#pragma omp for collapse(2)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++) functor(i0, i1);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 3>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++) functor(i0, i1, i2);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 4>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; i3++)
            functor(i0, i1, i2, i3);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 5>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; i3++)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; i4++)
              functor(i0, i1, i2, i3, i4);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 6>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; i3++)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; i4++)
              for (ptrdiff_t i5 = begin_5; i5 < end_5; i5++)
                functor(i0, i1, i2, i3, i4, i5);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 7>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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

    const ptrdiff_t begin_6 = offset[6];
    ptrdiff_t end_6         = begin_6 + policy.m_tile[6];
    end_6 = end_6 < policy.m_upper[6] ? end_6 : policy.m_upper[6];

#pragma omp for collapse(7)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; i3++)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; i4++)
              for (ptrdiff_t i5 = begin_5; i5 < end_5; i5++)
                for (ptrdiff_t i6 = begin_6; i6 < end_6; i6++)
                  functor(i0, i1, i2, i3, i4, i5, i6);
  }

  template <int Rank>
  inline typename std::enable_if<Rank == 8>::type execute_tile(
      typename Policy::point_type offset, const FunctorType& functor,
      const Policy& policy) const {
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

    const ptrdiff_t begin_6 = offset[6];
    ptrdiff_t end_6         = begin_6 + policy.m_tile[6];
    end_6 = end_6 < policy.m_upper[6] ? end_6 : policy.m_upper[6];

    const ptrdiff_t begin_7 = offset[7];
    ptrdiff_t end_7         = begin_7 + policy.m_tile[7];
    end_7 = end_7 < policy.m_upper[7] ? end_7 : policy.m_upper[7];

#pragma omp for collapse(8)
    for (ptrdiff_t i0 = begin_0; i0 < end_0; i0++)
      for (ptrdiff_t i1 = begin_1; i1 < end_1; i1++)
        for (ptrdiff_t i2 = begin_2; i2 < end_2; i2++)
          for (ptrdiff_t i3 = begin_3; i3 < end_3; i3++)
            for (ptrdiff_t i4 = begin_4; i4 < end_4; i4++)
              for (ptrdiff_t i5 = begin_5; i5 < end_5; i5++)
                for (ptrdiff_t i6 = begin_6; i6 < end_6; i6++)
                  for (ptrdiff_t i7 = begin_7; i7 < end_7; i7++)
                    functor(i0, i1, i2, i3, i4, i5, i6, i7);
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

template <class FunctorType, class ReducerType, class PointerType,
          class ValueType, class... PolicyArgs>
struct ParallelReduceSpecialize<FunctorType,
                                Kokkos::MDRangePolicy<PolicyArgs...>,
                                ReducerType, PointerType, ValueType, 0, 0> {
  typedef Kokkos::RangePolicy<PolicyArgs...> PolicyType;
  template <class TagType>
  inline static
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      execute_impl(const FunctorType& f, const PolicyType& p,
                   PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const typename PolicyType::member_type begin = p.begin();
    const typename PolicyType::member_type end   = p.end();

    ValueType result = ValueType();
#pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom:result) reduction(+: result)
    for (int i = begin; i < end; i++) f(i, result);

    *result_ptr = result;
  }

  template <class TagType>
  inline static
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      execute_impl(const FunctorType& f, const PolicyType& p,
                   PointerType result_ptr) {
    OpenMPTargetExec::verify_is_process(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    OpenMPTargetExec::verify_initialized(
        "Kokkos::Experimental::OpenMPTarget parallel_for");
    const typename PolicyType::member_type begin = p.begin();
    const typename PolicyType::member_type end   = p.end();

    ValueType result = ValueType();
#pragma omp target teams distribute parallel for num_teams(512) map(to:f) map(tofrom: result) reduction(+: result)
    for (int i = begin; i < end; i++) f(TagType(), i, result);

    *result_ptr = result;
  }

  inline static void execute(const FunctorType& f, const PolicyType& p,
                             PointerType ptr) {
    execute_impl<typename PolicyType::work_tag>(f, p, ptr);
  }
};
/*
template<class FunctorType, class PolicyType, class ReducerType, class
PointerType, class ValueType> struct ParallelReduceSpecialize<FunctorType,
PolicyType, ReducerType, PointerType, ValueType, 0,1> {

  #pragma omp declare reduction(custom: ValueType : ReducerType::join(omp_out,
omp_in)) initializer ( ReducerType::init(omp_priv) )

  template< class TagType >
  inline static
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType
result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget
parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget
parallel_for"); const typename PolicyType::member_type begin = p.begin(); const
typename PolicyType::member_type end = p.end();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f)
map(tofrom:result) reduction(custom: result) for(int i=begin; i<end; i++)
        f(i,result);

      *result_ptr=result;
    }


  template< class TagType >
  inline static
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  execute_impl(const FunctorType& f, const PolicyType& p, PointerType
result_ptr)
    {
      OpenMPTargetExec::verify_is_process("Kokkos::Experimental::OpenMPTarget
parallel_for");
      OpenMPTargetExec::verify_initialized("Kokkos::Experimental::OpenMPTarget
parallel_for"); const typename PolicyType::member_type begin = p.begin(); const
typename PolicyType::member_type end = p.end();

      ValueType result = ValueType();
      #pragma omp target teams distribute parallel for num_teams(512) map(to:f)
map(tofrom: result) reduction(custom: result) for(int i=begin; i<end; i++)
        f(TagType(),i,result);

      *result_ptr=result;
    }


    inline static
    void execute(const FunctorType& f, const PolicyType& p, PointerType ptr) {
      execute_impl<typename PolicyType::work_tag>(f,p,ptr);
    }
};


template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  typedef Kokkos::MDRangePolicy<Traits...> Policy;

  typedef typename Policy::work_tag WorkTag;
  typedef typename Policy::WorkRange WorkRange;
  typedef typename Policy::member_type Member;

  typedef Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                             FunctorType, ReducerType>
      ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type WorkTagFwd;

  // Static Assert WorkTag void if ReducerType not InvalidType

  typedef Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>
      ValueTraits;
  typedef Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd> ValueInit;
  typedef Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd> ValueJoin;

  enum { HasJoin = ReduceFunctorHasJoin<FunctorType>::value };
  enum { UseReducer = is_reducer_type<ReducerType>::value };

  typedef typename ValueTraits::pointer_type pointer_type;
  typedef typename ValueTraits::reference_type reference_type;

  typedef ParallelReduceSpecialize<
      FunctorType, Policy, ReducerType, pointer_type,
      typename ValueTraits::value_type, HasJoin, UseReducer>
      ParForSpecialize;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;

 public:
  inline void execute() const {
    ParForSpecialize::execute(m_functor, m_policy, m_result_ptr);
  }

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& arg_functor, Policy arg_policy,
      const ViewType& arg_result_view,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = NULL)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result_view.data()) {
    //static_assert( std::is_same< typename ViewType::memory_space
    //                                , Kokkos::HostSpace >::value
    //  , "Reduction result on Kokkos::Experimental::OpenMPTarget must be a
    //  Kokkos::View in HostSpace" );
  }

  inline ParallelReduce(const FunctorType& arg_functor, Policy arg_policy,
                        const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {
    //static_assert( std::is_same< typename ViewType::memory_space
    //                                , Kokkos::HostSpace >::value
    //  , "Reduction result on Kokkos::Experimental::OpenMPTarget must be a
    //  Kokkos::View in HostSpace" );
  }
};*/

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */
