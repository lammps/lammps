/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifndef KOKKOS_IMPL_ANALYZE_POLICY_HPP
#define KOKKOS_IMPL_ANALYZE_POLICY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>
#include <impl/Kokkos_Tags.hpp>

namespace Kokkos { namespace Impl {

template < typename ExecutionSpace   = void
         , typename Schedule         = void
         , typename WorkTag          = void
         , typename IndexType        = void
         , typename IterationPattern = void
         , typename LaunchBounds     = void
         >
struct PolicyTraitsBase
{
  using type = PolicyTraitsBase< ExecutionSpace, Schedule, WorkTag, IndexType, 
               IterationPattern, LaunchBounds>;

  using execution_space   = ExecutionSpace;
  using schedule_type     = Schedule;
  using work_tag          = WorkTag;
  using index_type        = IndexType;
  using iteration_pattern = IterationPattern;
  using launch_bounds     = LaunchBounds;
};


template <typename PolicyBase, typename ExecutionSpace>
struct SetExecutionSpace
{
  static_assert( is_void<typename PolicyBase::execution_space>::value
               , "Kokkos Error: More than one execution space given" );
  using type = PolicyTraitsBase< ExecutionSpace
                               , typename PolicyBase::schedule_type
                               , typename PolicyBase::work_tag
                               , typename PolicyBase::index_type
                               , typename PolicyBase::iteration_pattern
                               , typename PolicyBase::launch_bounds
                               >;
};

template <typename PolicyBase, typename Schedule>
struct SetSchedule
{
  static_assert( is_void<typename PolicyBase::schedule_type>::value
               , "Kokkos Error: More than one schedule type given" );
  using type = PolicyTraitsBase< typename PolicyBase::execution_space
                               , Schedule
                               , typename PolicyBase::work_tag
                               , typename PolicyBase::index_type
                               , typename PolicyBase::iteration_pattern
                               , typename PolicyBase::launch_bounds
                               >;
};

template <typename PolicyBase, typename WorkTag>
struct SetWorkTag
{
  static_assert( is_void<typename PolicyBase::work_tag>::value
               , "Kokkos Error: More than one work tag given" );
  using type = PolicyTraitsBase< typename PolicyBase::execution_space
                               , typename PolicyBase::schedule_type
                               , WorkTag
                               , typename PolicyBase::index_type
                               , typename PolicyBase::iteration_pattern
                               , typename PolicyBase::launch_bounds
                               >;
};

template <typename PolicyBase, typename IndexType>
struct SetIndexType
{
  static_assert( is_void<typename PolicyBase::index_type>::value
               , "Kokkos Error: More than one index type given" );
  using type = PolicyTraitsBase< typename PolicyBase::execution_space
                               , typename PolicyBase::schedule_type
                               , typename PolicyBase::work_tag
                               , IndexType
                               , typename PolicyBase::iteration_pattern
                               , typename PolicyBase::launch_bounds
                               >;
};


template <typename PolicyBase, typename IterationPattern>
struct SetIterationPattern
{
  static_assert( is_void<typename PolicyBase::iteration_pattern>::value
               , "Kokkos Error: More than one iteration_pattern given" );
  using type = PolicyTraitsBase< typename PolicyBase::execution_space
                               , typename PolicyBase::schedule_type
                               , typename PolicyBase::work_tag
                               , typename PolicyBase::index_type
                               , IterationPattern
                               , typename PolicyBase::launch_bounds
                               >;
};


template <typename PolicyBase, typename LaunchBounds>
struct SetLaunchBounds
{
  static_assert( is_void<typename PolicyBase::launch_bounds>::value
               , "Kokkos Error: More than one launch_bounds given" );
  using type = PolicyTraitsBase< typename PolicyBase::execution_space
                               , typename PolicyBase::schedule_type
                               , typename PolicyBase::work_tag
                               , typename PolicyBase::index_type
                               , typename PolicyBase::iteration_pattern
                               , LaunchBounds
                               >;
};


template <typename Base, typename... Traits>
struct AnalyzePolicy;

template <typename Base, typename T, typename... Traits>
struct AnalyzePolicy<Base, T, Traits...> : public
  AnalyzePolicy<
      typename std::conditional< is_execution_space<T>::value  , SetExecutionSpace<Base,T>
    , typename std::conditional< is_schedule_type<T>::value    , SetSchedule<Base,T>
    , typename std::conditional< is_index_type<T>::value       , SetIndexType<Base,T>
    , typename std::conditional< std::is_integral<T>::value    , SetIndexType<Base, IndexType<T> >
    , typename std::conditional< is_iteration_pattern<T>::value, SetIterationPattern<Base,T>
    , typename std::conditional< is_launch_bounds<T>::value    , SetLaunchBounds<Base,T>
    , SetWorkTag<Base,T>
    >::type >::type >::type >::type >::type>::type::type
  , Traits...
  >
{};

template <typename Base>
struct AnalyzePolicy<Base>
{
  using execution_space = typename std::conditional< is_void< typename Base::execution_space >::value
                                                   , DefaultExecutionSpace
                                                   , typename Base::execution_space
                                                   >::type;

  using schedule_type = typename std::conditional< is_void< typename Base::schedule_type >::value
                                                 , Schedule< Static >
                                                 , typename Base::schedule_type
                                                 >::type;

  using work_tag = typename Base::work_tag;

  using index_type = typename std::conditional< is_void< typename Base::index_type >::value
                                              , IndexType< typename execution_space::size_type >
                                              , typename Base::index_type
                                              >::type
                                               ::type // nasty hack to make index_type into an integral_type
                                              ;       // instead of the wrapped IndexType<T> for backwards compatibility

  using iteration_pattern = typename std::conditional< is_void< typename Base::iteration_pattern >::value
                                                     , void // TODO set default iteration pattern
                                                     , typename Base::iteration_pattern
                                                     >::type;

  using launch_bounds = typename std::conditional< is_void< typename Base::launch_bounds >::value
                                                     , LaunchBounds<>
                                                     , typename Base::launch_bounds
                                                     >::type;

  using type = PolicyTraitsBase< execution_space
                               , schedule_type
                               , work_tag
                               , index_type
                               , iteration_pattern
                               , launch_bounds
                               >;
};

template <typename... Traits>
struct PolicyTraits
  : public AnalyzePolicy< PolicyTraitsBase<>, Traits... >::type
{};

}} // namespace Kokkos::Impl


#endif //KOKKOS_IMPL_ANALYZE_POLICY_HPP

