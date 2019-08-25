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

#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {
  template<class T,int N>
  class MyArray {
  public:
    T values[N];
    KOKKOS_INLINE_FUNCTION
    void operator+= (const MyArray& src) { for(int i=0; i<N; i++) values[i] += src.values[i]; }
    KOKKOS_INLINE_FUNCTION
    void operator=  (const MyArray& src) { for(int i=0; i<N; i++) values[i] =  src.values[i]; }
    KOKKOS_INLINE_FUNCTION
    void operator+= (const volatile MyArray& src) volatile { for(int i=0; i<N; i++) values[i] += src.values[i]; }
    KOKKOS_INLINE_FUNCTION
    void operator=  (const volatile MyArray& src) volatile { for(int i=0; i<N; i++) values[i] =  src.values[i]; }
  };

  template<class T,int N, class PolicyType, int S>
  struct FunctorFor {
    double static_array[S];
    KOKKOS_INLINE_FUNCTION
    void operator() (const typename PolicyType::member_type& team) const {
    }
  };
  template<class T,int N, class PolicyType, int S>
  struct FunctorReduce {
    double static_array[S];
    KOKKOS_INLINE_FUNCTION
    void operator() (const typename PolicyType::member_type& team, MyArray<T,N>& lval) const {
      for(int j=0; j<N; j++)
        lval.values[j] += 1 + lval.values[0];
    }
  };
}


typedef Kokkos::TeamPolicy<TEST_EXECSPACE> policy_type;
typedef Kokkos::TeamPolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128,8> > policy_type_128_8;
typedef Kokkos::TeamPolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<1024,2> > policy_type_1024_2;

template<class T, int N, class PolicyType, int S>
void test_team_policy_max_recommended_static_size(int scratch_size) {
  PolicyType p = PolicyType(10000, Kokkos::AUTO, 4).set_scratch_size(0,Kokkos::PerTeam(scratch_size));
  int team_size_max_for = p.team_size_max(FunctorFor<T,N,PolicyType,S>(),Kokkos::ParallelForTag());
  int team_size_rec_for = p.team_size_recommended(FunctorFor<T,N,PolicyType,S>(),Kokkos::ParallelForTag());
  int team_size_max_reduce = p.team_size_max(FunctorReduce<T,N,PolicyType,S>(),Kokkos::ParallelReduceTag());
  int team_size_rec_reduce = p.team_size_recommended(FunctorReduce<T,N,PolicyType,S>(),Kokkos::ParallelReduceTag());

  ASSERT_TRUE( team_size_max_for >= team_size_rec_for );
  ASSERT_TRUE( team_size_max_reduce >= team_size_rec_reduce );
  ASSERT_TRUE( team_size_max_for >= team_size_max_reduce );

  Kokkos::parallel_for(PolicyType(10000, team_size_max_for, 4).set_scratch_size(0,Kokkos::PerTeam(scratch_size)),
                       FunctorFor<T,N,PolicyType,S>());
  Kokkos::parallel_for(PolicyType(10000, team_size_rec_for, 4).set_scratch_size(0,Kokkos::PerTeam(scratch_size)),
                       FunctorFor<T,N,PolicyType,S>());
  MyArray<T,N> val;
  Kokkos::parallel_reduce(PolicyType(10000, team_size_max_reduce, 4).set_scratch_size(0,Kokkos::PerTeam(scratch_size)),
                       FunctorReduce<T,N,PolicyType,S>(),val);
  Kokkos::parallel_reduce(PolicyType(10000, team_size_rec_reduce, 4).set_scratch_size(0,Kokkos::PerTeam(scratch_size)),
                       FunctorReduce<T,N,PolicyType,S>(),val);
  Kokkos::fence();
}

template<class T, int N, class PolicyType>
void test_team_policy_max_recommended(int scratch_size) {
  test_team_policy_max_recommended_static_size<T,N,PolicyType,1>(scratch_size);
  test_team_policy_max_recommended_static_size<T,N,PolicyType,1000>(scratch_size);
}

TEST_F( TEST_CATEGORY, team_policy_max_recommended )
{
  int max_scratch_size = policy_type::scratch_size_max(0);
  test_team_policy_max_recommended<double,2,policy_type>(0);
  test_team_policy_max_recommended<double,2,policy_type>(max_scratch_size/3);
  test_team_policy_max_recommended<double,2,policy_type>(max_scratch_size);
  test_team_policy_max_recommended<double,2,policy_type_128_8>(0);
  test_team_policy_max_recommended<double,2,policy_type_128_8>(max_scratch_size/3/8);
  test_team_policy_max_recommended<double,2,policy_type_128_8>(max_scratch_size/8);
  test_team_policy_max_recommended<double,2,policy_type_1024_2>(0);
  test_team_policy_max_recommended<double,2,policy_type_1024_2>(max_scratch_size/3/2);
  test_team_policy_max_recommended<double,2,policy_type_1024_2>(max_scratch_size/2);

  test_team_policy_max_recommended<double,16,policy_type>(0);
  test_team_policy_max_recommended<double,16,policy_type>(max_scratch_size/3);
  test_team_policy_max_recommended<double,16,policy_type>(max_scratch_size);
  test_team_policy_max_recommended<double,16,policy_type_128_8>(0);
  test_team_policy_max_recommended<double,16,policy_type_128_8>(max_scratch_size/3/8);
  test_team_policy_max_recommended<double,16,policy_type_128_8>(max_scratch_size/8);
  test_team_policy_max_recommended<double,16,policy_type_1024_2>(0);
  test_team_policy_max_recommended<double,16,policy_type_1024_2>(max_scratch_size/3/2);
  test_team_policy_max_recommended<double,16,policy_type_1024_2>(max_scratch_size/2);
}


} // namespace Test
