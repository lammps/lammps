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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

namespace Impl {

  template <typename Scalar, class Device>
  struct test_dualview_combinations
  {
    typedef test_dualview_combinations<Scalar,Device> self_type;

    typedef Scalar scalar_type;
    typedef Device execution_space;

    Scalar reference;
    Scalar result;

    template <typename ViewType>
    Scalar run_me(unsigned int n,unsigned int m){
      if(n<10) n = 10;
      if(m<3) m = 3;
      ViewType a("A",n,m);

      Kokkos::deep_copy( a.d_view , 1 );

      a.template modify<typename ViewType::execution_space>();
      a.template sync<typename ViewType::host_mirror_space>();

      a.h_view(5,1) = 3;
      a.h_view(6,1) = 4;
      a.h_view(7,2) = 5;
      a.template modify<typename ViewType::host_mirror_space>();
      ViewType b = Kokkos::subview(a,std::pair<unsigned int, unsigned int>(6,9),std::pair<unsigned int, unsigned int>(0,1));
      a.template sync<typename ViewType::execution_space>();
      b.template modify<typename ViewType::execution_space>();

      Kokkos::deep_copy( b.d_view , 2 );

      a.template sync<typename ViewType::host_mirror_space>();
      Scalar count = 0;
      for(unsigned int i = 0; i<a.d_view.dimension_0(); i++)
        for(unsigned int j = 0; j<a.d_view.dimension_1(); j++)
          count += a.h_view(i,j);
      return count -  a.d_view.dimension_0()*a.d_view.dimension_1()-2-4-3*2;
    }


    test_dualview_combinations(unsigned int size)
    {
      result = run_me< Kokkos::DualView<Scalar**,Kokkos::LayoutLeft,Device> >(size,3);
    }

   };

} // namespace Impl




template <typename Scalar, typename Device>
void test_dualview_combinations(unsigned int size)
{
  Impl::test_dualview_combinations<Scalar,Device> test(size);
  ASSERT_EQ( test.result,0);

}


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP

