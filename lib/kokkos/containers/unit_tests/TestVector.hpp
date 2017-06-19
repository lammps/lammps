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

#ifndef KOKKOS_TEST_VECTOR_HPP
#define KOKKOS_TEST_VECTOR_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

namespace Impl {

  template <typename Scalar, class Device>
  struct test_vector_combinations
  {
    typedef test_vector_combinations<Scalar,Device> self_type;

    typedef Scalar scalar_type;
    typedef Device execution_space;

    Scalar reference;
    Scalar result;

    template <typename Vector>
    Scalar run_me(unsigned int n){
      Vector a(n,1);


      a.push_back(2);
      a.resize(n+4);
      a[n+1] = 3;
      a[n+2] = 4;
      a[n+3] = 5;


      Scalar temp1 = a[2];
      Scalar temp2 = a[n];
      Scalar temp3 = a[n+1];

      a.assign(n+2,-1);

      a[2] = temp1;
      a[n] = temp2;
      a[n+1] = temp3;

      Scalar test1 = 0;
      for(unsigned int i=0; i<a.size(); i++)
        test1+=a[i];

      a.assign(n+1,-2);
      Scalar test2 = 0;
      for(unsigned int i=0; i<a.size(); i++)
        test2+=a[i];

      a.reserve(n+10);

      Scalar test3 = 0;
      for(unsigned int i=0; i<a.size(); i++)
        test3+=a[i];


      return (test1*test2+test3)*test2+test1*test3;
    }


    test_vector_combinations(unsigned int size)
    {
      reference = run_me<std::vector<Scalar> >(size);
      result = run_me<Kokkos::vector<Scalar,Device> >(size);
    }

   };

} // namespace Impl




template <typename Scalar, typename Device>
void test_vector_combinations(unsigned int size)
{
  Impl::test_vector_combinations<Scalar,Device> test(size);
  ASSERT_EQ( test.reference, test.result);
}


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP

