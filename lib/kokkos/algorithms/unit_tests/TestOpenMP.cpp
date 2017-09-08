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


#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_OPENMP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------
#include <TestRandom.hpp>
#include <TestSort.hpp>
#include <iomanip>

namespace Test {

class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;

    int threads_count = 0;
    #pragma omp parallel
    {
      #pragma omp atomic
      ++threads_count;
    }

    if (threads_count > 3) {
      threads_count /= 2;
    }

    Kokkos::OpenMP::initialize( threads_count );
    Kokkos::OpenMP::print_configuration( std::cout );
  }

  static void TearDownTestCase()
  {
    Kokkos::OpenMP::finalize();
  }
};

#define OPENMP_RANDOM_XORSHIFT64( num_draws )                                \
  TEST_F( openmp, Random_XorShift64 ) {   \
      Impl::test_random<Kokkos::Random_XorShift64_Pool<Kokkos::OpenMP> >(num_draws);                                   \
  }

#define OPENMP_RANDOM_XORSHIFT1024( num_draws )                                \
  TEST_F( openmp, Random_XorShift1024 ) {   \
      Impl::test_random<Kokkos::Random_XorShift1024_Pool<Kokkos::OpenMP> >(num_draws);                                   \
  }

#define OPENMP_SORT_UNSIGNED( size )                                \
  TEST_F( openmp, SortUnsigned ) {   \
      Impl::test_sort< Kokkos::OpenMP, unsigned >(size);                                   \
  }

OPENMP_RANDOM_XORSHIFT64( 10240000 )
OPENMP_RANDOM_XORSHIFT1024( 10130144 )
OPENMP_SORT_UNSIGNED(171)

#undef OPENMP_RANDOM_XORSHIFT64
#undef OPENMP_RANDOM_XORSHIFT1024
#undef OPENMP_SORT_UNSIGNED
} // namespace test
#else
void KOKKOS_ALGORITHMS_UNITTESTS_TESTOPENMP_PREVENT_LINK_ERROR() {}
#endif

