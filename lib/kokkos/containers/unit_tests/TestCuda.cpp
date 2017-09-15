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
#ifdef KOKKOS_ENABLE_CUDA

#include <iostream>
#include <iomanip>
#include <cstdint>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <Kokkos_Bitset.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_Vector.hpp>

#include <TestBitset.hpp>
#include <TestUnorderedMap.hpp>
#include <TestStaticCrsGraph.hpp>
#include <TestVector.hpp>
#include <TestDualView.hpp>
#include <TestDynamicView.hpp>

#include <Kokkos_DynRankView.hpp>
#include <TestDynViewAPI.hpp>

#include <Kokkos_ErrorReporter.hpp>
#include <TestErrorReporter.hpp>

#include <TestViewCtorPropEmbeddedDim.hpp>

//----------------------------------------------------------------------------



namespace Test {

class cuda : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
  }
  static void TearDownTestCase()
  {
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }
};

TEST_F( cuda , dyn_view_api) {
  TestDynViewAPI< double , Kokkos::Cuda >();
}

TEST_F( cuda, viewctorprop_embedded_dim ) {
  TestViewCtorProp_EmbeddedDim< Kokkos::Cuda >::test_vcpt( 2, 3 );
}

TEST_F( cuda , staticcrsgraph )
{
  TestStaticCrsGraph::run_test_graph< Kokkos::Cuda >();
  TestStaticCrsGraph::run_test_graph2< Kokkos::Cuda >();
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(1, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(1, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(1, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(1, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(3, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(3, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(3, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(3, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(75, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(75, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(75, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Cuda >(75, 100000);
}


void cuda_test_insert_close(  uint32_t num_nodes
                            , uint32_t num_inserts
                            , uint32_t num_duplicates
                           )
{
  test_insert< Kokkos::Cuda >( num_nodes, num_inserts, num_duplicates, true);
}

void cuda_test_insert_far(  uint32_t num_nodes
                          , uint32_t num_inserts
                          , uint32_t num_duplicates
                         )
{
  test_insert< Kokkos::Cuda >( num_nodes, num_inserts, num_duplicates, false);
}

void cuda_test_failed_insert(  uint32_t num_nodes )
{
  test_failed_insert< Kokkos::Cuda >( num_nodes );
}

void cuda_test_deep_copy(  uint32_t num_nodes )
{
  test_deep_copy< Kokkos::Cuda >( num_nodes );
}

void cuda_test_vector_combinations(unsigned int size)
{
  test_vector_combinations<int,Kokkos::Cuda>(size);
}

void cuda_test_dualview_combinations(unsigned int size)
{
  test_dualview_combinations<int,Kokkos::Cuda>(size);
}

void cuda_test_bitset()
{
  test_bitset<Kokkos::Cuda>();
}



/*TEST_F( cuda, bitset )
{
  cuda_test_bitset();
}*/

#define CUDA_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat )                                \
  TEST_F( cuda, UnorderedMap_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat##x) {   \
    for (int i=0; i<repeat; ++i)                                                                                \
      cuda_test_insert_##name(num_nodes,num_inserts,num_duplicates);                                            \
  }

#define CUDA_FAILED_INSERT_TEST( num_nodes, repeat )                           \
  TEST_F( cuda, UnorderedMap_failed_insert_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      cuda_test_failed_insert(num_nodes);                                      \
  }

#define CUDA_ASSIGNEMENT_TEST( num_nodes, repeat )                               \
  TEST_F( cuda, UnorderedMap_assignment_operators_##num_nodes##_##repeat##x) {  \
    for (int i=0; i<repeat; ++i)                                                 \
      cuda_test_assignment_operators(num_nodes);                                 \
  }

#define CUDA_DEEP_COPY( num_nodes, repeat )                             \
  TEST_F( cuda, UnorderedMap_deep_copy##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      cuda_test_deep_copy(num_nodes);                     \
  }

#define CUDA_VECTOR_COMBINE_TEST( size )                             \
  TEST_F( cuda, vector_combination##size##x) {       \
      cuda_test_vector_combinations(size);                     \
  }

#define CUDA_DUALVIEW_COMBINE_TEST( size )                             \
  TEST_F( cuda, dualview_combination##size##x) {       \
      cuda_test_dualview_combinations(size);                     \
  }

CUDA_DUALVIEW_COMBINE_TEST( 10 )
CUDA_VECTOR_COMBINE_TEST( 10 )
CUDA_VECTOR_COMBINE_TEST( 3057 )


CUDA_INSERT_TEST(close,               100000, 90000, 100, 500)
CUDA_INSERT_TEST(far,                 100000, 90000, 100, 500)
CUDA_DEEP_COPY( 10000, 1 )
CUDA_FAILED_INSERT_TEST( 10000, 1000 )


#undef CUDA_INSERT_TEST
#undef CUDA_FAILED_INSERT_TEST
#undef CUDA_ASSIGNEMENT_TEST
#undef CUDA_DEEP_COPY
#undef CUDA_VECTOR_COMBINE_TEST
#undef CUDA_DUALVIEW_COMBINE_TEST


TEST_F( cuda , dynamic_view )
{
  typedef TestDynamicView< double , Kokkos::CudaUVMSpace >
    TestDynView ;

  for ( int i = 0 ; i < 10 ; ++i ) {
    TestDynView::run( 100000 + 100 * i );
  }
}


#if defined(KOKKOS_CLASS_LAMBDA)
TEST_F(cuda, ErrorReporterViaLambda)
{
  TestErrorReporter<ErrorReporterDriverUseLambda<Kokkos::Cuda>>();
}
#endif

TEST_F(cuda, ErrorReporter)
{
  TestErrorReporter<ErrorReporterDriver<Kokkos::Cuda>>();
}

}

#else
void KOKKOS_CONTAINERS_UNIT_TESTS_TESTCUDA_PREVENT_EMPTY_LINK_ERROR() {}
#endif  /* #ifdef KOKKOS_ENABLE_CUDA */

