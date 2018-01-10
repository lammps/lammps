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
#ifdef KOKKOS_ENABLE_ROCM

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

class rocm : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Experimental::ROCm::initialize( Kokkos::Experimental::ROCm::SelectDevice(0) );
  }
  static void TearDownTestCase()
  {
    Kokkos::Experimental::ROCm::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }
};

#if !defined(KOKKOS_ENABLE_ROCM)
//issue 964
TEST_F( rocm , dyn_view_api) {
  TestDynViewAPI< double , Kokkos::Experimental::ROCm >();
}
#endif 

TEST_F( rocm, viewctorprop_embedded_dim ) {
  TestViewCtorProp_EmbeddedDim< Kokkos::Experimental::ROCm >::test_vcpt( 2, 3 );
}

TEST_F( rocm , staticcrsgraph )
{
  TestStaticCrsGraph::run_test_graph< Kokkos::Experimental::ROCm >();
  TestStaticCrsGraph::run_test_graph2< Kokkos::Experimental::ROCm >();
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(1, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(1, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(1, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(1, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(3, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(3, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(3, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(3, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(75, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(75, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(75, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Experimental::ROCm >(75, 100000);
}


#if !defined(KOKKOS_ENABLE_ROCM)
// issue 1089
// same as 130203 (MemPool, static member function link issue
void rocm_test_insert_close(  uint32_t num_nodes
                            , uint32_t num_inserts
                            , uint32_t num_duplicates
                           )
{
  test_insert< Kokkos::Experimental::ROCm >( num_nodes, num_inserts, num_duplicates, true);
}

// hcc link error , Referencing function in another module!
void rocm_test_insert_far(  uint32_t num_nodes
                          , uint32_t num_inserts
                          , uint32_t num_duplicates
                         )
{
  test_insert< Kokkos::Experimental::ROCm >( num_nodes, num_inserts, num_duplicates, false);
}

void rocm_test_failed_insert(  uint32_t num_nodes )
{
  test_failed_insert< Kokkos::Experimental::ROCm >( num_nodes );
}

void rocm_test_deep_copy(  uint32_t num_nodes )
{
  test_deep_copy< Kokkos::Experimental::ROCm >( num_nodes );
}

void rocm_test_vector_combinations(unsigned int size)
{
  test_vector_combinations<int,Kokkos::Experimental::ROCm>(size);
}

void rocm_test_dualview_combinations(unsigned int size)
{
  test_dualview_combinations<int,Kokkos::Experimental::ROCm>(size);
}

void rocm_test_bitset()
{
  test_bitset<Kokkos::Experimental::ROCm>();
}



/*TEST_F( rocm, bitset )
{
  rocm_test_bitset();
}*/

#define ROCM_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat )                                \
  TEST_F( rocm, UnorderedMap_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat##x) {   \
    for (int i=0; i<repeat; ++i)                                                                                \
      rocm_test_insert_##name(num_nodes,num_inserts,num_duplicates);                                            \
  }

#define ROCM_FAILED_INSERT_TEST( num_nodes, repeat )                           \
  TEST_F( rocm, UnorderedMap_failed_insert_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      rocm_test_failed_insert(num_nodes);                                      \
  }

#define ROCM_ASSIGNEMENT_TEST( num_nodes, repeat )                               \
  TEST_F( rocm, UnorderedMap_assignment_operators_##num_nodes##_##repeat##x) {  \
    for (int i=0; i<repeat; ++i)                                                 \
      rocm_test_assignment_operators(num_nodes);                                 \
  }

#define ROCM_DEEP_COPY( num_nodes, repeat )                             \
  TEST_F( rocm, UnorderedMap_deep_copy##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      rocm_test_deep_copy(num_nodes);                     \
  }

#define ROCM_VECTOR_COMBINE_TEST( size )                             \
  TEST_F( rocm, vector_combination##size##x) {       \
      rocm_test_vector_combinations(size);                     \
  }

#define ROCM_DUALVIEW_COMBINE_TEST( size )                             \
  TEST_F( rocm, dualview_combination##size##x) {       \
      rocm_test_dualview_combinations(size);                     \
  }

//ROCM_DUALVIEW_COMBINE_TEST( 10 )
//ROCM_VECTOR_COMBINE_TEST( 10 )
//ROCM_VECTOR_COMBINE_TEST( 3057 )


//ROCM_INSERT_TEST(close,               100000, 90000, 100, 500)
//ROCM_INSERT_TEST(far,                 100000, 90000, 100, 500)
//ROCM_DEEP_COPY( 10000, 1 )
//ROCM_FAILED_INSERT_TEST( 10000, 1000 )


#undef ROCM_INSERT_TEST
#undef ROCM_FAILED_INSERT_TEST
#undef ROCM_ASSIGNEMENT_TEST
#undef ROCM_DEEP_COPY
#undef ROCM_VECTOR_COMBINE_TEST
#undef ROCM_DUALVIEW_COMBINE_TEST


#endif
#if !defined(KOKKOS_ENABLE_ROCM)
//static member function issue 
TEST_F( rocm , dynamic_view )
{
//  typedef TestDynamicView< double , Kokkos::ROCmUVMSpace >
  typedef TestDynamicView< double , Kokkos::Experimental::ROCmSpace >
    TestDynView ;

  for ( int i = 0 ; i < 10 ; ++i ) {
    TestDynView::run( 100000 + 100 * i );
  }
}
#endif


#if defined(KOKKOS_CLASS_LAMBDA)
TEST_F(rocm, ErrorReporterViaLambda)
{
  TestErrorReporter<ErrorReporterDriverUseLambda<Kokkos::Experimental::ROCm>>();
}
#endif

TEST_F(rocm, ErrorReporter)
{
  TestErrorReporter<ErrorReporterDriver<Kokkos::Experimental::ROCm>>();
}

}

#else
void KOKKOS_CONTAINERS_UNIT_TESTS_TESTROCM_PREVENT_EMPTY_LINK_ERROR() {}
#endif  /* #ifdef KOKKOS_ENABLE_ROCM */

