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
#ifdef KOKKOS_ENABLE_THREADS

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <Kokkos_Bitset.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <Kokkos_Vector.hpp>
#include <iomanip>


//----------------------------------------------------------------------------
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

namespace Test {

class threads : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;

    unsigned num_threads = 4;

    if (Kokkos::hwloc::available()) {
      num_threads = Kokkos::hwloc::get_available_numa_count()
                    * Kokkos::hwloc::get_available_cores_per_numa()
                 // * Kokkos::hwloc::get_available_threads_per_core()
                    ;

    }

    std::cout << "Threads: " << num_threads << std::endl;

    Kokkos::Threads::initialize( num_threads );
  }

  static void TearDownTestCase()
  {
    Kokkos::Threads::finalize();
  }
};

TEST_F( threads , dyn_view_api) {
  TestDynViewAPI< double , Kokkos::Threads >();
}

TEST_F( threads, viewctorprop_embedded_dim ) {
  TestViewCtorProp_EmbeddedDim< Kokkos::Threads >::test_vcpt( 2, 3 );
}

TEST_F( threads , staticcrsgraph )
{
  TestStaticCrsGraph::run_test_graph< Kokkos::Threads >();
  TestStaticCrsGraph::run_test_graph2< Kokkos::Threads >();
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(1, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(1, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(1, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(1, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(3, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(3, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(3, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(3, 100000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(75, 0);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(75, 1000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(75, 10000);
  TestStaticCrsGraph::run_test_graph3< Kokkos::Threads >(75, 100000);
}

/*TEST_F( threads, bitset )
{
  test_bitset<Kokkos::Threads>();
}*/

#define THREADS_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat, near )                                \
  TEST_F( threads, UnorderedMap_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat##x) {   \
    for (int i=0; i<repeat; ++i)                                                                                \
      test_insert<Kokkos::Threads>(num_nodes,num_inserts,num_duplicates, near);                                   \
  }

#define THREADS_FAILED_INSERT_TEST( num_nodes, repeat )                            \
  TEST_F( threads, UnorderedMap_failed_insert_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_failed_insert<Kokkos::Threads>(num_nodes);                             \
  }

#define THREADS_ASSIGNEMENT_TEST( num_nodes, repeat )                             \
  TEST_F( threads, UnorderedMap_assignment_operators_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_assignement_operators<Kokkos::Threads>(num_nodes);                     \
  }

#define THREADS_DEEP_COPY( num_nodes, repeat )                             \
  TEST_F( threads, UnorderedMap_deep_copy##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_deep_copy<Kokkos::Threads>(num_nodes);                     \
  }

#define THREADS_VECTOR_COMBINE_TEST( size )                             \
  TEST_F( threads, vector_combination##size##x) {       \
      test_vector_combinations<int,Kokkos::Threads>(size);                     \
  }

#define THREADS_DUALVIEW_COMBINE_TEST( size )                             \
  TEST_F( threads, dualview_combination##size##x) {       \
      test_dualview_combinations<int,Kokkos::Threads>(size);                     \
  }

THREADS_INSERT_TEST(far, 100000, 90000, 100, 500, false)
THREADS_FAILED_INSERT_TEST( 10000, 1000 )
THREADS_DEEP_COPY( 10000, 1 )

THREADS_VECTOR_COMBINE_TEST( 10 )
THREADS_VECTOR_COMBINE_TEST( 3057 )
THREADS_DUALVIEW_COMBINE_TEST( 10 )


#undef THREADS_INSERT_TEST
#undef THREADS_FAILED_INSERT_TEST
#undef THREADS_ASSIGNEMENT_TEST
#undef THREADS_DEEP_COPY
#undef THREADS_VECTOR_COMBINE_TEST
#undef THREADS_DUALVIEW_COMBINE_TEST


TEST_F( threads , dynamic_view )
{
  typedef TestDynamicView< double , Kokkos::Threads >
    TestDynView ;

  for ( int i = 0 ; i < 10 ; ++i ) {
    TestDynView::run( 100000 + 100 * i );
  }
}


#if defined(KOKKOS_CLASS_LAMBDA)
TEST_F(threads, ErrorReporterViaLambda)
{
  TestErrorReporter<ErrorReporterDriverUseLambda<Kokkos::Threads>>();
}
#endif

TEST_F(threads, ErrorReporter)
{
  TestErrorReporter<ErrorReporterDriver<Kokkos::Threads>>();
}

} // namespace Test

#else
void KOKKOS_CONTAINERS_UNIT_TESTS_TESTTHREADS_PREVENT_EMPTY_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_THREADS ) */

