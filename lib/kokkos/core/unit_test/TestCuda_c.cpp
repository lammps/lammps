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

#include <gtest/gtest.h>

#include <iostream>

#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------

#include <Cuda/Kokkos_Cuda_TaskPolicy.hpp>
#include <impl/Kokkos_ViewTileLeft.hpp>
#include <TestTile.hpp>

//----------------------------------------------------------------------------

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>
#include <TestAtomicOperations.hpp>

#include <TestViewAPI.hpp>
#include <TestViewSubview.hpp>
#include <TestViewOfClass.hpp>

#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestMemoryPool.hpp>
#include <TestTeamVector.hpp>
#include <TestTemplateMetaFunctions.hpp>
#include <TestCXX11Deduction.hpp>

#include <TestTaskPolicy.hpp>
#include <TestPolicyConstruction.hpp>

//----------------------------------------------------------------------------

class cuda : public ::testing::Test {
protected:
  static void SetUpTestCase();
  static void TearDownTestCase();
};

//----------------------------------------------------------------------------

namespace Test {

TEST_F( cuda, atomic )
{
  const int loop_count = 1e3 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Cuda>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Cuda>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Cuda>(100,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<TestAtomic::SuperScalar<4> ,Kokkos::Cuda>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<TestAtomic::SuperScalar<4> ,Kokkos::Cuda>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<TestAtomic::SuperScalar<4> ,Kokkos::Cuda>(100,3) ) );

}

TEST_F( cuda , atomic_operations )
{
  const int start = 1; //Avoid zero for division
  const int end = 11;
  for (int i = start; i < end; ++i)
  {
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<int,Kokkos::Cuda>(start, end-i, 9 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned int,Kokkos::Cuda>(start, end-i, 9 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long int,Kokkos::Cuda>(start, end-i, 9 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<unsigned long int,Kokkos::Cuda>(start, end-i, 9 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType<long long int,Kokkos::Cuda>(start, end-i, 9 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<double,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<double,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<double,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<double,Kokkos::Cuda>(start, end-i, 4 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<float,Kokkos::Cuda>(start, end-i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<float,Kokkos::Cuda>(start, end-i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<float,Kokkos::Cuda>(start, end-i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType<float,Kokkos::Cuda>(start, end-i, 4 ) ) );
  }

}

//----------------------------------------------------------------------------

TEST_F( cuda, tile_layout)
{
  TestTile::test< Kokkos::Cuda , 1 , 1 >( 1 , 1 );
  TestTile::test< Kokkos::Cuda , 1 , 1 >( 2 , 3 );
  TestTile::test< Kokkos::Cuda , 1 , 1 >( 9 , 10 );

  TestTile::test< Kokkos::Cuda , 2 , 2 >( 1 , 1 );
  TestTile::test< Kokkos::Cuda , 2 , 2 >( 2 , 3 );
  TestTile::test< Kokkos::Cuda , 2 , 2 >( 4 , 4 );
  TestTile::test< Kokkos::Cuda , 2 , 2 >( 9 , 9 );

  TestTile::test< Kokkos::Cuda , 2 , 4 >( 9 , 9 );
  TestTile::test< Kokkos::Cuda , 4 , 4 >( 9 , 9 );

  TestTile::test< Kokkos::Cuda , 4 , 4 >( 1 , 1 );
  TestTile::test< Kokkos::Cuda , 4 , 4 >( 4 , 4 );
  TestTile::test< Kokkos::Cuda , 4 , 4 >( 9 , 9 );
  TestTile::test< Kokkos::Cuda , 4 , 4 >( 9 , 11 );

  TestTile::test< Kokkos::Cuda , 8 , 8 >( 1 , 1 );
  TestTile::test< Kokkos::Cuda , 8 , 8 >( 4 , 4 );
  TestTile::test< Kokkos::Cuda , 8 , 8 >( 9 , 9 );
  TestTile::test< Kokkos::Cuda , 8 , 8 >( 9 , 11 );
}

TEST_F( cuda , view_aggregate )
{
  TestViewAggregate< Kokkos::Cuda >();
  TestViewAggregateReduction< Kokkos::Cuda >();
}

TEST_F( cuda , scan )
{
  TestScan< Kokkos::Cuda >::test_range( 1 , 1000 );
  TestScan< Kokkos::Cuda >( 1000000 );
  TestScan< Kokkos::Cuda >( 10000000 );

  TestScan< Kokkos::Cuda >( 0 );
  TestScan< Kokkos::Cuda >( 0 , 0 );

  Kokkos::Cuda::fence();
}

TEST_F( cuda , team_scan )
{
  TestScanTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 10 );
  TestScanTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 10 );
  TestScanTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 10000 );
  TestScanTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 10000 );
}

TEST_F( cuda , memory_pool )
{
//  typedef Kokkos::CudaUVMSpace  device_type;
  typedef Kokkos::Cuda          device_type;

  bool val = TestMemoryPool::test_mempool< device_type >( 128, 128000000 );
  ASSERT_TRUE( val );

  Kokkos::Cuda::fence();

  TestMemoryPool::test_mempool2< device_type >( 64, 4, 100000, 200000 );

  Kokkos::Cuda::fence();

  TestMemoryPool::test_memory_exhaustion< Kokkos::Cuda >();

  Kokkos::Cuda::fence();
}

}

//----------------------------------------------------------------------------

TEST_F( cuda , template_meta_functions )
{
  TestTemplateMetaFunctions<int, Kokkos::Cuda >();
}

//----------------------------------------------------------------------------

namespace Test {

TEST_F( cuda , reduction_deduction )
{
  TestCXX11::test_reduction_deduction< Kokkos::Cuda >();
}

TEST_F( cuda , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(0) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(3) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(4) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(5) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(6) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(7) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(8) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(9) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >(10) ) );
}

TEST_F( cuda, triple_nested_parallelism )
{
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048 , 32 , 32 );
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048 , 32 , 16 );
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048 , 16 , 16 );
}

}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_TASKPOLICY )

TEST_F( cuda , task_fib )
{
  for ( int i = 0 ; i < 25 ; ++i ) {
    TestTaskPolicy::TestFib< Kokkos::Cuda >::run(i, (i+1)*1000000 );
  }
}

TEST_F( cuda , task_depend )
{
  for ( int i = 0 ; i < 25 ; ++i ) {
    TestTaskPolicy::TestTaskDependence< Kokkos::Cuda >::run(i);
  }
}

TEST_F( cuda , task_team )
{
  //TestTaskPolicy::TestTaskTeam< Kokkos::Cuda >::run(1000);
  TestTaskPolicy::TestTaskTeam< Kokkos::Cuda >::run(104);
  TestTaskPolicy::TestTaskTeamValue< Kokkos::Cuda >::run(1000);
}

//----------------------------------------------------------------------------

TEST_F( cuda , old_task_policy )
{
  TestTaskPolicy::test_task_dep< Kokkos::Cuda >( 10 );

  for ( long i = 0 ; i < 15 ; ++i ) {
      // printf("TestTaskPolicy::test_fib< Kokkos::Cuda >(%d);\n",i);
    TestTaskPolicy::test_fib< Kokkos::Cuda >(i,4096);
  }
  for ( long i = 0 ; i < 35 ; ++i ) {
      // printf("TestTaskPolicy::test_fib2< Kokkos::Cuda >(%d);\n",i);
    TestTaskPolicy::test_fib2< Kokkos::Cuda >(i,4096);
  }
}

TEST_F( cuda , old_task_team )
{
  TestTaskPolicy::test_task_team< Kokkos::Cuda >(1000);
}

TEST_F( cuda , old_task_latch )
{
  TestTaskPolicy::test_latch< Kokkos::Cuda >(10);
  TestTaskPolicy::test_latch< Kokkos::Cuda >(1000);
}

#endif // #if defined( KOKKOS_ENABLE_TASKPOLICY )

