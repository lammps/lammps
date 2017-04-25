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

#include <cuda/TestCuda.hpp>

namespace Test {

TEST_F( cuda, long_reduce )
{
  TestReduce< long, Kokkos::Cuda >( 0 );
  TestReduce< long, Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, double_reduce )
{
  TestReduce< double, Kokkos::Cuda >( 0 );
  TestReduce< double, Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, long_reduce_dynamic )
{
  TestReduceDynamic< long, Kokkos::Cuda >( 0 );
  TestReduceDynamic< long, Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, double_reduce_dynamic )
{
  TestReduceDynamic< double, Kokkos::Cuda >( 0 );
  TestReduceDynamic< double, Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, long_reduce_dynamic_view )
{
  TestReduceDynamicView< long, Kokkos::Cuda >( 0 );
  TestReduceDynamicView< long, Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, scan )
{
  TestScan< Kokkos::Cuda >::test_range( 1, 1000 );
  TestScan< Kokkos::Cuda >( 0 );
  TestScan< Kokkos::Cuda >( 100000 );
  TestScan< Kokkos::Cuda >( 10000000 );
  Kokkos::Cuda::fence();
}

#if 0
TEST_F( cuda, scan_small )
{
  typedef TestScan< Kokkos::Cuda, Kokkos::Impl::CudaExecUseScanSmall > TestScanFunctor;

  for ( int i = 0; i < 1000; ++i ) {
    TestScanFunctor( 10 );
    TestScanFunctor( 10000 );
  }
  TestScanFunctor( 1000000 );
  TestScanFunctor( 10000000 );

  Kokkos::Cuda::fence();
}
#endif

TEST_F( cuda, team_scan )
{
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 0 );
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 0 );
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 10 );
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 10 );
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 10000 );
  TestScanTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 10000 );
}

TEST_F( cuda, team_long_reduce )
{
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 0 );
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 0 );
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< long, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
}

TEST_F( cuda, team_double_reduce )
{
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 0 );
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 0 );
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< double, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
}

TEST_F( cuda, reduction_deduction )
{
  TestCXX11::test_reduction_deduction< Kokkos::Cuda >();
}

} // namespace Test
