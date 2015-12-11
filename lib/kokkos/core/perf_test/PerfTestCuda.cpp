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

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#if defined( KOKKOS_HAVE_CUDA )

#include <impl/Kokkos_Timer.hpp>

#include <PerfTestHexGrad.hpp>
#include <PerfTestBlasKernels.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>


namespace Test {

class cuda : public ::testing::Test {
  protected:
    static void SetUpTestCase() {
      Kokkos::HostSpace::execution_space::initialize();
      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
    }
    static void TearDownTestCase() {
      Kokkos::Cuda::finalize();
      Kokkos::HostSpace::execution_space::finalize();
    }
};

TEST_F( cuda, hexgrad )
{
  EXPECT_NO_THROW( run_test_hexgrad< Kokkos::Cuda >( 10 , 20, "Kokkos::Cuda" ) );
}

TEST_F( cuda, gramschmidt )
{
  EXPECT_NO_THROW( run_test_gramschmidt< Kokkos::Cuda >( 10 , 20, "Kokkos::Cuda" ) );
}

namespace {

template <typename T>
struct TextureFetch
{
  typedef Kokkos::View< T *, Kokkos::CudaSpace> array_type;
  typedef Kokkos::View< const T *, Kokkos::CudaSpace, Kokkos::MemoryRandomAccess> const_array_type;
  typedef Kokkos::View< int *, Kokkos::CudaSpace> index_array_type;
  typedef Kokkos::View< const int *, Kokkos::CudaSpace> const_index_array_type;

  struct FillArray
  {
    array_type m_array;
    FillArray( const array_type & array )
      : m_array(array)
    {}

    void apply() const
    {
      Kokkos::parallel_for( Kokkos::RangePolicy<Kokkos::Cuda,int>(0,m_array.dimension_0()), *this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int i) const { m_array(i) = i; }
  };

  struct RandomIndexes
  {
    index_array_type m_indexes;
    typename index_array_type::HostMirror m_host_indexes;
    RandomIndexes( const index_array_type & indexes)
      : m_indexes(indexes)
      , m_host_indexes(Kokkos::create_mirror(m_indexes))
    {}

    void apply() const
    {
      Kokkos::parallel_for( Kokkos::RangePolicy<Kokkos::HostSpace::execution_space,int>(0,m_host_indexes.dimension_0()), *this);
      //random shuffle
      Kokkos::HostSpace::execution_space::fence();
      std::random_shuffle(m_host_indexes.ptr_on_device(), m_host_indexes.ptr_on_device() + m_host_indexes.dimension_0());
      Kokkos::deep_copy(m_indexes,m_host_indexes);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int i) const { m_host_indexes(i) = i; }
  };

  struct RandomReduce
  {
    const_array_type       m_array;
    const_index_array_type m_indexes;
    RandomReduce( const const_array_type & array, const const_index_array_type & indexes)
      : m_array(array)
      , m_indexes(indexes)
    {}

    void apply(T & reduce) const
    {
      Kokkos::parallel_reduce( Kokkos::RangePolicy<Kokkos::Cuda,int>(0,m_array.dimension_0()), *this, reduce);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int i, T & reduce) const
    { reduce += m_array(m_indexes(i)); }
  };

  static void run(int size, double & reduce_time, T &reduce)
  {
    array_type array("array",size);
    index_array_type indexes("indexes",size);

    { FillArray f(array); f.apply(); }
    { RandomIndexes f(indexes); f.apply(); }

    Kokkos::Cuda::fence();

    Kokkos::Impl::Timer timer;
    for (int j=0; j<10; ++j) {
      RandomReduce f(array,indexes);
      f.apply(reduce);
    }
    Kokkos::Cuda::fence();
    reduce_time = timer.seconds();
  }
};

} // unnamed namespace

TEST_F( cuda, texture_double )
{
  printf("Random reduce of double through texture fetch\n");
  for (int i=1; i<=26; ++i) {
    int size = 1<<i;
    double time = 0;
    double reduce = 0;
    TextureFetch<double>::run(size,time,reduce);
    printf("   time = %1.3e   size = 2^%d\n", time, i);
  }
}

} // namespace Test

#endif /* #if defined( KOKKOS_HAVE_CUDA ) */

