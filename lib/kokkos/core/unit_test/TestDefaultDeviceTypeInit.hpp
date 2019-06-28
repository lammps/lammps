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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_OPENMP
#include <omp.h>
#endif

#if !defined( KOKKOS_ENABLE_CUDA ) || defined( __CUDACC__ )

namespace Test {

namespace Impl {

char** init_kokkos_args( bool do_threads, bool do_numa, bool do_device, bool do_other, int & nargs, Kokkos::InitArguments & init_args ) {
  nargs = ( do_threads ? 1 : 0 ) +
          ( do_numa ? 1 : 0 ) +
          ( do_device ? 1 : 0 ) +
          ( do_other ? 4 : 0 );

  char** args_kokkos = new char*[nargs];
  for ( int i = 0; i < nargs; i++ ) {
    args_kokkos[i] = new char[20];
  }

  int threads_idx = do_other ? 1 : 0;
  int numa_idx = ( do_other ? 3 : 0 ) + ( do_threads ? 1 : 0 );
  int device_idx = ( do_other ? 3 : 0 ) + ( do_threads ? 1 : 0 ) + ( do_numa ? 1 : 0 );

  if ( do_threads ) {
    int nthreads = 3;

#ifdef KOKKOS_ENABLE_OPENMP
    if ( omp_get_max_threads() < 3 )
      nthreads = omp_get_max_threads();
#endif

    if ( Kokkos::hwloc::available() ) {
      if ( Kokkos::hwloc::get_available_threads_per_core() < 3 )
        nthreads =   Kokkos::hwloc::get_available_threads_per_core()
                   * Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::Serial, Kokkos::DefaultExecutionSpace >::value ||
         std::is_same< Kokkos::Serial, Kokkos::DefaultHostExecutionSpace >::value ) {
      nthreads = 1;
    }
#endif

    init_args.num_threads = nthreads;
    sprintf( args_kokkos[threads_idx], "--threads=%i", nthreads );
  }

  if ( do_numa ) {
    int numa = 1;
    if ( Kokkos::hwloc::available() ) {
      numa = Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::Serial, Kokkos::DefaultExecutionSpace >::value ||
         std::is_same< Kokkos::Serial, Kokkos::DefaultHostExecutionSpace >::value ) {
      numa = 1;
    }
#endif

    init_args.num_numa = numa;
    sprintf( args_kokkos[numa_idx], "--numa=%i", numa );
  }

  if ( do_device ) {
    init_args.device_id = 0;
    sprintf( args_kokkos[device_idx], "--device=%i", 0 );
  }

  if ( do_other ) {
    sprintf( args_kokkos[0], "--dummyarg=1" );
    sprintf( args_kokkos[ threads_idx + ( do_threads ? 1 : 0 ) ], "--dummy2arg" );
    sprintf( args_kokkos[ threads_idx + ( do_threads ? 1 : 0 ) + 1 ], "dummy3arg" );
    sprintf( args_kokkos[ device_idx + ( do_device ? 1 : 0 ) ], "dummy4arg=1" );
  }

  return args_kokkos;
}

Kokkos::InitArguments init_initstruct( bool do_threads, bool do_numa, bool do_device ) {
  Kokkos::InitArguments args;

  if ( do_threads ) {
    int nthreads = 3;

#ifdef KOKKOS_ENABLE_OPENMP
    if ( omp_get_max_threads() < 3 ) {
      nthreads = omp_get_max_threads();
    }
#endif

    if ( Kokkos::hwloc::available() ) {
      if ( Kokkos::hwloc::get_available_threads_per_core() < 3 ) {
        nthreads =   Kokkos::hwloc::get_available_threads_per_core()
                   * Kokkos::hwloc::get_available_numa_count();
      }
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::Serial, Kokkos::DefaultExecutionSpace >::value ||
         std::is_same< Kokkos::Serial, Kokkos::DefaultHostExecutionSpace >::value ) {
      nthreads = 1;
    }
#endif

    args.num_threads = nthreads;
  }

  if ( do_numa ) {
    int numa = 1;
    if ( Kokkos::hwloc::available() ) {
      numa = Kokkos::hwloc::get_available_numa_count();
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::Serial, Kokkos::DefaultExecutionSpace >::value ||
         std::is_same< Kokkos::Serial, Kokkos::DefaultHostExecutionSpace >::value ) {
      numa = 1;
    }
#endif

    args.num_numa = numa;
  }

  if ( do_device ) {
    args.device_id = 0;
  }

  return args;
}

void check_correct_initialization( const Kokkos::InitArguments & argstruct ) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  ASSERT_EQ( Kokkos::DefaultExecutionSpace::is_initialized(), 1 );
  ASSERT_EQ( Kokkos::HostSpace::execution_space::is_initialized(), 1 );
#else
  ASSERT_EQ( Kokkos::DefaultExecutionSpace::impl_is_initialized(), 1 );
  ASSERT_EQ( Kokkos::HostSpace::execution_space::impl_is_initialized(), 1 );
#endif

  // Figure out the number of threads the HostSpace ExecutionSpace should have initialized to.
  int expected_nthreads = argstruct.num_threads;

#ifdef KOKKOS_ENABLE_OPENMP
  if ( std::is_same< Kokkos::HostSpace::execution_space, Kokkos::OpenMP >::value ) {
    // use openmp default num threads
    if ( expected_nthreads < 0 || ( expected_nthreads == 0 && !Kokkos::hwloc::available() ) ) {
      expected_nthreads = omp_get_max_threads();
    }
    // use hwloc if available
    else if ( expected_nthreads == 0 && Kokkos::hwloc::available() ) {
      expected_nthreads = Kokkos::hwloc::get_available_numa_count()
                        * Kokkos::hwloc::get_available_cores_per_numa()
                        * Kokkos::hwloc::get_available_threads_per_core();
    }
  }
#endif

  if ( expected_nthreads < 1 ) {
    if ( Kokkos::hwloc::available() ) {
      expected_nthreads = Kokkos::hwloc::get_available_numa_count()
                        * Kokkos::hwloc::get_available_cores_per_numa()
                        * Kokkos::hwloc::get_available_threads_per_core();
    }
    else {
        expected_nthreads = 1;
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Serial >::value ||
         std::is_same< Kokkos::DefaultHostExecutionSpace, Kokkos::Serial >::value ) {
      expected_nthreads = 1;
    }
#endif

#ifdef KOKKOS_ENABLE_HPX
    // HPX uses all cores on machine by default. Skip this test.
    if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Experimental::HPX >::value ||
         std::is_same< Kokkos::DefaultHostExecutionSpace, Kokkos::Experimental::HPX >::value ) {
      return;
    }
#endif
  }

  int expected_numa = argstruct.num_numa;

  if ( expected_numa < 1 ) {
    if ( Kokkos::hwloc::available() ) {
      expected_numa = Kokkos::hwloc::get_available_numa_count();
    }
    else {
      expected_numa = 1;
    }

#ifdef KOKKOS_ENABLE_SERIAL
    if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Serial >::value ||
         std::is_same< Kokkos::DefaultHostExecutionSpace, Kokkos::Serial >::value )
      expected_numa = 1;
#endif
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  ASSERT_EQ( Kokkos::HostSpace::execution_space::thread_pool_size(), expected_nthreads );
#else
  ASSERT_EQ( Kokkos::HostSpace::execution_space::impl_thread_pool_size(), expected_nthreads );
#endif

#ifdef KOKKOS_ENABLE_CUDA
  if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Cuda >::value ) {
    int device;
    cudaGetDevice( &device );

    int expected_device = argstruct.device_id;
    if ( argstruct.device_id < 0 ) {
      expected_device = 0;
    }

    ASSERT_EQ( expected_device, device );
  }
#endif
}

// TODO: Add check whether correct number of threads are actually started.
void test_no_arguments() {
  Kokkos::initialize();
  check_correct_initialization( Kokkos::InitArguments() );
  Kokkos::finalize();
}

void test_commandline_args( int nargs, char** args, const Kokkos::InitArguments & argstruct ) {
  Kokkos::initialize( nargs, args );
  check_correct_initialization( argstruct );
  Kokkos::finalize();
}

void test_initstruct_args( const Kokkos::InitArguments & args ) {
  Kokkos::initialize( args );
  check_correct_initialization( args );
  Kokkos::finalize();
}

} // namespace Impl

class defaultdevicetypeinit : public ::testing::Test {
protected:
  static void SetUpTestCase() {}

  static void TearDownTestCase() {}
};

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_01
TEST_F( defaultdevicetypeinit, no_args )
{
  Impl::test_no_arguments();
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_02
TEST_F( defaultdevicetypeinit, commandline_args_empty )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( false, false, false, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_03
TEST_F( defaultdevicetypeinit, commandline_args_other )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( false, false, false, true, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_04
TEST_F( defaultdevicetypeinit, commandline_args_nthreads )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( true, false, false, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_05
TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( true, true, false, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_06
TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa_device )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( true, true, true, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_07
TEST_F( defaultdevicetypeinit, commandline_args_nthreads_device )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( true, false, true, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_08
TEST_F( defaultdevicetypeinit, commandline_args_numa_device )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( false, true, true, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_09
TEST_F( defaultdevicetypeinit, commandline_args_device )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( false, false, true, false, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_10
TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa_device_other )
{
  Kokkos::InitArguments argstruct;
  int nargs = 0;
  char** args = Impl::init_kokkos_args( true, true, true, true, nargs, argstruct );
  Impl::test_commandline_args( nargs, args, argstruct );

  for ( int i = 0; i < nargs; i++ ) {
    delete [] args[i];
  }
  delete [] args;
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_11
TEST_F( defaultdevicetypeinit, initstruct_default )
{
  Kokkos::InitArguments args;
  Impl::test_initstruct_args( args );
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_12
TEST_F( defaultdevicetypeinit, initstruct_nthreads )
{
  Kokkos::InitArguments args = Impl::init_initstruct( true, false, false );
  Impl::test_initstruct_args( args );
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_13
TEST_F( defaultdevicetypeinit, initstruct_nthreads_numa )
{
  Kokkos::InitArguments args = Impl::init_initstruct( true, true, false );
  Impl::test_initstruct_args( args );
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_14
TEST_F( defaultdevicetypeinit, initstruct_device )
{
  Kokkos::InitArguments args = Impl::init_initstruct( false, false, true );
  Impl::test_initstruct_args( args );
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_15
TEST_F( defaultdevicetypeinit, initstruct_nthreads_device )
{
  Kokkos::InitArguments args = Impl::init_initstruct( true, false, true );
  Impl::test_initstruct_args( args );
}
#endif

#ifdef KOKKOS_DEFAULTDEVICETYPE_INIT_TEST_16
TEST_F( defaultdevicetypeinit, initstruct_nthreads_numa_device )
{
  Kokkos::InitArguments args = Impl::init_initstruct( true, true, true );
  Impl::test_initstruct_args( args );
}
#endif

} // namespace Test

#endif
