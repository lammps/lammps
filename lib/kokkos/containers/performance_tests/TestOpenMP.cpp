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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_OPENMP )

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <Kokkos_UnorderedMap.hpp>

#include <TestGlobal2LocalIds.hpp>
#include <TestUnorderedMapPerformance.hpp>

#include <TestDynRankView.hpp>
#include <TestScatterView.hpp>

#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>


namespace Performance {

class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;

    Kokkos::initialize();
    Kokkos::OpenMP::print_configuration( std::cout );
  }

  static void TearDownTestCase()
  {
    Kokkos::finalize();
  }
};

TEST_F( openmp, dynrankview_perf )
{
  std::cout << "OpenMP" << std::endl;
  std::cout << " DynRankView vs View: Initialization Only " << std::endl;
  test_dynrankview_op_perf<Kokkos::OpenMP>( 8192 );
}

TEST_F( openmp, global_2_local)
{
  std::cout << "OpenMP" << std::endl;
  std::cout << "size, create, generate, fill, find" << std::endl;
  for (unsigned i=Performance::begin_id_size; i<=Performance::end_id_size; i *= Performance::id_step)
    test_global_to_local_ids<Kokkos::OpenMP>(i);
}

TEST_F( openmp, unordered_map_performance_near)
{
  unsigned num_openmp = 4;
  if (Kokkos::hwloc::available()) {
    num_openmp = Kokkos::hwloc::get_available_numa_count() *
                  Kokkos::hwloc::get_available_cores_per_numa() *
                  Kokkos::hwloc::get_available_threads_per_core();

  }
  std::ostringstream base_file_name;
  base_file_name << "openmp-" << num_openmp << "-near";
  Perf::run_performance_tests<Kokkos::OpenMP,true>(base_file_name.str());
}

TEST_F( openmp, unordered_map_performance_far)
{
  unsigned num_openmp = 4;
  if (Kokkos::hwloc::available()) {
    num_openmp = Kokkos::hwloc::get_available_numa_count() *
                  Kokkos::hwloc::get_available_cores_per_numa() *
                  Kokkos::hwloc::get_available_threads_per_core();

  }
  std::ostringstream base_file_name;
  base_file_name << "openmp-" << num_openmp << "-far";
  Perf::run_performance_tests<Kokkos::OpenMP,false>(base_file_name.str());
}

TEST_F( openmp, scatter_view)
{
  std::cout << "ScatterView data-duplicated test:\n";
  Perf::test_scatter_view<Kokkos::OpenMP, Kokkos::LayoutRight,
    Kokkos::Experimental::ScatterDuplicated,
    Kokkos::Experimental::ScatterNonAtomic>(10, 1000 * 1000);
//std::cout << "ScatterView atomics test:\n";
//Perf::test_scatter_view<Kokkos::OpenMP, Kokkos::LayoutRight,
//  Kokkos::Experimental::ScatterNonDuplicated,
//  Kokkos::Experimental::ScatterAtomic>(10, 1000 * 1000);
}

} // namespace test
#else
void KOKKOS_CONTAINERS_PERFORMANCE_TESTS_TESTOPENMP_PREVENT_EMPTY_LINK_ERROR() {}
#endif

