
/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#include <openmp/TestOpenMP_Category.hpp>
#include <Kokkos_Core.hpp>

#include <mutex>

namespace Test {

TEST(openmp, partition_master) {
  using Mutex = Kokkos::Experimental::MasterLock<Kokkos::OpenMP>;

  Mutex mtx;
  int errors = 0;

  auto master = [&errors, &mtx](int /*partition_id*/, int /*num_partitions*/) {
    const int pool_size = Kokkos::OpenMP::impl_thread_pool_size();

    {
      std::unique_lock<Mutex> lock(mtx);
      if (Kokkos::OpenMP::in_parallel()) {
        ++errors;
      }
      if (Kokkos::OpenMP::impl_thread_pool_rank() != 0) {
        ++errors;
      }
    }

    {
      int local_errors = 0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<Kokkos::OpenMP>(0, 1000),
          [pool_size](const int, int& errs) {
            if (Kokkos::OpenMP::impl_thread_pool_size() != pool_size) {
              ++errs;
            }
          },
          local_errors);
      Kokkos::atomic_add(&errors, local_errors);
    }

    Kokkos::Experimental::UniqueToken<Kokkos::OpenMP> token;

    Kokkos::View<int*, Kokkos::OpenMP> count("", token.size());

    Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0, 1000),
                         [=](const int) {
                           int i = token.acquire();
                           ++count[i];
                           token.release(i);
                         });

    Kokkos::View<int, Kokkos::OpenMP> sum("");
    Kokkos::parallel_for(
        Kokkos::RangePolicy<Kokkos::OpenMP>(0, token.size()),
        [=](const int i) { Kokkos::atomic_add(sum.data(), count[i]); });

    if (sum() != 1000) {
      Kokkos::atomic_add(&errors, 1);
    }
  };

  master(0, 1);

  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 4, 0);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 0, 4);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 2, 2);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 8, 0);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 0, 8);
  ASSERT_EQ(errors, 0);

  Kokkos::OpenMP::partition_master(master, 8, 8);
  ASSERT_EQ(errors, 0);
}

}  // namespace Test
