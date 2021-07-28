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

#include <Kokkos_Core.hpp>
#include <TestHPX_Category.hpp>

#ifdef KOKKOS_ENABLE_HPX_ASYNC_DISPATCH

namespace Test {

TEST(hpx, instance_ids) {
  Kokkos::InitArguments arguments{-1, -1, -1, false};
  Kokkos::initialize(arguments);

  {
    Kokkos::Experimental::HPX hpx_global1;
    Kokkos::Experimental::HPX hpx_global2 = hpx_global1;
    Kokkos::Experimental::HPX hpx_global3{hpx_global1};
    Kokkos::Experimental::HPX hpx_global4(
        Kokkos::Experimental::HPX::instance_mode::global);

    ASSERT_EQ(0, hpx_global1.impl_instance_id());
    ASSERT_EQ(0, hpx_global2.impl_instance_id());
    ASSERT_EQ(0, hpx_global3.impl_instance_id());
    ASSERT_EQ(0, hpx_global4.impl_instance_id());

    Kokkos::Experimental::HPX hpx_independent1(
        Kokkos::Experimental::HPX::instance_mode::independent);
    Kokkos::Experimental::HPX hpx_independent2 = hpx_independent1;
    Kokkos::Experimental::HPX hpx_independent3{hpx_independent1};

    ASSERT_NE(hpx_global1.impl_instance_id(),
              hpx_independent1.impl_instance_id());
    ASSERT_EQ(hpx_independent1.impl_instance_id(),
              hpx_independent2.impl_instance_id());
    ASSERT_EQ(hpx_independent1.impl_instance_id(),
              hpx_independent3.impl_instance_id());

    hpx::shared_future<void> f = hpx::make_ready_future<void>();
    Kokkos::Experimental::HPX hpx_independent_future1(f);
    Kokkos::Experimental::HPX hpx_independent_future2 = hpx_independent_future1;
    Kokkos::Experimental::HPX hpx_independent_future3{hpx_independent_future1};

    ASSERT_NE(hpx_global1.impl_instance_id(),
              hpx_independent1.impl_instance_id());
    ASSERT_NE(hpx_independent1.impl_instance_id(),
              hpx_independent_future1.impl_instance_id());
    ASSERT_EQ(hpx_independent_future1.impl_instance_id(),
              hpx_independent_future2.impl_instance_id());
    ASSERT_EQ(hpx_independent_future1.impl_instance_id(),
              hpx_independent_future3.impl_instance_id());
  }

  Kokkos::finalize();
}
}  // namespace Test

#endif
