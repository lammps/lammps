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

#include <hpx/local/future.hpp>

#ifdef KOKKOS_ENABLE_HPX_ASYNC_DISPATCH

namespace Test {

TEST(hpx, delayed_execution) {
  Kokkos::InitArguments arguments{-1, -1, -1, false};
  Kokkos::initialize(arguments);

  {
    Kokkos::View<bool, Kokkos::Experimental::HPX> ran("ran");
    hpx::lcos::local::promise<void> p;
    hpx::shared_future<void> f = p.get_future();

    Kokkos::Experimental::HPX hpx(f);
    Kokkos::parallel_for(
        "Test::hpx::independent_instances::delay_execution",
        Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx, 0, 1),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        KOKKOS_LAMBDA(int) { ran() = true; });

    ASSERT_EQ(false, ran());
    ASSERT_EQ(false, hpx.impl_get_future().is_ready());

    p.set_value();

    hpx.fence();
    ASSERT_EQ(true, hpx.impl_get_future().is_ready());
  }

  Kokkos::finalize();
}
}  // namespace Test

#endif
