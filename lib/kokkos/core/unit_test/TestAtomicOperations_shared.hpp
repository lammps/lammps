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

namespace Test {

// FIXME_SYCL This doesn't work yet for SYCL+CUDA
#if !defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ARCH_INTEL_GPU)
template <typename ExecutionSpace>
struct TestSharedAtomicsFunctor {
  Kokkos::View<int, typename ExecutionSpace::memory_space> m_view;

  TestSharedAtomicsFunctor(
      Kokkos::View<int, typename ExecutionSpace::memory_space>& view)
      : m_view(view) {}

  KOKKOS_INLINE_FUNCTION void operator()(
      const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type t) const {
    int* x = (int*)t.team_shmem().get_shmem(sizeof(int));
    Kokkos::single(Kokkos::PerTeam(t), [&]() { *x = 0; });
    t.team_barrier();
    Kokkos::atomic_add(x, 1);
    t.team_barrier();
    Kokkos::single(Kokkos::PerTeam(t), [&]() { m_view() = *x; });
  }
};

TEST(TEST_CATEGORY, atomic_shared) {
  TEST_EXECSPACE exec;
  Kokkos::View<int, typename TEST_EXECSPACE::memory_space> view("ref_value");
  auto team_size =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, Kokkos::AUTO)
          .team_size_recommended(TestSharedAtomicsFunctor<TEST_EXECSPACE>(view),
                                 Kokkos::ParallelForTag{});
  Kokkos::parallel_for(Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, team_size)
                           .set_scratch_size(0, Kokkos::PerTeam(8)),
                       TestSharedAtomicsFunctor<TEST_EXECSPACE>(view));
  exec.fence("Fence after test kernel");
  int i = 0;
  Kokkos::deep_copy(i, view);
  ASSERT_EQ(i, team_size);
}
#endif
}  // namespace Test
