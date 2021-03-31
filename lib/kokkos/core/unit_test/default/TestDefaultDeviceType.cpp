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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <default/TestDefaultDeviceType_Category.hpp>
#include <TestHalfConversion.hpp>
#include <TestHalfOperators.hpp>

#if !defined(KOKKOS_ENABLE_CUDA) || defined(__CUDACC__)

namespace Test {

TEST(TEST_CATEGORY, host_space_access) {
  using host_exec_space = Kokkos::HostSpace::execution_space;
  using device_space    = Kokkos::Device<host_exec_space, Kokkos::HostSpace>;
  using mirror_space =
      Kokkos::Impl::HostMirror<Kokkos::DefaultExecutionSpace>::Space;

  static_assert(Kokkos::Impl::SpaceAccessibility<host_exec_space,
                                                 Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::Impl::SpaceAccessibility<device_space,
                                                 Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::Impl::SpaceAccessibility<mirror_space,
                                                 Kokkos::HostSpace>::accessible,
                "");
}

}  // namespace Test

#endif
