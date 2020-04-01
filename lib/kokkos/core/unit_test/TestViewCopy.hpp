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

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template <typename ExecSpace>
struct TestViewCopy {
  using InExecSpace = ExecSpace;

  static void test_view_copy(const int dim0, const int dim1, const int dim2) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_ROCM)
    // ExecSpace = CudaUVM, CudaHostPinned
    // This test will fail at runtime with an illegal memory access if something
    // goes wrong Test 1: deep_copy from host_mirror_space to ExecSpace and
    // ExecSpace back to host_mirror_space
    {
      typedef Kokkos::View<double****, InExecSpace> Rank4ViewType;
      Rank4ViewType view_4;
      view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

      typedef typename Kokkos::Impl::is_space<
          InExecSpace>::host_mirror_space::execution_space host_space_type;
      Kokkos::View<double**, Kokkos::LayoutLeft, host_space_type> srcView(
          "srcView", dim2, dim2);

      // Strided dst view
      auto dstView =
          Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

      // host_mirror_space to ExecSpace
      Kokkos::deep_copy(dstView, srcView);
      Kokkos::fence();

      // ExecSpace to host_mirror_space
      Kokkos::deep_copy(srcView, dstView);
      Kokkos::fence();
    }

    // Test 2: deep_copy from Cuda to ExecSpace and ExecSpace back to Cuda
    {
      typedef Kokkos::View<double****, InExecSpace> Rank4ViewType;
      Rank4ViewType view_4;
      view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

#if defined(KOKKOS_ENABLE_CUDA)
      typedef typename std::conditional<
          Kokkos::Impl::MemorySpaceAccess<
              Kokkos::CudaSpace,
              typename InExecSpace::memory_space>::accessible,
          Kokkos::CudaSpace, InExecSpace>::type space_type;
#endif
#if defined(KOKKOS_ENABLE_ROCM)
      typedef typename std::conditional<
          Kokkos::Impl::MemorySpaceAccess<
              Kokkos::ROCmSpace,
              typename InExecSpace::memory_space>::accessible,
          Kokkos::ROCmSpace, InExecSpace>::type space_type;
#endif
      Kokkos::View<double**, Kokkos::LayoutLeft, space_type> srcView(
          "srcView", dim2, dim2);

      // Strided dst view
      auto dstView =
          Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

      // Cuda to ExecSpace
      Kokkos::deep_copy(dstView, srcView);
      Kokkos::fence();

      // ExecSpace to Cuda
      Kokkos::deep_copy(srcView, dstView);
      Kokkos::fence();
    }

    // Test 3: deep_copy from host_space to ExecSpace and ExecSpace back to
    // host_space
    {
      typedef Kokkos::View<double****, InExecSpace> Rank4ViewType;
      Rank4ViewType view_4;
      view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

      typedef Kokkos::HostSpace host_space_type;
      Kokkos::View<double**, Kokkos::LayoutLeft, host_space_type> srcView(
          "srcView", dim2, dim2);

      // Strided dst view
      auto dstView =
          Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

      // host_space to ExecSpace
      Kokkos::deep_copy(dstView, srcView);
      Kokkos::fence();

      // ExecSpace to host_space
      Kokkos::deep_copy(srcView, dstView);
      Kokkos::fence();
    }
#endif
  }  // end test_view_copy

};  // end struct

}  // namespace

TEST(TEST_CATEGORY, view_copy_tests) {
  // Only include this file to be compiled with CudaUVM and CudaHostPinned
  TestViewCopy<TEST_EXECSPACE>::test_view_copy(4, 2, 3);
  TestViewCopy<TEST_EXECSPACE>::test_view_copy(4, 2, 0);
}

TEST(TEST_CATEGORY, view_copy_degenerated) {
  // Only include this file to be compiled with CudaUVM and CudaHostPinned
  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v_um_def_1;
  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v_um_1(
      reinterpret_cast<int*>(-1), 0);
  Kokkos::View<int*> v_m_def_1;
  Kokkos::View<int*> v_m_1("v_m_1", 0);

  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v_um_def_2;
  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v_um_2(
      reinterpret_cast<int*>(-1), 0);
  Kokkos::View<int*> v_m_def_2;
  Kokkos::View<int*> v_m_2("v_m_2", 0);

  Kokkos::deep_copy(v_um_def_1, v_um_def_2);
  Kokkos::deep_copy(v_um_def_1, v_um_2);
  Kokkos::deep_copy(v_um_def_1, v_m_def_2);
  Kokkos::deep_copy(v_um_def_1, v_m_2);

  Kokkos::deep_copy(v_um_1, v_um_def_2);
  Kokkos::deep_copy(v_um_1, v_um_2);
  Kokkos::deep_copy(v_um_1, v_m_def_2);
  Kokkos::deep_copy(v_um_1, v_m_2);

  Kokkos::deep_copy(v_m_def_1, v_um_def_2);
  Kokkos::deep_copy(v_m_def_1, v_um_2);
  Kokkos::deep_copy(v_m_def_1, v_m_def_2);
  Kokkos::deep_copy(v_m_def_1, v_m_2);

  Kokkos::deep_copy(v_m_1, v_um_def_2);
  Kokkos::deep_copy(v_m_1, v_um_2);
  Kokkos::deep_copy(v_m_1, v_m_def_2);
  Kokkos::deep_copy(v_m_1, v_m_2);
}

}  // namespace Test
