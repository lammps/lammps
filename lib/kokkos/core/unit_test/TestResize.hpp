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
#ifndef TESTVIEWSUBVIEW_HPP_
#define TESTVIEWSUBVIEW_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace TestViewResize {

template<class DeviceType>
void testResize ()
{
  const int sizes[8] = {2, 3, 4, 5, 6, 7, 8, 9};

  // Check #904 fix (no reallocation if dimensions didn't change).
  {
    typedef Kokkos::View<int*, DeviceType> view_type;
    view_type view_1d ("view_1d", sizes[0]);
    const int* oldPointer = view_1d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_1d, sizes[0]);
    const int* newPointer = view_1d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int**, DeviceType> view_type;
    view_type view_2d ("view_2d", sizes[0], sizes[1]);
    const int* oldPointer = view_2d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_2d, sizes[0], sizes[1]);
    const int* newPointer = view_2d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int***, DeviceType> view_type;
    view_type view_3d ("view_3d", sizes[0], sizes[1], sizes[2]);
    const int* oldPointer = view_3d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_3d, sizes[0], sizes[1], sizes[2]);
    const int* newPointer = view_3d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int****, DeviceType> view_type;
    view_type view_4d ("view_4d", sizes[0], sizes[1], sizes[2], sizes[3]);
    const int* oldPointer = view_4d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_4d, sizes[0], sizes[1], sizes[2], sizes[3]);
    const int* newPointer = view_4d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int*****, DeviceType> view_type;
    view_type view_5d ("view_5d", sizes[0], sizes[1], sizes[2], sizes[3],
                       sizes[4]);
    const int* oldPointer = view_5d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_5d, sizes[0], sizes[1], sizes[2], sizes[3], sizes[4]);
    const int* newPointer = view_5d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int******, DeviceType> view_type;
    view_type view_6d ("view_6d", sizes[0], sizes[1], sizes[2], sizes[3],
                       sizes[4], sizes[5]);
    const int* oldPointer = view_6d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_6d, sizes[0], sizes[1], sizes[2], sizes[3], sizes[4],
                    sizes[5]);
    const int* newPointer = view_6d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int*******, DeviceType> view_type;
    view_type view_7d ("view_7d", sizes[0], sizes[1], sizes[2], sizes[3],
                       sizes[4], sizes[5], sizes[6]);
    const int* oldPointer = view_7d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_7d, sizes[0], sizes[1], sizes[2], sizes[3], sizes[4],
                    sizes[5], sizes[6]);
    const int* newPointer = view_7d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
  {
    typedef Kokkos::View<int********, DeviceType> view_type;
    view_type view_8d ("view_8d", sizes[0], sizes[1], sizes[2], sizes[3],
                       sizes[4], sizes[5], sizes[6], sizes[7]);
    const int* oldPointer = view_8d.data ();
    EXPECT_TRUE( oldPointer != NULL );
    Kokkos::resize (view_8d, sizes[0], sizes[1], sizes[2], sizes[3], sizes[4],
                    sizes[5], sizes[6], sizes[7]);
    const int* newPointer = view_8d.data ();
    EXPECT_TRUE( oldPointer == newPointer );
  }
}

} // namespace TestViewSubview

#endif // TESTVIEWSUBVIEW_HPP_
