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

/*
 * FIXME the OffsetView class is really not very well tested.
 */
#ifndef CONTAINERS_UNIT_TESTS_TESTOFFSETVIEW_HPP_
#define CONTAINERS_UNIT_TESTS_TESTOFFSETVIEW_HPP_



#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_OffsetView.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

using std::endl;
using std::cout;

namespace Test{

   template <typename Scalar, typename Device>
   void test_offsetview_construction(unsigned int size)
   {

      typedef Kokkos::Experimental::OffsetView<Scalar**, Device> offset_view_type;
      typedef Kokkos::View<Scalar**, Device> view_type;

      Kokkos::Experimental::index_list_type range0 = {-1, 3};
      Kokkos::Experimental::index_list_type range1 = {-2, 2};

      offset_view_type ov("firstOV", range0, range1);

      ASSERT_EQ("firstOV", ov.label());
      ASSERT_EQ(2, ov.Rank);

      ASSERT_EQ(ov.begin(0), -1);
      ASSERT_EQ(ov.end(0), 4);

      ASSERT_EQ(ov.begin(1), -2);
      ASSERT_EQ(ov.end(1), 3);

      ASSERT_EQ(ov.extent(0), 5);
      ASSERT_EQ(ov.extent(1), 5);

      const int ovmin0 = ov.begin(0);
      const int ovend0 = ov.end(0);
      const int ovmin1 = ov.begin(1);
      const int ovend1 = ov.end(1);

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
      {
         Kokkos::Experimental::OffsetView<Scalar*, Device> offsetV1("OneDOffsetView", range0);

         Kokkos::RangePolicy<Device, int> rangePolicy1(offsetV1.begin(0), offsetV1.end(0));
         Kokkos::parallel_for(rangePolicy1, KOKKOS_LAMBDA (const int i){
            offsetV1(i) = 1;
         }
         );
	 Kokkos::fence();

         int OVResult = 0;
         Kokkos::parallel_reduce(rangePolicy1, KOKKOS_LAMBDA(const int i, int & updateMe){
            updateMe += offsetV1(i);
         }, OVResult);
	 
	 Kokkos::fence();
         ASSERT_EQ(OVResult, offsetV1.end(0) - offsetV1.begin(0)) << "found wrong number of elements in OffsetView that was summed.";

      }
      {  //test deep copy of scalar const value into mirro
         const int constVal = 6;
         typename offset_view_type::HostMirror hostOffsetView =
               Kokkos::Experimental::create_mirror_view(ov);

         Kokkos::Experimental::deep_copy(hostOffsetView, constVal);

         for(int i = hostOffsetView.begin(0); i < hostOffsetView.end(0); ++i) {
            for(int j = hostOffsetView.begin(1); j < hostOffsetView.end(1); ++j) {
               ASSERT_EQ(hostOffsetView(i,j),  constVal) << "Bad data found in OffsetView";
            }
         }
      }

      typedef Kokkos::MDRangePolicy<Device, Kokkos::Rank<2>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::point_type point_type;

      range_type rangePolicy2D(point_type{ {ovmin0, ovmin1 } },
            point_type{ { ovend0, ovend1 } });

      const int constValue = 9;
      Kokkos::parallel_for(rangePolicy2D, KOKKOS_LAMBDA (const int i, const int j) {
         ov(i,j) =  constValue;
      }
      );
      
      //test offsetview to offsetviewmirror deep copy
      typename offset_view_type::HostMirror hostOffsetView =
            Kokkos::Experimental::create_mirror_view(ov);

      Kokkos::Experimental::deep_copy(hostOffsetView, ov);

      for(int i = hostOffsetView.begin(0); i < hostOffsetView.end(0); ++i) {
         for(int j = hostOffsetView.begin(1); j < hostOffsetView.end(1); ++j) {
            ASSERT_EQ(hostOffsetView(i,j),  constValue) << "Bad data found in OffsetView";
         }
      }
      
     int OVResult = 0;
      Kokkos::parallel_reduce(rangePolicy2D, KOKKOS_LAMBDA(const int i, const int j, int & updateMe){
         updateMe += ov(i, j);
      }, OVResult);

      int answer = 0;
      for(int i = ov.begin(0); i < ov.end(0); ++i) {
         for(int j = ov.begin(1); j < ov.end(1); ++j) {
            answer += constValue;
         }
      }
      
      ASSERT_EQ(OVResult, answer) << "Bad data found in OffsetView";
#endif

      {
         offset_view_type ovCopy(ov);
         ASSERT_EQ(ovCopy==ov, true) <<
               "Copy constructor or equivalence operator broken";
      }
      
      {
         offset_view_type ovAssigned = ov;
         ASSERT_EQ(ovAssigned==ov, true) <<
               "Assignment operator or equivalence operator broken";
      }
      
      {  //construct OffsetView from a View plus begins array
         const int extent0 = 100;
         const int extent1 = 200;
         const int extent2 = 300;
         Kokkos::View<Scalar***, Device> view3D("view3D", extent0, extent1, extent2);

         Kokkos::deep_copy(view3D, 1);

         Kokkos::Array<int64_t,3> begins = {{-10, -20, -30}};
         Kokkos::Experimental::OffsetView<Scalar***, Device> offsetView3D(view3D, begins);

         typedef Kokkos::MDRangePolicy<Device, Kokkos::Rank<3>, Kokkos::IndexType<int64_t> > range3_type;
         typedef typename range3_type::point_type point3_type;

         range3_type rangePolicy3DZero(point3_type{ {0, 0, 0 } },
               point3_type{ { extent0, extent1, extent2 } });

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
        int view3DSum = 0;
         Kokkos::parallel_reduce(rangePolicy3DZero, KOKKOS_LAMBDA(const int i, const int j, int k, int & updateMe){
            updateMe += view3D(i, j, k);
         }, view3DSum);

         range3_type rangePolicy3D(point3_type{ {begins[0], begins[1], begins[2] } },
               point3_type{ { begins[0] + extent0, begins[1] + extent1, begins[2] + extent2 } });
         int offsetView3DSum = 0;

         Kokkos::parallel_reduce(rangePolicy3D, KOKKOS_LAMBDA(const int i, const int j, int k, int & updateMe){
            updateMe += offsetView3D(i, j, k);
         }, offsetView3DSum);

         ASSERT_EQ(view3DSum, offsetView3DSum) << "construction of OffsetView from View and begins array broken.";
#endif
      }
      view_type viewFromOV = ov.view();

      ASSERT_EQ(viewFromOV == ov, true) <<
            "OffsetView::view() or equivalence operator View == OffsetView broken";

      {
         offset_view_type ovFromV(viewFromOV, {-1, -2});

         ASSERT_EQ(ovFromV == viewFromOV , true) <<
               "Construction of OffsetView from View or equivalence operator OffsetView == View broken";
      }
      {
         offset_view_type ovFromV = viewFromOV;
         ASSERT_EQ(ovFromV == viewFromOV , true) <<
               "Construction of OffsetView from View by assignment (implicit conversion) or equivalence operator OffsetView == View broken";
      }

      {// test offsetview to view deep copy
         view_type aView("aView", ov.extent(0), ov.extent(1));
         Kokkos::Experimental::deep_copy(aView, ov);

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
         int sum = 0;
         Kokkos::parallel_reduce(rangePolicy2D, KOKKOS_LAMBDA(const int i, const int j, int & updateMe){
            updateMe += ov(i, j) - aView(i- ov.begin(0), j-ov.begin(1));
         }, sum);

         ASSERT_EQ(sum, 0) << "deep_copy(view, offsetView) broken.";
#endif
      }

      {// test view to  offsetview deep copy
         view_type aView("aView", ov.extent(0), ov.extent(1));

         Kokkos::deep_copy(aView, 99);
         Kokkos::Experimental::deep_copy(ov, aView);
	 

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
         int sum = 0;
         Kokkos::parallel_reduce(rangePolicy2D, KOKKOS_LAMBDA(const int i, const int j, int & updateMe){
            updateMe += ov(i, j) - aView(i- ov.begin(0), j-ov.begin(1));
         }, sum);

         ASSERT_EQ(sum, 0) << "deep_copy(offsetView, view) broken.";
#endif
      }
   }
   template <typename Scalar, typename Device>
   void test_offsetview_subview(unsigned int size)
   {
      {//test subview 1
          Kokkos::Experimental::OffsetView<Scalar*, Device> sliceMe("offsetToSlice", {-10, 20});
          {
             auto offsetSubviewa = Kokkos::Experimental::subview(sliceMe, 0);
             ASSERT_EQ(offsetSubviewa.Rank, 0) << "subview of offset is broken.";
          }

       }
      {//test subview 2
         Kokkos::Experimental::OffsetView<Scalar**, Device> sliceMe("offsetToSlice", {-10,20}, {-20,30});
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(),-2);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }

         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
      }


      {//test subview rank 3

         Kokkos::Experimental::OffsetView<Scalar***, Device> sliceMe("offsetToSlice", {-10,20}, {-20,30}, {-30,40});

         //slice 1
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe,Kokkos::ALL(),Kokkos::ALL(), 0);
            ASSERT_EQ(offsetSubview.Rank, 2) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe,Kokkos::ALL(), 0,Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 2) << "subview of offset is broken.";
         }

         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe,0, Kokkos::ALL(),Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 2) << "subview of offset is broken.";

         }
         {
	   auto offsetSubview = Kokkos::Experimental::subview(sliceMe,0, Kokkos::ALL(), Kokkos::make_pair(-30, -21));
            ASSERT_EQ(offsetSubview.Rank, 2) << "subview of offset is broken.";

            ASSERT_EQ(offsetSubview.begin(0) , -20);
            ASSERT_EQ(offsetSubview.end(0) , 31);
            ASSERT_EQ(offsetSubview.begin(1) , 0);
            ASSERT_EQ(offsetSubview.end(1) , 9);

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
            typedef Kokkos::MDRangePolicy<Device, Kokkos::Rank<2>, Kokkos::IndexType<int> > range_type;
            typedef typename range_type::point_type point_type;

            const int b0 = offsetSubview.begin(0);
            const int b1 = offsetSubview.begin(1);

            const int e0 = offsetSubview.end(0);
            const int e1 = offsetSubview.end(1);

            range_type rangeP2D(point_type{ {b0, b1 } }, point_type{ { e0, e1} });

            Kokkos::parallel_for(rangeP2D, KOKKOS_LAMBDA(const int i, const int j) {
               offsetSubview(i,j) =  6;
            }
            );

            int sum = 0;
             Kokkos::parallel_reduce(rangeP2D, KOKKOS_LAMBDA(const int i, const int j, int & updateMe){
                updateMe += offsetSubview(i, j);
             }, sum);

            ASSERT_EQ(sum, 6*(e0-b0)*(e1-b1));
#endif
         }

         // slice 2
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), 0, 0);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, 0, 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }

         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, 0, Kokkos::ALL(), 0);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
      }

      {//test subview rank 4

         Kokkos::Experimental::OffsetView<Scalar****, Device> sliceMe("offsetToSlice", {-10,20}, {-20,30}, {-30,40}, {-40, 50});

         //slice 1
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(),Kokkos::ALL(), Kokkos::ALL(), 0);
            ASSERT_EQ(offsetSubview.Rank, 3) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), Kokkos::ALL(), 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 3) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe ,Kokkos::ALL(), 0, Kokkos::ALL(),Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 3) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe , 0, Kokkos::ALL(), Kokkos::ALL(),  Kokkos::ALL() );
            ASSERT_EQ(offsetSubview.Rank, 3) << "subview of offset is broken.";
         }

         // slice 2
         auto offsetSubview2a = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), Kokkos::ALL(), 0, 0);
         ASSERT_EQ(offsetSubview2a.Rank, 2) << "subview of offset is broken.";
         {
            auto offsetSubview2b = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), 0, Kokkos::ALL(), 0);
            ASSERT_EQ(offsetSubview2b.Rank, 2) << "subview of offset is broken.";
         }
         {
            auto offsetSubview2b = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), 0, 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview2b.Rank, 2) << "subview of offset is broken.";
         }
         {
            auto offsetSubview2b = Kokkos::Experimental::subview(sliceMe,  0, Kokkos::ALL(), 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview2b.Rank, 2) << "subview of offset is broken.";
         }
         {
            auto offsetSubview2b = Kokkos::Experimental::subview(sliceMe,  0, 0, Kokkos::ALL(), Kokkos::ALL());
            ASSERT_EQ(offsetSubview2b.Rank, 2) << "subview of offset is broken.";
         }
         // slice 3
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, Kokkos::ALL(), 0, 0, 0);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe, 0, Kokkos::ALL(), 0, 0);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe,  0, 0, Kokkos::ALL(), 0);
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }
         {
            auto offsetSubview = Kokkos::Experimental::subview(sliceMe,  0, 0, 0, Kokkos::ALL());
            ASSERT_EQ(offsetSubview.Rank, 1) << "subview of offset is broken.";
         }

      }

   }

   TEST_F( TEST_CATEGORY, offsetview_construction) {
      test_offsetview_construction<int,TEST_EXECSPACE>(10);
   }
   TEST_F( TEST_CATEGORY, offsetview_subview) {
      test_offsetview_subview<int,TEST_EXECSPACE>(10);
   }

} // namespace Test

#endif /* CONTAINERS_UNIT_TESTS_TESTOFFSETVIEW_HPP_ */
