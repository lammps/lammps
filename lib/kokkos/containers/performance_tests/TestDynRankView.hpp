
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

#ifndef KOKKOS_TEST_DYNRANKVIEW_HPP
#define KOKKOS_TEST_DYNRANKVIEW_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>
#include <vector>

#include <impl/Kokkos_Timer.hpp>

// Compare performance of DynRankView to View, specific focus on the parenthesis operators

namespace Performance {

//View functor
template <typename DeviceType>
struct InitViewFunctor {
  typedef Kokkos::View<double***, DeviceType> inviewtype;
  inviewtype _inview;

  InitViewFunctor( inviewtype &inview_ ) : _inview(inview_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (unsigned j = 0; j < _inview.extent(1); ++j) {
      for (unsigned k = 0; k < _inview.extent(2); ++k) {
        _inview(i,j,k) = i/2 -j*j + k/3;
      }
    }
  }

  struct SumComputationTest
  {
    typedef Kokkos::View<double***, DeviceType> inviewtype;
    inviewtype _inview;

    typedef Kokkos::View<double*, DeviceType> outviewtype;
    outviewtype _outview;

    KOKKOS_INLINE_FUNCTION
    SumComputationTest(inviewtype &inview_ , outviewtype &outview_) : _inview(inview_), _outview(outview_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      for (unsigned j = 0; j < _inview.extent(1); ++j) {
        for (unsigned k = 0; k < _inview.extent(2); ++k) {
          _outview(i) += _inview(i,j,k) ;
        }
      }
    }
  };

};

template <typename DeviceType>
struct InitStrideViewFunctor {
  typedef Kokkos::View<double***, Kokkos::LayoutStride, DeviceType> inviewtype;
  inviewtype _inview;

  InitStrideViewFunctor( inviewtype &inview_ ) : _inview(inview_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (unsigned j = 0; j < _inview.extent(1); ++j) {
      for (unsigned k = 0; k < _inview.extent(2); ++k) {
        _inview(i,j,k) = i/2 -j*j + k/3;
      }
    }
  }

};

template <typename DeviceType>
struct InitViewRank7Functor {
  typedef Kokkos::View<double*******, DeviceType> inviewtype;
  inviewtype _inview;

  InitViewRank7Functor( inviewtype &inview_ ) : _inview(inview_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (unsigned j = 0; j < _inview.extent(1); ++j) {
      for (unsigned k = 0; k < _inview.extent(2); ++k) {
        _inview(i,j,k,0,0,0,0) = i/2 -j*j + k/3;
      }
    }
  }

};

//DynRankView functor
template <typename DeviceType>
struct InitDynRankViewFunctor {
  typedef Kokkos::DynRankView<double, DeviceType> inviewtype;
  inviewtype _inview;

  InitDynRankViewFunctor( inviewtype &inview_ ) : _inview(inview_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (unsigned j = 0; j < _inview.extent(1); ++j) {
      for (unsigned k = 0; k < _inview.extent(2); ++k) {
        _inview(i,j,k) = i/2 -j*j + k/3;
      }
    }
  }

  struct SumComputationTest
  {
    typedef Kokkos::DynRankView<double, DeviceType> inviewtype;
    inviewtype _inview;

    typedef Kokkos::DynRankView<double, DeviceType> outviewtype;
    outviewtype _outview;

    KOKKOS_INLINE_FUNCTION
    SumComputationTest(inviewtype &inview_ , outviewtype &outview_) : _inview(inview_), _outview(outview_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      for (unsigned j = 0; j < _inview.extent(1); ++j) {
        for (unsigned k = 0; k < _inview.extent(2); ++k) {
          _outview(i) += _inview(i,j,k) ;
        }
      }
    }
  };

};


template <typename DeviceType>
void test_dynrankview_op_perf( const int par_size )
{

  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;
  const size_type dim_2 = 90;
  const size_type dim_3 = 30;

  double elapsed_time_view = 0;
  double elapsed_time_compview = 0;
  double elapsed_time_strideview = 0;
  double elapsed_time_view_rank7 = 0;
  double elapsed_time_drview = 0;
  double elapsed_time_compdrview = 0;
  Kokkos::Timer timer;
  {
    Kokkos::View<double***,DeviceType> testview("testview",par_size,dim_2,dim_3);
    typedef InitViewFunctor<DeviceType> FunctorType;

    timer.reset();
    Kokkos::RangePolicy<DeviceType> policy(0,par_size);
    Kokkos::parallel_for( policy , FunctorType(testview) );
    DeviceType::fence();
    elapsed_time_view = timer.seconds();
    std::cout << " View time (init only): " << elapsed_time_view << std::endl;


    timer.reset();
    Kokkos::View<double*,DeviceType> sumview("sumview",par_size);
    Kokkos::parallel_for( policy , typename FunctorType::SumComputationTest(testview, sumview) );
    DeviceType::fence();
    elapsed_time_compview = timer.seconds();
    std::cout << " View sum computation time: " << elapsed_time_view << std::endl;


    Kokkos::View<double***,Kokkos::LayoutStride, DeviceType> teststrideview = Kokkos::subview(testview, Kokkos::ALL, Kokkos::ALL,Kokkos::ALL);
    typedef InitStrideViewFunctor<DeviceType> FunctorStrideType;

    timer.reset();
    Kokkos::parallel_for( policy , FunctorStrideType(teststrideview) );
    DeviceType::fence();
    elapsed_time_strideview = timer.seconds();
    std::cout << " Strided View time (init only): " << elapsed_time_strideview << std::endl;
  }
  {
    Kokkos::View<double*******,DeviceType> testview("testview",par_size,dim_2,dim_3,1,1,1,1);
    typedef InitViewRank7Functor<DeviceType> FunctorType;

    timer.reset();
    Kokkos::RangePolicy<DeviceType> policy(0,par_size);
    Kokkos::parallel_for( policy , FunctorType(testview) );
    DeviceType::fence();
    elapsed_time_view_rank7 = timer.seconds();
    std::cout << " View Rank7 time (init only): " << elapsed_time_view_rank7 << std::endl;
  }
  {
    Kokkos::DynRankView<double,DeviceType> testdrview("testdrview",par_size,dim_2,dim_3);
    typedef InitDynRankViewFunctor<DeviceType> FunctorType;

    timer.reset();
    Kokkos::RangePolicy<DeviceType> policy(0,par_size);
    Kokkos::parallel_for( policy , FunctorType(testdrview) );
    DeviceType::fence();
    elapsed_time_drview = timer.seconds();
    std::cout << " DynRankView time (init only): " << elapsed_time_drview << std::endl;

    timer.reset();
    Kokkos::DynRankView<double,DeviceType> sumview("sumview",par_size);
    Kokkos::parallel_for( policy , typename FunctorType::SumComputationTest(testdrview, sumview) );
    DeviceType::fence();
    elapsed_time_compdrview = timer.seconds();
    std::cout << " DynRankView sum computation time: " << elapsed_time_compdrview << std::endl;

  }

  std::cout << " Ratio of View to DynRankView time: " << elapsed_time_view / elapsed_time_drview << std::endl; //expect < 1
  std::cout << " Ratio of View to DynRankView sum computation time: " << elapsed_time_compview / elapsed_time_compdrview << std::endl; //expect < 1
  std::cout << " Ratio of View to View Rank7  time: " << elapsed_time_view / elapsed_time_view_rank7 << std::endl; //expect < 1
  std::cout << " Ratio of StrideView to DynRankView time: " << elapsed_time_strideview / elapsed_time_drview << std::endl; //expect < 1
  std::cout << " Ratio of DynRankView to View Rank7  time: " << elapsed_time_drview / elapsed_time_view_rank7 << std::endl; //expect ?

  timer.reset();

} //end test_dynrankview


} //end Performance
#endif

