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

#ifndef KOKKOS_TEST_SCATTER_VIEW_HPP
#define KOKKOS_TEST_SCATTER_VIEW_HPP

#include <Kokkos_ScatterView.hpp>

namespace Test {

template <typename ExecSpace, typename Layout, int duplication, int contribution, int op>
struct test_scatter_view_impl_cls;

template <typename ExecSpace, typename Layout, int duplication, int contribution>
struct test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, Kokkos::Experimental::ScatterSum>   
{
public:   

   typedef Kokkos::Experimental::ScatterView
       < double*[3]
       , Layout
       , ExecSpace
       , Kokkos::Experimental::ScatterSum
       , duplication
       , contribution
        > scatter_view_type;

   typedef Kokkos::View<double *[3], Layout, ExecSpace> orig_view_type; 


   scatter_view_type scatter_view;
   int scatterSize;

   test_scatter_view_impl_cls(const scatter_view_type& view){
      scatter_view = view;
      scatterSize = 0;
   }

   void initialize(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        host_view(i, 0) = 0.0;
        host_view(i, 1) = 0.0;
        host_view(i, 2) = 0.0;
      }
      Kokkos::fence();
      Kokkos::deep_copy(orig, host_view);
   }

   void run_parallel(int n) {
        scatterSize = n;
        auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
        Kokkos::parallel_for(policy, *this, "scatter_view_test: Sum");
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(int i) const {
      auto scatter_access = scatter_view.access();
      auto scatter_access_atomic = scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
      for (int j = 0; j < 10; ++j) {
        auto k = (i + j) % scatterSize;
        scatter_access(k, 0) += 4.2;
        scatter_access_atomic(k, 1) += 2.0;
        scatter_access(k, 2) += 1.0;
      }
    }

    void validateResults(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        auto val0 = host_view(i, 0);
        auto val1 = host_view(i, 1);
        auto val2 = host_view(i, 2);
        EXPECT_TRUE(std::fabs((val0 - 84.0) / 84.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val1 - 40.0) / 40.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val2 - 20.0) / 20.0) < 1e-14);
      }
    }
};


template <typename ExecSpace, typename Layout, int duplication, int contribution>
struct test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, Kokkos::Experimental::ScatterProd>   
{
public:   

   typedef Kokkos::Experimental::ScatterView
       < double*[3]
       , Layout
       , ExecSpace
       , Kokkos::Experimental::ScatterProd
       , duplication
       , contribution
        > scatter_view_type;

   typedef Kokkos::View<double *[3], Layout, ExecSpace> orig_view_type; 


   scatter_view_type scatter_view;
   int scatterSize;

   test_scatter_view_impl_cls(const scatter_view_type& view){
      scatter_view = view;
      scatterSize = 0;
   }

   void initialize(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        host_view(i, 0) = 1.0;
        host_view(i, 1) = 1.0;
        host_view(i, 2) = 1.0;
      }
      Kokkos::fence();
      Kokkos::deep_copy(orig, host_view);
   }

   void run_parallel(int n) {
        scatterSize = n;
        auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
        Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(int i) const {
      auto scatter_access = scatter_view.access();
      auto scatter_access_atomic = scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
      for (int j = 0; j < 4; ++j) {
        auto k = (i + j) % scatterSize;
        scatter_access(k, 0) *= 4.0;
        scatter_access_atomic(k, 1) *= 2.0;
        scatter_access(k, 2) *= 1.0;
      }
    }

    void validateResults(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        auto val0 = host_view(i, 0);
        auto val1 = host_view(i, 1);
        auto val2 = host_view(i, 2);
        EXPECT_TRUE(std::fabs((val0 - 65536.0) / 65536.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val1 - 256.0) / 256.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14);
      }
    }
};


template <typename ExecSpace, typename Layout, int duplication, int contribution>
struct test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, Kokkos::Experimental::ScatterMin>   
{
public:   

   typedef Kokkos::Experimental::ScatterView
       < double*[3]
       , Layout
       , ExecSpace
       , Kokkos::Experimental::ScatterMin
       , duplication
       , contribution
        > scatter_view_type;

   typedef Kokkos::View<double *[3], Layout, ExecSpace> orig_view_type; 


   scatter_view_type scatter_view;
   int scatterSize;

   test_scatter_view_impl_cls(const scatter_view_type& view){
      scatter_view = view;
      scatterSize = 0;
   }

   void initialize(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        host_view(i, 0) = 999999.0;
        host_view(i, 1) = 999999.0;
        host_view(i, 2) = 999999.0;
      }
      Kokkos::fence();
      Kokkos::deep_copy(orig, host_view);
   }

   void run_parallel(int n) {
        scatterSize = n;
        auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
        Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(int i) const {
      auto scatter_access = scatter_view.access();
      auto scatter_access_atomic = scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
      for (int j = 0; j < 4; ++j) {
        auto k = (i + j) % scatterSize;
        scatter_access(k, 0).update((double)(j+1)*4);
        scatter_access_atomic(k, 1).update((double)(j+1)*2.0);
        scatter_access(k, 2).update((double)(j+1)*1.0);
      }
    }

    void validateResults(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        auto val0 = host_view(i, 0);
        auto val1 = host_view(i, 1);
        auto val2 = host_view(i, 2);
        EXPECT_TRUE(std::fabs((val0 - 4.0) / 4.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val1 - 2.0) / 2.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14);
      }
    }
};


template <typename ExecSpace, typename Layout, int duplication, int contribution>
struct test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, Kokkos::Experimental::ScatterMax>   
{
public:   

   typedef Kokkos::Experimental::ScatterView
       < double*[3]
       , Layout
       , ExecSpace
       , Kokkos::Experimental::ScatterMax
       , duplication
       , contribution
        > scatter_view_type;

   typedef Kokkos::View<double *[3], Layout, ExecSpace> orig_view_type; 


   scatter_view_type scatter_view;
   int scatterSize;

   test_scatter_view_impl_cls(const scatter_view_type& view){
      scatter_view = view;
      scatterSize = 0;
   }

   void initialize(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        host_view(i, 0) = 0.0;
        host_view(i, 1) = 0.0;
        host_view(i, 2) = 0.0;
      }
      Kokkos::fence();
      Kokkos::deep_copy(orig, host_view);
   }

   void run_parallel(int n) {
        scatterSize = n;
        auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
        Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
   }

   KOKKOS_INLINE_FUNCTION
   void operator()(int i) const {
      auto scatter_access = scatter_view.access();
      auto scatter_access_atomic = scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
      for (int j = 0; j < 4; ++j) {
        auto k = (i + j) % scatterSize;
        scatter_access(k, 0).update((double)(j+1)*4);
        scatter_access_atomic(k, 1).update((double)(j+1)*2.0);
        scatter_access(k, 2).update((double)(j+1)*1.0);
      }
    }

    void validateResults(orig_view_type orig) {
      auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
      Kokkos::fence();
      for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
        auto val0 = host_view(i, 0);
        auto val1 = host_view(i, 1);
        auto val2 = host_view(i, 2);
        EXPECT_TRUE(std::fabs((val0 - 16.0) / 16.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val1 - 8.0) / 8.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val2 - 4.0) / 4.0) < 1e-14);
      }
    }
};



template <typename ExecSpace, typename Layout, int duplication, int contribution, int op>
struct test_scatter_view_config
{
 public:
   typedef typename test_scatter_view_impl_cls<ExecSpace, Layout, 
         duplication, contribution, op>::scatter_view_type scatter_view_def;
   typedef typename test_scatter_view_impl_cls<ExecSpace, Layout, 
         duplication, contribution, op>::orig_view_type orig_view_def;

   test_scatter_view_config() {
   }

   void run_test(int n)
   {
     //Test creation via create_scatter_view
     {
     orig_view_def original_view("original_view", n);
     scatter_view_def scatter_view = Kokkos::Experimental::create_scatter_view
       < op
       , duplication
       , contribution
       > (original_view);

     test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, op> scatter_view_test_impl(scatter_view);
     scatter_view_test_impl.initialize(original_view);
     scatter_view_test_impl.run_parallel(n);

     Kokkos::Experimental::contribute(original_view, scatter_view);
     scatter_view.reset_except(original_view);

     scatter_view_test_impl.run_parallel(n);

     Kokkos::Experimental::contribute(original_view, scatter_view);
     Kokkos::fence();

     scatter_view_test_impl.validateResults(original_view);

     {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
     }
     }
     //Test creation via constructor
     {
     orig_view_def original_view("original_view", n);
     scatter_view_def scatter_view(original_view);

     test_scatter_view_impl_cls<ExecSpace, Layout, duplication, contribution, op> scatter_view_test_impl(scatter_view);
     scatter_view_test_impl.initialize(original_view);
     scatter_view_test_impl.run_parallel(n);

     Kokkos::Experimental::contribute(original_view, scatter_view);
     scatter_view.reset_except(original_view);

     scatter_view_test_impl.run_parallel(n);

     Kokkos::Experimental::contribute(original_view, scatter_view);
     Kokkos::fence();

     scatter_view_test_impl.validateResults(original_view);

     {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
     }
     }
   }

};


template <typename ExecSpace, int ScatterType>
struct TestDuplicatedScatterView {
  TestDuplicatedScatterView(int n) {
    // ScatterSum test
    test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
      Kokkos::Experimental::ScatterDuplicated,
      Kokkos::Experimental::ScatterNonAtomic,
      ScatterType> test_sv_right_config;
    test_sv_right_config.run_test(n);
    test_scatter_view_config<ExecSpace, Kokkos::LayoutLeft,
      Kokkos::Experimental::ScatterDuplicated,
      Kokkos::Experimental::ScatterNonAtomic,
      ScatterType> test_sv_left_config;
    test_sv_left_config.run_test(n);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
// disable duplicated instantiation with CUDA until
// UniqueToken can support it
template <int ScatterType>
struct TestDuplicatedScatterView<Kokkos::Cuda, ScatterType> {
  TestDuplicatedScatterView(int) {
  }
};
#endif

#ifdef KOKKOS_ENABLE_ROCM
// disable duplicated instantiation with ROCm until
// UniqueToken can support it
template <int ScatterType>
struct TestDuplicatedScatterView<Kokkos::Experimental::ROCm, ScatterType> {
  TestDuplicatedScatterView(int) {
  }
};
#endif

template <typename ExecSpace, int ScatterType>
void test_scatter_view(int n)
{
  // all of these configurations should compile okay, but only some of them are
  // correct and/or sensible in terms of memory use
  Kokkos::Experimental::UniqueToken<ExecSpace> unique_token{ExecSpace()};

  // no atomics or duplication is only sensible if the execution space
  // is running essentially in serial (doesn't have to be Serial though,
  // we also test OpenMP with one thread: LAMMPS cares about that)
  if (unique_token.size() == 1) {
    test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
      Kokkos::Experimental::ScatterNonDuplicated,
      Kokkos::Experimental::ScatterNonAtomic,
      ScatterType> test_sv_config;
    test_sv_config.run_test(n);
  }
#ifdef KOKKOS_ENABLE_SERIAL
  if (!std::is_same<ExecSpace, Kokkos::Serial>::value) {
#endif
  test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
    Kokkos::Experimental::ScatterNonDuplicated,
    Kokkos::Experimental::ScatterAtomic,
    ScatterType> test_sv_config;
  test_sv_config.run_test(n);
#ifdef KOKKOS_ENABLE_SERIAL
  }
#endif
  // with hundreds of threads we were running out of memory.
  // limit (n) so that duplication doesn't exceed 8GB
  constexpr std::size_t maximum_allowed_total_bytes = 8ull * 1024ull * 1024ull * 1024ull;
  std::size_t const maximum_allowed_copy_bytes = maximum_allowed_total_bytes / std::size_t(unique_token.size());
  constexpr std::size_t bytes_per_value = sizeof(double) * 3;
  std::size_t const maximum_allowed_copy_values = maximum_allowed_copy_bytes / bytes_per_value;
  n = std::min(n, int(maximum_allowed_copy_values));
  TestDuplicatedScatterView<ExecSpace, ScatterType> duptest(n);
}

TEST_F( TEST_CATEGORY, scatterview) {
#ifndef KOKKOS_ENABLE_ROCM
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterSum>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterProd>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMin>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMax>(10);
  // tests were timing out in DEBUG mode, reduce the amount of work
#ifdef KOKKOS_ENABLE_DEBUG
  int big_n = 100 * 1000;
#else
  int big_n = 10 * 1000 * 1000;
#endif
  test_scatter_view<TEST_EXECSPACE,Kokkos::Experimental::ScatterSum>(big_n);
  test_scatter_view<TEST_EXECSPACE,Kokkos::Experimental::ScatterProd>(big_n);
  test_scatter_view<TEST_EXECSPACE,Kokkos::Experimental::ScatterMin>(big_n);
  test_scatter_view<TEST_EXECSPACE,Kokkos::Experimental::ScatterMax>(big_n);
#endif
}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP


