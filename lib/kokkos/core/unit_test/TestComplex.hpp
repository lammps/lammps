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

#include<Kokkos_Core.hpp>
#include<cstdio>

namespace Test {

// Test construction and assignment

template<class ExecSpace>
struct TestComplexConstruction {
  Kokkos::View<Kokkos::complex<double>*,ExecSpace> d_results;
  typename Kokkos::View<Kokkos::complex<double>*,ExecSpace>::HostMirror h_results;
  
  void testit () {
    d_results = Kokkos::View<Kokkos::complex<double>*,ExecSpace>("TestComplexConstruction",10);
    h_results = Kokkos::create_mirror_view(d_results);
   
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,1), *this);
    Kokkos::fence();
    Kokkos::deep_copy(h_results,d_results);

    ASSERT_FLOAT_EQ(h_results(0).real(),1.5);  ASSERT_FLOAT_EQ(h_results(0).imag(),2.5);
    ASSERT_FLOAT_EQ(h_results(1).real(),1.5);  ASSERT_FLOAT_EQ(h_results(1).imag(),2.5);
    ASSERT_FLOAT_EQ(h_results(2).real(),0.0);  ASSERT_FLOAT_EQ(h_results(2).imag(),0.0);
    ASSERT_FLOAT_EQ(h_results(3).real(),3.5);  ASSERT_FLOAT_EQ(h_results(3).imag(),0.0);
    ASSERT_FLOAT_EQ(h_results(4).real(),4.5);  ASSERT_FLOAT_EQ(h_results(4).imag(),5.5);
    ASSERT_FLOAT_EQ(h_results(5).real(),1.5);  ASSERT_FLOAT_EQ(h_results(5).imag(),2.5);
    ASSERT_FLOAT_EQ(h_results(6).real(),4.5);  ASSERT_FLOAT_EQ(h_results(6).imag(),5.5);
    ASSERT_FLOAT_EQ(h_results(7).real(),7.5);  ASSERT_FLOAT_EQ(h_results(7).imag(),0.0);
    ASSERT_FLOAT_EQ(h_results(8).real(),double(8));  ASSERT_FLOAT_EQ(h_results(8).imag(),0.0);

#ifndef KOKKOS_ENABLE_ROCM
    Kokkos::complex<double> a(1.5,2.5),b(3.25,5.25),r_kk;
    std::complex<double> sa(a),sb(3.25,5.25),r;
    r = a; r_kk = a;         ASSERT_FLOAT_EQ(r.real(),r_kk.real()); ASSERT_FLOAT_EQ(r.imag(),r_kk.imag());
    r = sb*a; r_kk = b*a;    ASSERT_FLOAT_EQ(r.real(),r_kk.real()); ASSERT_FLOAT_EQ(r.imag(),r_kk.imag());
    r = sa; r_kk = a;        ASSERT_FLOAT_EQ(r.real(),r_kk.real()); ASSERT_FLOAT_EQ(r.imag(),r_kk.imag());
#endif

  }

  KOKKOS_INLINE_FUNCTION 
  void operator() (const int &i ) const {
    Kokkos::complex<double> a(1.5,2.5);
    d_results(0) = a;
    Kokkos::complex<double> b(a);
    d_results(1) = b;
    Kokkos::complex<double> c = Kokkos::complex<double>();
    d_results(2) = c;
    Kokkos::complex<double> d(3.5);
    d_results(3) = d; 
    volatile Kokkos::complex<double> a_v(4.5,5.5);
    d_results(4) = a_v;
    volatile Kokkos::complex<double> b_v(a);
    d_results(5) = b_v;
    Kokkos::complex<double> e(a_v);
    d_results(6) = e;

    d_results(7) = double(7.5);
    d_results(8) = int(8);
  } 
};

TEST_F(TEST_CATEGORY, complex_construction) {
  TestComplexConstruction<TEST_EXECSPACE> test;
  test.testit();
} 

// Test Math FUnction

template<class ExecSpace>
struct TestComplexBasicMath {
  Kokkos::View<Kokkos::complex<double>*,ExecSpace> d_results;
  typename Kokkos::View<Kokkos::complex<double>*,ExecSpace>::HostMirror h_results;

  void testit () {
    d_results = Kokkos::View<Kokkos::complex<double>*,ExecSpace>("TestComplexBasicMath",24);
    h_results = Kokkos::create_mirror_view(d_results);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,1), *this);
    Kokkos::fence();
    Kokkos::deep_copy(h_results,d_results);

    std::complex<double> a(1.5,2.5);
    std::complex<double> b(3.25,5.75);
    std::complex<double> d(1.0,2.0);
    double c = 9.3;
    int e = 2;

    std::complex<double> r;
    r = a+b; ASSERT_FLOAT_EQ(h_results(0).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(0).imag(),  r.imag());
    r = a-b; ASSERT_FLOAT_EQ(h_results(1).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(1).imag(),  r.imag());
    r = a*b; ASSERT_FLOAT_EQ(h_results(2).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(2).imag(),  r.imag());
    r = a/b; ASSERT_FLOAT_EQ(h_results(3).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(3).imag(),  r.imag());
    r = d+a; ASSERT_FLOAT_EQ(h_results(4).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(4).imag(),  r.imag());
    r = d-a; ASSERT_FLOAT_EQ(h_results(5).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(5).imag(),  r.imag());
    r = d*a; ASSERT_FLOAT_EQ(h_results(6).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(6).imag(),  r.imag());
    r = d/a; ASSERT_FLOAT_EQ(h_results(7).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(7).imag(),  r.imag());
    r = a+c; ASSERT_FLOAT_EQ(h_results(8).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(8).imag(),  r.imag());
    r = a-c; ASSERT_FLOAT_EQ(h_results(9).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(9).imag(),  r.imag());
    r = a*c; ASSERT_FLOAT_EQ(h_results(10).real(), r.real()); ASSERT_FLOAT_EQ(h_results(10).imag(), r.imag());
    r = a/c; ASSERT_FLOAT_EQ(h_results(11).real(), r.real()); ASSERT_FLOAT_EQ(h_results(11).imag(), r.imag());
    r = d+c; ASSERT_FLOAT_EQ(h_results(12).real(), r.real()); ASSERT_FLOAT_EQ(h_results(12).imag(), r.imag());
    r = d-c; ASSERT_FLOAT_EQ(h_results(13).real(), r.real()); ASSERT_FLOAT_EQ(h_results(13).imag(), r.imag());
    r = d*c; ASSERT_FLOAT_EQ(h_results(14).real(), r.real()); ASSERT_FLOAT_EQ(h_results(14).imag(), r.imag());
    r = d/c; ASSERT_FLOAT_EQ(h_results(15).real(), r.real()); ASSERT_FLOAT_EQ(h_results(15).imag(), r.imag());
    r = c+a; ASSERT_FLOAT_EQ(h_results(16).real(), r.real()); ASSERT_FLOAT_EQ(h_results(16).imag(), r.imag());
    r = c-a; ASSERT_FLOAT_EQ(h_results(17).real(), r.real()); ASSERT_FLOAT_EQ(h_results(17).imag(), r.imag());
    r = c*a; ASSERT_FLOAT_EQ(h_results(18).real(), r.real()); ASSERT_FLOAT_EQ(h_results(18).imag(), r.imag());
    r = c/a; ASSERT_FLOAT_EQ(h_results(19).real(), r.real()); ASSERT_FLOAT_EQ(h_results(19).imag(), r.imag());

    r = a; 
    /* r = a+e; */ ASSERT_FLOAT_EQ(h_results(20).real(),  r.real()+e); ASSERT_FLOAT_EQ(h_results(20).imag(),  r.imag());
    /* r = a-e; */ ASSERT_FLOAT_EQ(h_results(21).real(),  r.real()-e); ASSERT_FLOAT_EQ(h_results(21).imag(),  r.imag());
    /* r = a*e; */ ASSERT_FLOAT_EQ(h_results(22).real(),  r.real()*e); ASSERT_FLOAT_EQ(h_results(22).imag(),  r.imag()*e);
    /* r = a/e; */ ASSERT_FLOAT_EQ(h_results(23).real(),  r.real()/2); ASSERT_FLOAT_EQ(h_results(23).imag(),  r.imag()/e);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i ) const {
    Kokkos::complex<double> a(1.5,2.5);
    Kokkos::complex<double> b(3.25,5.75);
    // Basic math complex / complex
    d_results(0) = a+b;
    d_results(1) = a-b;
    d_results(2) = a*b;
    d_results(3) = a/b;
    d_results(4).real(1.0);
    d_results(4).imag(2.0);
    d_results(4) += a;
    d_results(5) = Kokkos::complex<double>(1.0,2.0);
    d_results(5) -= a;
    d_results(6) = Kokkos::complex<double>(1.0,2.0);
    d_results(6) *= a;
    d_results(7) = Kokkos::complex<double>(1.0,2.0);
    d_results(7) /= a;

    // Basic math complex / scalar
    double c = 9.3;
    d_results(8) = a+c;
    d_results(9) = a-c;
    d_results(10) = a*c;
    d_results(11) = a/c;
    d_results(12).real(1.0);
    d_results(12).imag(2.0);
    d_results(12) += c;
    d_results(13) = Kokkos::complex<double>(1.0,2.0);
    d_results(13) -= c;
    d_results(14) = Kokkos::complex<double>(1.0,2.0);
    d_results(14) *= c;
    d_results(15) = Kokkos::complex<double>(1.0,2.0);
    d_results(15) /= c;


    // Basic math scalar / complex
    d_results(16) = c+a;
    d_results(17) = c-a;
    d_results(18) = c*a;
    d_results(19) = c/a;

    int e = 2;
    d_results(20) = a+e;
    d_results(21) = a-e;
    d_results(22) = a*e;
    d_results(23) = a/e;
  }
};

TEST_F(TEST_CATEGORY, complex_basic_math) {
  TestComplexBasicMath<TEST_EXECSPACE> test;
  test.testit();
}


template<class ExecSpace>
struct TestComplexSpecialFunctions {
  Kokkos::View<Kokkos::complex<double>*,ExecSpace> d_results;
  typename Kokkos::View<Kokkos::complex<double>*,ExecSpace>::HostMirror h_results;

  void testit () {
    d_results = Kokkos::View<Kokkos::complex<double>*,ExecSpace>("TestComplexSpecialFunctions",20);
    h_results = Kokkos::create_mirror_view(d_results);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,1), *this);
    Kokkos::fence();
    Kokkos::deep_copy(h_results,d_results);

    std::complex<double> a(1.5,2.5);
    double c = 9.3;

    std::complex<double> r;
    r = a;             ASSERT_FLOAT_EQ(h_results(0).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(0).imag(),  r.imag());
    r = std::sqrt(a);  ASSERT_FLOAT_EQ(h_results(1).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(1).imag(),  r.imag());
    r = std::pow(a,c); ASSERT_FLOAT_EQ(h_results(2).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(2).imag(),  r.imag());
    r = std::abs(a);   ASSERT_FLOAT_EQ(h_results(3).real(),  r.real()); ASSERT_FLOAT_EQ(h_results(3).imag(),  r.imag());
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i ) const {
    Kokkos::complex<double> a(1.5,2.5);
    Kokkos::complex<double> b(3.25,5.75);
    double c = 9.3;

    d_results(0) = Kokkos::complex<double>(Kokkos::real(a),Kokkos::imag(a));
    d_results(1) = Kokkos::sqrt(a);
    d_results(2) = Kokkos::pow(a,c);
    d_results(3) = Kokkos::abs(a);

  }
};

TEST_F(TEST_CATEGORY, complex_special_funtions) {
  TestComplexSpecialFunctions<TEST_EXECSPACE> test;
  test.testit();
}
} // namespace Test


