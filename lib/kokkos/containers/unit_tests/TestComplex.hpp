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


#ifndef KOKKOS_TEST_COMPLEX_HPP
#define KOKKOS_TEST_COMPLEX_HPP

#include <Kokkos_Complex.hpp>
#include <gtest/gtest.h>
#include <iostream>

namespace Test {

namespace Impl {
  template <typename RealType>
  void testComplexConstructors () {
    typedef Kokkos::complex<RealType> complex_type;

    complex_type z1;
    complex_type z2 (0.0, 0.0);
    complex_type z3 (1.0, 0.0);
    complex_type z4 (0.0, 1.0);
    complex_type z5 (-1.0, -2.0);

    ASSERT_TRUE( z1 == z2 );
    ASSERT_TRUE( z1 != z3 );
    ASSERT_TRUE( z1 != z4 );
    ASSERT_TRUE( z1 != z5 );

    ASSERT_TRUE( z2 != z3 );
    ASSERT_TRUE( z2 != z4 );
    ASSERT_TRUE( z2 != z5 );

    ASSERT_TRUE( z3 != z4 );
    ASSERT_TRUE( z3 != z5 );

    complex_type z6 (-1.0, -2.0);
    ASSERT_TRUE( z5 == z6 );

    // Make sure that complex has value semantics, in particular, that
    // equality tests use values and not pointers, so that
    // reassignment actually changes the value.
    z1 = complex_type (-3.0, -4.0);
    ASSERT_TRUE( z1.real () == -3.0 );
    ASSERT_TRUE( z1.imag () == -4.0 );
    ASSERT_TRUE( z1 != z2 );

    complex_type z7 (1.0);
    ASSERT_TRUE( z3 == z7 );
    ASSERT_TRUE( z7 == 1.0 );
    ASSERT_TRUE( z7 != -1.0 );

    z7 = complex_type (5.0);
    ASSERT_TRUE( z7.real () == 5.0 );
    ASSERT_TRUE( z7.imag () == 0.0 );
  }

  template <typename RealType>
  void testPlus () {
    typedef Kokkos::complex<RealType> complex_type;

    complex_type z1 (1.0, -1.0);
    complex_type z2 (-1.0, 1.0);
    complex_type z3 = z1 + z2;
    ASSERT_TRUE( z3 == complex_type (0.0, 0.0) );
  }

  template <typename RealType>
  void testMinus () {
    typedef Kokkos::complex<RealType> complex_type;

    // Test binary minus.
    complex_type z1 (1.0, -1.0);
    complex_type z2 (-1.0, 1.0);
    complex_type z3 = z1 - z2;
    ASSERT_TRUE( z3 == complex_type (2.0, -2.0) );

    // Test unary minus.
    complex_type z4 (3.0, -4.0);
    ASSERT_TRUE( -z1 == complex_type (-3.0, 4.0) );
  }

  template <typename RealType>
  void testTimes () {
    typedef Kokkos::complex<RealType> complex_type;

    complex_type z1 (1.0, -1.0);
    complex_type z2 (-1.0, 1.0);
    complex_type z3 = z1 * z2;
    ASSERT_TRUE( z3 == complex_type (0.0, 2.0) );

    // Make sure that std::complex * Kokkos::complex works too.
    std::complex<RealType> z4 (-1.0, 1.0);
    complex_type z5 = z4 * z1;
    ASSERT_TRUE( z5 == complex_type (0.0, 2.0) );
  }

  template <typename RealType>
  void testDivide () {
    typedef Kokkos::complex<RealType> complex_type;

    // Test division of a complex number by a real number.
    complex_type z1 (1.0, -1.0);
    complex_type z2 (1.0 / 2.0, -1.0 / 2.0);
    ASSERT_TRUE( z1 / 2.0 == z2 );

    // (-1+2i)/(1-i) == ((-1+2i)(1+i)) / ((1-i)(1+i))
    // (-1+2i)(1+i) == -3 + i
    complex_type z3 (-1.0, 2.0);
    complex_type z4 (1.0, -1.0);
    complex_type z5 (-3.0, 1.0);
    ASSERT_TRUE(z3 * Kokkos::conj (z4) == z5 );

    // Test division of a complex number by a complex number.
    // This assumes that RealType is a floating-point type.
    complex_type z6 (Kokkos::real (z5) / 2.0,
                     Kokkos::imag (z5) / 2.0);

    complex_type z7 = z3 / z4;
    ASSERT_TRUE( z7 == z6 );
  }

  template <typename RealType>
  void testOutsideKernel () {
    testComplexConstructors<RealType> ();
    testPlus<RealType> ();
    testTimes<RealType> ();
    testDivide<RealType> ();
  }


  template<typename RealType, typename Device>
  void testCreateView () {
    typedef Kokkos::complex<RealType> complex_type;
    Kokkos::View<complex_type*, Device> x ("x", 10);
    ASSERT_TRUE( x.dimension_0 () == 10 );

    // Test that View assignment works.
    Kokkos::View<complex_type*, Device> x_nonconst = x;
    Kokkos::View<const complex_type*, Device> x_const = x;
  }

  template<typename RealType, typename Device>
  class Fill {
  public:
    typedef typename Device::execution_space execution_space;

    typedef Kokkos::View<Kokkos::complex<RealType>*, Device> view_type;
    typedef typename view_type::size_type size_type;

    KOKKOS_INLINE_FUNCTION
    void operator () (const size_type i) const {
      x_(i) = val_;
    }

    Fill (const view_type& x, const Kokkos::complex<RealType>& val) :
      x_ (x), val_ (val)
    {}

  private:
    view_type x_;
    const Kokkos::complex<RealType> val_;
  };

  template<typename RealType, typename Device>
  class Sum {
  public:
    typedef typename Device::execution_space execution_space;

    typedef Kokkos::View<const Kokkos::complex<RealType>*, Device> view_type;
    typedef typename view_type::size_type size_type;
    typedef Kokkos::complex<RealType> value_type;

    KOKKOS_INLINE_FUNCTION
    void operator () (const size_type i, Kokkos::complex<RealType>& sum) const {
      sum += x_(i);
    }

    Sum (const view_type& x) : x_ (x) {}

  private:
    view_type x_;
  };

  template<typename RealType, typename Device>
  void testInsideKernel () {
    typedef Kokkos::complex<RealType> complex_type;
    typedef Kokkos::View<complex_type*, Device> view_type;
    typedef typename view_type::size_type size_type;

    const size_type N = 1000;
    view_type x ("x", N);
    ASSERT_TRUE( x.dimension_0 () == N );

    // Kokkos::parallel_reduce (N, [=] (const size_type i, complex_type& result) {
    //     result += x[i];
    //   });

    Kokkos::parallel_for (N, Fill<RealType, Device> (x, complex_type (1.0, -1.0)));

    complex_type sum;
    Kokkos::parallel_reduce (N, Sum<RealType, Device> (x), sum);

    ASSERT_TRUE( sum.real () == 1000.0 && sum.imag () == -1000.0 );
  }
} // namespace Impl


template <typename Device>
void testComplex ()
{
  Impl::testOutsideKernel<float> ();
  Impl::testOutsideKernel<double> ();

  Impl::testCreateView<float, Device> ();
  Impl::testCreateView<double, Device> ();

  Impl::testInsideKernel<float, Device> ();
  Impl::testInsideKernel<double, Device> ();
}


} // namespace Test

#endif // KOKKOS_TEST_COMPLEX_HPP
