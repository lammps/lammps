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
#ifndef KOKKOS_COMPLEX_HPP
#define KOKKOS_COMPLEX_HPP

#include <Kokkos_Atomic.hpp>
#include <complex>
#include <iostream>

namespace Kokkos {

/// \class complex
/// \brief Partial reimplementation of std::complex that works as the
///   result of a Kokkos::parallel_reduce.
/// \tparam RealType The type of the real and imaginary parts of the
///   complex number.  As with std::complex, this is only defined for
///   \c float, \c double, and <tt>long double</tt>.  The latter is
///   currently forbidden in CUDA device kernels.
template<class RealType>
class complex {
private:
  RealType re_, im_;

public:
  //! The type of the real or imaginary parts of this complex number.
  typedef RealType value_type;

  //! Default constructor (initializes both real and imaginary parts to zero).
  KOKKOS_INLINE_FUNCTION complex () :
    re_ (0.0), im_ (0.0)
  {}

  //! Copy constructor.
  KOKKOS_INLINE_FUNCTION complex (const complex<RealType>& src) :
    re_ (src.re_), im_ (src.im_)
  {}

  //! Copy constructor from volatile.
  KOKKOS_INLINE_FUNCTION complex (const volatile complex<RealType>& src) :
    re_ (src.re_), im_ (src.im_)
  {}

  /// \brief Conversion constructor from std::complex.
  ///
  /// This constructor cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  template<class InputRealType>
  complex (const std::complex<InputRealType>& src) :
    re_ (std::real (src)), im_ (std::imag (src))
  {}

  /// \brief Conversion operator to std::complex.
  ///
  /// This operator cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  operator std::complex<RealType> () const {
    return std::complex<RealType> (re_, im_);
  }

  /// \brief Constructor that takes just the real part, and sets the
  ///   imaginary part to zero.
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION complex (const InputRealType& val) :
    re_ (val), im_ (0.0)
  {}

  //! Constructor that takes the real and imaginary parts.
  template<class RealType1, class RealType2>
  KOKKOS_INLINE_FUNCTION complex (const RealType1& re, const RealType2& im) :
    re_ (re), im_ (im)
  {}

  //! Assignment operator.
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const complex<InputRealType>& src) {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  //! Assignment operator.
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  volatile complex<RealType>& operator= (const complex<InputRealType>& src) volatile {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  //! Assignment operator.
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  volatile complex<RealType>& operator= (const volatile complex<InputRealType>& src) volatile {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  //! Assignment operator.
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const volatile complex<InputRealType>& src) {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  //! Assignment operator (from a real number).
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator= (const InputRealType& val) {
    re_ = val;
    im_ = static_cast<RealType> (0.0);
    return *this;
  }

  //! Assignment operator (from a real number).
  template<class InputRealType>
  KOKKOS_INLINE_FUNCTION
  void operator= (const InputRealType& val) volatile {
    re_ = val;
    im_ = static_cast<RealType> (0.0);
  }

  /// \brief Assignment operator from std::complex.
  ///
  /// This constructor cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  template<class InputRealType>
  complex<RealType>& operator= (const std::complex<InputRealType>& src) {
    re_ = std::real (src);
    im_ = std::imag (src);
    return *this;
  }

  //! The imaginary part of this complex number.
  KOKKOS_INLINE_FUNCTION RealType& imag () {
    return im_;
  }

  //! The real part of this complex number.
  KOKKOS_INLINE_FUNCTION RealType& real () {
    return re_;
  }

  //! The imaginary part of this complex number.
  KOKKOS_INLINE_FUNCTION const RealType imag () const {
    return im_;
  }

  //! The real part of this complex number.
  KOKKOS_INLINE_FUNCTION const RealType real () const {
    return re_;
  }

  //! The imaginary part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION volatile RealType& imag () volatile {
    return im_;
  }

  //! The real part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION volatile RealType& real () volatile {
    return re_;
  }

  //! The imaginary part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION const RealType imag () const volatile {
    return im_;
  }

  //! The real part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION const RealType real () const volatile {
    return re_;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator += (const complex<RealType>& src) {
    re_ += src.re_;
    im_ += src.im_;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile complex<RealType>& src) volatile {
    re_ += src.re_;
    im_ += src.im_;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator += (const RealType& src) {
    re_ += src;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile RealType& src) volatile {
    re_ += src;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator -= (const complex<RealType>& src) {
    re_ -= src.re_;
    im_ -= src.im_;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator -= (const RealType& src) {
    re_ -= src;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator *= (const complex<RealType>& src) {
    const RealType realPart = re_ * src.re_ - im_ * src.im_;
    const RealType imagPart = re_ * src.im_ + im_ * src.re_;
    re_ = realPart;
    im_ = imagPart;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator *= (const volatile complex<RealType>& src) volatile {
    const RealType realPart = re_ * src.re_ - im_ * src.im_;
    const RealType imagPart = re_ * src.im_ + im_ * src.re_;
    re_ = realPart;
    im_ = imagPart;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator *= (const RealType& src) {
    re_ *= src;
    im_ *= src;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator *= (const volatile RealType& src) volatile {
    re_ *= src;
    im_ *= src;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator /= (const complex<RealType>& y) {
    // Scale (by the "1-norm" of y) to avoid unwarranted overflow.
    // If the real part is +/-Inf and the imaginary part is -/+Inf,
    // this won't change the result.
    const RealType s = ::fabs (y.real ()) + ::fabs (y.imag ());

    // If s is 0, then y is zero, so x/y == real(x)/0 + i*imag(x)/0.
    // In that case, the relation x/y == (x/s) / (y/s) doesn't hold,
    // because y/s is NaN.
    if (s == 0.0) {
      this->re_ /= s;
      this->im_ /= s;
    }
    else {
      const complex<RealType> x_scaled (this->re_ / s, this->im_ / s);
      const complex<RealType> y_conj_scaled (y.re_ / s, -(y.im_) / s);
      const RealType y_scaled_abs = y_conj_scaled.re_ * y_conj_scaled.re_ +
        y_conj_scaled.im_ * y_conj_scaled.im_; // abs(y) == abs(conj(y))
      *this = x_scaled * y_conj_scaled;
      *this /= y_scaled_abs;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  complex<RealType>& operator /= (const RealType& src) {
    re_ /= src;
    im_ /= src;
    return *this;
  }
};

//! Binary + operator for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator + (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () + y.real (), x.imag () + y.imag ());
}

//! Unary + operator for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator + (const complex<RealType>& x) {
  return x;
}

//! Binary - operator for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator - (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () - y.real (), x.imag () - y.imag ());
}

//! Unary - operator for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator - (const complex<RealType>& x) {
  return complex<RealType> (-x.real (), -x.imag ());
}

//! Binary * operator for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator * (const complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () * y.real () - x.imag () * y.imag (),
                            x.real () * y.imag () + x.imag () * y.real ());
}

/// \brief Binary * operator for std::complex and complex.
///
/// This function exists because GCC 4.7.2 (and perhaps other
/// compilers) are not able to deduce that they can multiply
/// std::complex by Kokkos::complex, by first converting std::complex
/// to Kokkos::complex.
///
/// This function cannot be called in a CUDA device function, because
/// std::complex's methods and nonmember functions are not marked as
/// CUDA device functions.
template<class RealType>
complex<RealType>
operator * (const std::complex<RealType>& x, const complex<RealType>& y) {
  return complex<RealType> (x.real () * y.real () - x.imag () * y.imag (),
                            x.real () * y.imag () + x.imag () * y.real ());
}

/// \brief Binary * operator for RealType times complex.
///
/// This function exists because the compiler doesn't know that
/// RealType and complex<RealType> commute with respect to operator*.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator * (const RealType& x, const complex<RealType>& y) {
  return complex<RealType> (x * y.real (), x * y.imag ());
}


//! Imaginary part of a complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
RealType imag (const complex<RealType>& x) {
  return x.imag ();
}

//! Real part of a complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
RealType real (const complex<RealType>& x) {
  return x.real ();
}

//! Absolute value (magnitude) of a complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
RealType abs (const complex<RealType>& x) {
  // FIXME (mfh 31 Oct 2014) Scale to avoid unwarranted overflow.
  return ::sqrt (real (x) * real (x) + imag (x) * imag (x));
}

//! Conjugate of a complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType> conj (const complex<RealType>& x) {
  return complex<RealType> (real (x), -imag (x));
}


//! Binary operator / for complex and real numbers
template<class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
complex<RealType1>
operator / (const complex<RealType1>& x, const RealType2& y) {
  return complex<RealType1> (real (x) / y, imag (x) / y);
}

//! Binary operator / for complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
complex<RealType>
operator / (const complex<RealType>& x, const complex<RealType>& y) {
  // Scale (by the "1-norm" of y) to avoid unwarranted overflow.
  // If the real part is +/-Inf and the imaginary part is -/+Inf,
  // this won't change the result.
  const RealType s = ::fabs (real (y)) + ::fabs (imag (y));

  // If s is 0, then y is zero, so x/y == real(x)/0 + i*imag(x)/0.
  // In that case, the relation x/y == (x/s) / (y/s) doesn't hold,
  // because y/s is NaN.
  if (s == 0.0) {
    return complex<RealType> (real (x) / s, imag (x) / s);
  }
  else {
    const complex<RealType> x_scaled (real (x) / s, imag (x) / s);
    const complex<RealType> y_conj_scaled (real (y) / s, -imag (y) / s);
    const RealType y_scaled_abs = real (y_conj_scaled) * real (y_conj_scaled) +
      imag (y_conj_scaled) * imag (y_conj_scaled); // abs(y) == abs(conj(y))
    complex<RealType> result = x_scaled * y_conj_scaled;
    result /= y_scaled_abs;
    return result;
  }
}

//! Equality operator for two complex numbers.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator == (const complex<RealType>& x, const complex<RealType>& y) {
  return real (x) == real (y) && imag (x) == imag (y);
}

//! Equality operator for std::complex and Kokkos::complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator == (const std::complex<RealType>& x, const complex<RealType>& y) {
  return std::real (x) == real (y) && std::imag (x) == imag (y);
}

//! Equality operator for complex and real number.
template<class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
bool operator == (const complex<RealType1>& x, const RealType2& y) {
  return real (x) == y && imag (x) == static_cast<RealType1> (0.0);
}

//! Equality operator for real and complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator == (const RealType& x, const complex<RealType>& y) {
  return y == x;
}

//! Inequality operator for two complex numbers.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const complex<RealType>& x, const complex<RealType>& y) {
  return real (x) != real (y) || imag (x) != imag (y);
}

//! Inequality operator for std::complex and Kokkos::complex.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const std::complex<RealType>& x, const complex<RealType>& y) {
  return std::real (x) != real (y) || std::imag (x) != imag (y);
}

//! Inequality operator for complex and real number.
template<class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
bool operator != (const complex<RealType1>& x, const RealType2& y) {
  return real (x) != y || imag (x) != static_cast<RealType1> (0.0);
}

//! Inequality operator for real and complex number.
template<class RealType>
KOKKOS_INLINE_FUNCTION
bool operator != (const RealType& x, const complex<RealType>& y) {
  return y != x;
}

template<class RealType>
std::ostream& operator << (std::ostream& os, const complex<RealType>& x) {
  const std::complex<RealType> x_std (Kokkos::real (x), Kokkos::imag (x));
  os << x_std;
  return os;
}

template<class RealType>
std::ostream& operator >> (std::ostream& os, complex<RealType>& x) {
  std::complex<RealType> x_std;
  os >> x_std;
  x = x_std; // only assigns on success of above
  return os;
}


} // namespace Kokkos

#endif // KOKKOS_COMPLEX_HPP
