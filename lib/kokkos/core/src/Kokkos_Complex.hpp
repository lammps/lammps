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
#ifndef KOKKOS_COMPLEX_HPP
#define KOKKOS_COMPLEX_HPP

#include <Kokkos_Atomic.hpp>
#include <Kokkos_NumericTraits.hpp>
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
template <class RealType>
class
#ifdef KOKKOS_ENABLE_COMPLEX_ALIGN
    alignas(2 * sizeof(RealType))
#endif
        complex {
 private:
  RealType re_{};
  RealType im_{};

 public:
  //! The type of the real or imaginary parts of this complex number.
  using value_type = RealType;

  //! Default constructor (initializes both real and imaginary parts to zero).
  KOKKOS_DEFAULTED_FUNCTION
  complex() noexcept = default;

  //! Copy constructor.
  KOKKOS_DEFAULTED_FUNCTION
  complex(const complex&) noexcept = default;

  KOKKOS_DEFAULTED_FUNCTION
  complex& operator=(const complex&) noexcept = default;

  /// \brief Conversion constructor from compatible RType
  template <class RType,
            typename std::enable_if<std::is_convertible<RType, RealType>::value,
                                    int>::type = 0>
  KOKKOS_INLINE_FUNCTION complex(const complex<RType>& other) noexcept
      // Intentionally do the conversions implicitly here so that users don't
      // get any warnings about narrowing, etc., that they would expect to get
      // otherwise.
      : re_(other.real()), im_(other.imag()) {}

  /// \brief Conversion constructor from std::complex.
  ///
  /// This constructor cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  KOKKOS_INLINE_FUNCTION
  complex(const std::complex<RealType>& src) noexcept
      // We can use this aspect of the standard to avoid calling
      // non-device-marked functions `std::real` and `std::imag`: "For any
      // object z of type complex<T>, reinterpret_cast<T(&)[2]>(z)[0] is the
      // real part of z and reinterpret_cast<T(&)[2]>(z)[1] is the imaginary
      // part of z." Now we don't have to provide a whole bunch of the overloads
      // of things taking either Kokkos::complex or std::complex
      : re_(reinterpret_cast<const RealType (&)[2]>(src)[0]),
        im_(reinterpret_cast<const RealType (&)[2]>(src)[1]) {}

  /// \brief Conversion operator to std::complex.
  ///
  /// This operator cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  // TODO: make explicit.  DJS 2019-08-28
  operator std::complex<RealType>() const noexcept {
    return std::complex<RealType>(re_, im_);
  }

  /// \brief Constructor that takes just the real part, and sets the
  ///   imaginary part to zero.
  KOKKOS_INLINE_FUNCTION complex(const RealType& val) noexcept
      : re_(val), im_(static_cast<RealType>(0)) {}

  //! Constructor that takes the real and imaginary parts.
  KOKKOS_INLINE_FUNCTION
  complex(const RealType& re, const RealType& im) noexcept : re_(re), im_(im) {}

  //! Assignment operator (from a real number).
  KOKKOS_INLINE_FUNCTION complex& operator=(const RealType& val) noexcept {
    re_ = val;
    im_ = RealType(0);
    return *this;
  }

  /// \brief Assignment operator from std::complex.
  ///
  /// This constructor cannot be called in a CUDA device function,
  /// because std::complex's methods and nonmember functions are not
  /// marked as CUDA device functions.
  complex& operator=(const std::complex<RealType>& src) noexcept {
    *this = complex(src);
    return *this;
  }

  //! The imaginary part of this complex number.
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14 RealType& imag() noexcept { return im_; }

  //! The real part of this complex number.
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14 RealType& real() noexcept { return re_; }

  //! The imaginary part of this complex number.
  KOKKOS_INLINE_FUNCTION
  constexpr RealType imag() const noexcept { return im_; }

  //! The real part of this complex number.
  KOKKOS_INLINE_FUNCTION
  constexpr RealType real() const noexcept { return re_; }

  //! Set the imaginary part of this complex number.
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  void imag(RealType v) noexcept { im_ = v; }

  //! Set the real part of this complex number.
  KOKKOS_INLINE_FUNCTION
  KOKKOS_CONSTEXPR_14
  void real(RealType v) noexcept { re_ = v; }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator+=(
      const complex<RealType>& src) noexcept {
    re_ += src.re_;
    im_ += src.im_;
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator+=(
      const RealType& src) noexcept {
    re_ += src;
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator-=(
      const complex<RealType>& src) noexcept {
    re_ -= src.re_;
    im_ -= src.im_;
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator-=(
      const RealType& src) noexcept {
    re_ -= src;
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator*=(
      const complex<RealType>& src) noexcept {
    const RealType realPart = re_ * src.re_ - im_ * src.im_;
    const RealType imagPart = re_ * src.im_ + im_ * src.re_;
    re_                     = realPart;
    im_                     = imagPart;
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator*=(
      const RealType& src) noexcept {
    re_ *= src;
    im_ *= src;
    return *this;
  }

  // Conditional noexcept, just in case RType throws on divide-by-zero
  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator/=(
      const complex<RealType>& y) noexcept(noexcept(RealType{} / RealType{})) {
    // Scale (by the "1-norm" of y) to avoid unwarranted overflow.
    // If the real part is +/-Inf and the imaginary part is -/+Inf,
    // this won't change the result.
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
    using std::fabs;
#endif
    const RealType s = fabs(y.real()) + fabs(y.imag());

    // If s is 0, then y is zero, so x/y == real(x)/0 + i*imag(x)/0.
    // In that case, the relation x/y == (x/s) / (y/s) doesn't hold,
    // because y/s is NaN.
    // TODO mark this branch unlikely
    if (s == RealType(0)) {
      this->re_ /= s;
      this->im_ /= s;
    } else {
      const complex x_scaled(this->re_ / s, this->im_ / s);
      const complex y_conj_scaled(y.re_ / s, -(y.im_) / s);
      const RealType y_scaled_abs =
          y_conj_scaled.re_ * y_conj_scaled.re_ +
          y_conj_scaled.im_ * y_conj_scaled.im_;  // abs(y) == abs(conj(y))
      *this = x_scaled * y_conj_scaled;
      *this /= y_scaled_abs;
    }
    return *this;
  }

  KOKKOS_CONSTEXPR_14
  KOKKOS_INLINE_FUNCTION complex& operator/=(
      const std::complex<RealType>& y) noexcept(noexcept(RealType{} /
                                                         RealType{})) {
    // Scale (by the "1-norm" of y) to avoid unwarranted overflow.
    // If the real part is +/-Inf and the imaginary part is -/+Inf,
    // this won't change the result.
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
    using std::fabs;
#endif
    const RealType s = fabs(y.real()) + fabs(y.imag());

    // If s is 0, then y is zero, so x/y == real(x)/0 + i*imag(x)/0.
    // In that case, the relation x/y == (x/s) / (y/s) doesn't hold,
    // because y/s is NaN.
    if (s == RealType(0)) {
      this->re_ /= s;
      this->im_ /= s;
    } else {
      const complex x_scaled(this->re_ / s, this->im_ / s);
      const complex y_conj_scaled(y.re_ / s, -(y.im_) / s);
      const RealType y_scaled_abs =
          y_conj_scaled.re_ * y_conj_scaled.re_ +
          y_conj_scaled.im_ * y_conj_scaled.im_;  // abs(y) == abs(conj(y))
      *this = x_scaled * y_conj_scaled;
      *this /= y_scaled_abs;
    }
    return *this;
  }

  KOKKOS_CONSTEXPR_14 KOKKOS_INLINE_FUNCTION complex& operator/=(
      const RealType& src) noexcept(noexcept(RealType{} / RealType{})) {
    re_ /= src;
    im_ /= src;
    return *this;
  }

  //---------------------------------------------------------------------------
  // TODO: refactor Kokkos reductions to remove dependency on
  // volatile member overloads since they are being deprecated in c++20
  //---------------------------------------------------------------------------

  //! Copy constructor from volatile.
  template <class RType,
            typename std::enable_if<std::is_convertible<RType, RealType>::value,
                                    int>::type = 0>
  KOKKOS_INLINE_FUNCTION complex(const volatile complex<RType>& src) noexcept
      // Intentionally do the conversions implicitly here so that users don't
      // get any warnings about narrowing, etc., that they would expect to get
      // otherwise.
      : re_(src.re_), im_(src.im_) {}

  /// \brief Assignment operator, for volatile <tt>*this</tt> and
  ///   nonvolatile input.
  ///
  /// \param src [in] Input; right-hand side of the assignment.
  ///
  /// This operator returns \c void instead of <tt>volatile
  /// complex& </tt>.  See Kokkos Issue #177 for the
  /// explanation.  In practice, this means that you should not chain
  /// assignments with volatile lvalues.
  //
  // Templated, so as not to be a copy assignment operator (Kokkos issue #2577)
  // Intended to behave as
  //    void operator=(const complex&) volatile noexcept
  //
  // Use cases:
  //    complex r;
  //    const complex cr;
  //    volatile complex vl;
  //    vl = r;
  //    vl = cr;
  template <class Complex,
            typename std::enable_if<std::is_same<Complex, complex>::value,
                                    int>::type = 0>
  KOKKOS_INLINE_FUNCTION void operator=(const Complex& src) volatile noexcept {
    re_ = src.re_;
    im_ = src.im_;
    // We deliberately do not return anything here.  See explanation
    // in public documentation above.
  }

  //! Assignment operator, volatile LHS and volatile RHS
  // TODO Should this return void like the other volatile assignment operators?
  //
  // Templated, so as not to be a copy assignment operator (Kokkos issue #2577)
  // Intended to behave as
  //    volatile complex& operator=(const volatile complex&) volatile noexcept
  //
  // Use cases:
  //    volatile complex vr;
  //    const volatile complex cvr;
  //    volatile complex vl;
  //    vl = vr;
  //    vl = cvr;
  template <class Complex,
            typename std::enable_if<std::is_same<Complex, complex>::value,
                                    int>::type = 0>
  KOKKOS_INLINE_FUNCTION volatile complex& operator=(
      const volatile Complex& src) volatile noexcept {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  //! Assignment operator, volatile RHS and non-volatile LHS
  //
  // Templated, so as not to be a copy assignment operator (Kokkos issue #2577)
  // Intended to behave as
  //    complex& operator=(const volatile complex&) noexcept
  //
  // Use cases:
  //    volatile complex vr;
  //    const volatile complex cvr;
  //    complex l;
  //    l = vr;
  //    l = cvr;
  //
  template <class Complex,
            typename std::enable_if<std::is_same<Complex, complex>::value,
                                    int>::type = 0>
  KOKKOS_INLINE_FUNCTION complex& operator=(
      const volatile Complex& src) noexcept {
    re_ = src.re_;
    im_ = src.im_;
    return *this;
  }

  // Mirroring the behavior of the assignment operators from complex RHS in the
  // RealType RHS versions.

  //! Assignment operator (from a volatile real number).
  KOKKOS_INLINE_FUNCTION void operator=(const volatile RealType& val) noexcept {
    re_ = val;
    im_ = RealType(0);
    // We deliberately do not return anything here.  See explanation
    // in public documentation above.
  }

  //! Assignment operator volatile LHS and non-volatile RHS
  KOKKOS_INLINE_FUNCTION complex& operator=(
      const RealType& val) volatile noexcept {
    re_ = val;
    im_ = RealType(0);
    return *this;
  }

  //! Assignment operator volatile LHS and volatile RHS
  // TODO Should this return void like the other volatile assignment operators?
  KOKKOS_INLINE_FUNCTION complex& operator=(
      const volatile RealType& val) volatile noexcept {
    re_ = val;
    im_ = RealType(0);
    return *this;
  }

  //! The imaginary part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION
  volatile RealType& imag() volatile noexcept { return im_; }

  //! The real part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION
  volatile RealType& real() volatile noexcept { return re_; }

  //! The imaginary part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION
  RealType imag() const volatile noexcept { return im_; }

  //! The real part of this complex number (volatile overload).
  KOKKOS_INLINE_FUNCTION
  RealType real() const volatile noexcept { return re_; }

  KOKKOS_INLINE_FUNCTION void operator+=(
      const volatile complex<RealType>& src) volatile noexcept {
    re_ += src.re_;
    im_ += src.im_;
  }

  KOKKOS_INLINE_FUNCTION void operator+=(
      const volatile RealType& src) volatile noexcept {
    re_ += src;
  }

  KOKKOS_INLINE_FUNCTION void operator*=(
      const volatile complex<RealType>& src) volatile noexcept {
    const RealType realPart = re_ * src.re_ - im_ * src.im_;
    const RealType imagPart = re_ * src.im_ + im_ * src.re_;

    re_ = realPart;
    im_ = imagPart;
  }

  KOKKOS_INLINE_FUNCTION void operator*=(
      const volatile RealType& src) volatile noexcept {
    re_ *= src;
    im_ *= src;
  }

  // TODO DSH 2019-10-7 why are there no volatile /= and friends?
};

//==============================================================================
// <editor-fold desc="Equality and inequality"> {{{1

// Note that this is not the same behavior as std::complex, which doesn't allow
// implicit conversions, but since this is the way we had it before, we have
// to do it this way now.

//! Binary == operator for complex complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION bool operator==(complex<RealType1> const& x,
                                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) == common_type(y.real()) &&
         common_type(x.imag()) == common_type(y.imag());
}

// TODO (here and elsewhere) decide if we should convert to a Kokkos::complex
//      and do the comparison in a device-marked function
//! Binary == operator for std::complex complex.
template <class RealType1, class RealType2>
inline bool operator==(std::complex<RealType1> const& x,
                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) == common_type(y.real()) &&
         common_type(x.imag()) == common_type(y.imag());
}

//! Binary == operator for complex std::complex.
template <class RealType1, class RealType2>
inline bool operator==(complex<RealType1> const& x,
                       std::complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) == common_type(y.real()) &&
         common_type(x.imag()) == common_type(y.imag());
}

//! Binary == operator for complex real.
template <
    class RealType1, class RealType2,
    // Constraints to avoid participation in oparator==() for every possible RHS
    typename std::enable_if<std::is_convertible<RealType2, RealType1>::value,
                            int>::type = 0>
KOKKOS_INLINE_FUNCTION bool operator==(complex<RealType1> const& x,
                                       RealType2 const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) == common_type(y) &&
         common_type(x.imag()) == common_type(0);
}

//! Binary == operator for real complex.
template <
    class RealType1, class RealType2,
    // Constraints to avoid participation in oparator==() for every possible RHS
    typename std::enable_if<std::is_convertible<RealType1, RealType2>::value,
                            int>::type = 0>
KOKKOS_INLINE_FUNCTION bool operator==(RealType1 const& x,
                                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x) == common_type(y.real()) &&
         common_type(0) == common_type(y.imag());
}

//! Binary != operator for complex complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION bool operator!=(complex<RealType1> const& x,
                                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) != common_type(y.real()) ||
         common_type(x.imag()) != common_type(y.imag());
}

//! Binary != operator for std::complex complex.
template <class RealType1, class RealType2>
inline bool operator!=(std::complex<RealType1> const& x,
                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) != common_type(y.real()) ||
         common_type(x.imag()) != common_type(y.imag());
}

//! Binary != operator for complex std::complex.
template <class RealType1, class RealType2>
inline bool operator!=(complex<RealType1> const& x,
                       std::complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) != common_type(y.real()) ||
         common_type(x.imag()) != common_type(y.imag());
}

//! Binary != operator for complex real.
template <
    class RealType1, class RealType2,
    // Constraints to avoid participation in oparator==() for every possible RHS
    typename std::enable_if<std::is_convertible<RealType2, RealType1>::value,
                            int>::type = 0>
KOKKOS_INLINE_FUNCTION bool operator!=(complex<RealType1> const& x,
                                       RealType2 const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x.real()) != common_type(y) ||
         common_type(x.imag()) != common_type(0);
}

//! Binary != operator for real complex.
template <
    class RealType1, class RealType2,
    // Constraints to avoid participation in oparator==() for every possible RHS
    typename std::enable_if<std::is_convertible<RealType1, RealType2>::value,
                            int>::type = 0>
KOKKOS_INLINE_FUNCTION bool operator!=(RealType1 const& x,
                                       complex<RealType2> const& y) noexcept {
  using common_type = typename std::common_type<RealType1, RealType2>::type;
  return common_type(x) != common_type(y.real()) ||
         common_type(0) != common_type(y.imag());
}

// </editor-fold> end Equality and inequality }}}1
//==============================================================================

//! Binary + operator for complex complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator+(const complex<RealType1>& x,
              const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() + y.real(), x.imag() + y.imag());
}

//! Binary + operator for complex scalar.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator+(const complex<RealType1>& x, const RealType2& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() + y, x.imag());
}

//! Binary + operator for scalar complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator+(const RealType1& x, const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x + y.real(), y.imag());
}

//! Unary + operator for complex.
template <class RealType>
KOKKOS_INLINE_FUNCTION complex<RealType> operator+(
    const complex<RealType>& x) noexcept {
  return complex<RealType>{+x.real(), +x.imag()};
}

//! Binary - operator for complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator-(const complex<RealType1>& x,
              const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() - y.real(), x.imag() - y.imag());
}

//! Binary - operator for complex scalar.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator-(const complex<RealType1>& x, const RealType2& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() - y, x.imag());
}

//! Binary - operator for scalar complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator-(const RealType1& x, const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x - y.real(), -y.imag());
}

//! Unary - operator for complex.
template <class RealType>
KOKKOS_INLINE_FUNCTION complex<RealType> operator-(
    const complex<RealType>& x) noexcept {
  return complex<RealType>(-x.real(), -x.imag());
}

//! Binary * operator for complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator*(const complex<RealType1>& x,
              const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() * y.real() - x.imag() * y.imag(),
      x.real() * y.imag() + x.imag() * y.real());
}

/// \brief Binary * operator for std::complex and complex.
///
/// This needs to exist because template parameters can't be deduced when
/// conversions occur.  We could probably fix this using hidden friends patterns
///
/// This function cannot be called in a CUDA device function, because
/// std::complex's methods and nonmember functions are not marked as
/// CUDA device functions.
template <class RealType1, class RealType2>
inline complex<typename std::common_type<RealType1, RealType2>::type> operator*(
    const std::complex<RealType1>& x, const complex<RealType2>& y) {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x.real() * y.real() - x.imag() * y.imag(),
      x.real() * y.imag() + x.imag() * y.real());
}

/// \brief Binary * operator for RealType times complex.
///
/// This function exists because the compiler doesn't know that
/// RealType and complex<RealType> commute with respect to operator*.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator*(const RealType1& x, const complex<RealType2>& y) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x * y.real(), x * y.imag());
}

/// \brief Binary * operator for RealType times complex.
///
/// This function exists because the compiler doesn't know that
/// RealType and complex<RealType> commute with respect to operator*.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator*(const complex<RealType1>& y, const RealType2& x) noexcept {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      x * y.real(), x * y.imag());
}

//! Imaginary part of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION RealType imag(const complex<RealType>& x) noexcept {
  return x.imag();
}

//! Real part of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION RealType real(const complex<RealType>& x) noexcept {
  return x.real();
}

//! Absolute value (magnitude) of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION RealType abs(const complex<RealType>& x) {
#if !defined(__CUDA_ARCH__) && \
    !defined(__HIP_DEVICE_COMPILE__)  // FIXME_CUDA FIXME_HIP
  using std::hypot;
#endif
  return hypot(x.real(), x.imag());
}

//! Power of a complex number
template <class RealType>
KOKKOS_INLINE_FUNCTION Kokkos::complex<RealType> pow(const complex<RealType>& x,
                                                     const RealType& e) {
  RealType r = abs(x);
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
  using std::atan;
  using std::cos;
  using std::pow;
  using std::sin;
#endif
  using ::pow;
  RealType phi = atan(x.imag() / x.real());
  return pow(r, e) * Kokkos::complex<RealType>(cos(phi * e), sin(phi * e));
}

//! Square root of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION Kokkos::complex<RealType> sqrt(
    const complex<RealType>& x) {
  RealType r = abs(x);
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
  using std::atan;
  using std::cos;
  using std::sin;
  using std::sqrt;
#endif
  using ::sqrt;
  RealType phi = atan(x.imag() / x.real());
  return sqrt(r) * Kokkos::complex<RealType>(cos(phi * 0.5), sin(phi * 0.5));
}

//! Conjugate of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION complex<RealType> conj(
    const complex<RealType>& x) noexcept {
  return complex<RealType>(real(x), -imag(x));
}

//! Exponential of a complex number.
template <class RealType>
KOKKOS_INLINE_FUNCTION complex<RealType> exp(const complex<RealType>& x) {
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
  using std::cos;
  using std::exp;
  using std::sin;
#else
  using ::exp;
#endif
  return exp(x.real()) * complex<RealType>(cos(x.imag()), sin(x.imag()));
}

/// This function cannot be called in a CUDA device function,
/// because std::complex's methods and nonmember functions are not
/// marked as CUDA device functions.
template <class RealType>
inline complex<RealType> exp(const std::complex<RealType>& c) {
  return complex<RealType>(std::exp(c.real()) * std::cos(c.imag()),
                           std::exp(c.real()) * std::sin(c.imag()));
}

//! Binary operator / for complex and real numbers
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator/(const complex<RealType1>& x,
              const RealType2& y) noexcept(noexcept(RealType1{} /
                                                    RealType2{})) {
  return complex<typename std::common_type<RealType1, RealType2>::type>(
      real(x) / y, imag(x) / y);
}

//! Binary operator / for complex.
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator/(const complex<RealType1>& x,
              const complex<RealType2>& y) noexcept(noexcept(RealType1{} /
                                                             RealType2{})) {
  // Scale (by the "1-norm" of y) to avoid unwarranted overflow.
  // If the real part is +/-Inf and the imaginary part is -/+Inf,
  // this won't change the result.
#if !defined(__HIP_DEVICE_COMPILE__)  // FIXME_HIP
  using std::fabs;
#endif
  typedef
      typename std::common_type<RealType1, RealType2>::type common_real_type;
  const common_real_type s = fabs(real(y)) + fabs(imag(y));

  // If s is 0, then y is zero, so x/y == real(x)/0 + i*imag(x)/0.
  // In that case, the relation x/y == (x/s) / (y/s) doesn't hold,
  // because y/s is NaN.
  if (s == 0.0) {
    return complex<common_real_type>(real(x) / s, imag(x) / s);
  } else {
    const complex<common_real_type> x_scaled(real(x) / s, imag(x) / s);
    const complex<common_real_type> y_conj_scaled(real(y) / s, -imag(y) / s);
    const RealType1 y_scaled_abs =
        real(y_conj_scaled) * real(y_conj_scaled) +
        imag(y_conj_scaled) * imag(y_conj_scaled);  // abs(y) == abs(conj(y))
    complex<common_real_type> result = x_scaled * y_conj_scaled;
    result /= y_scaled_abs;
    return result;
  }
}

//! Binary operator / for complex and real numbers
template <class RealType1, class RealType2>
KOKKOS_INLINE_FUNCTION
    complex<typename std::common_type<RealType1, RealType2>::type>
    operator/(const RealType1& x,
              const complex<RealType2>& y) noexcept(noexcept(RealType1{} /
                                                             RealType2{})) {
  return complex<typename std::common_type<RealType1, RealType2>::type>(x) / y;
}

template <class RealType>
std::ostream& operator<<(std::ostream& os, const complex<RealType>& x) {
  const std::complex<RealType> x_std(Kokkos::real(x), Kokkos::imag(x));
  os << x_std;
  return os;
}

template <class RealType>
std::istream& operator>>(std::istream& is, complex<RealType>& x) {
  std::complex<RealType> x_std;
  is >> x_std;
  x = x_std;  // only assigns on success of above
  return is;
}

template <class T>
struct reduction_identity<Kokkos::complex<T> > {
  typedef reduction_identity<T> t_red_ident;
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::complex<T>
  sum() noexcept {
    return Kokkos::complex<T>(t_red_ident::sum(), t_red_ident::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Kokkos::complex<T>
  prod() noexcept {
    return Kokkos::complex<T>(t_red_ident::prod(), t_red_ident::sum());
  }
};

}  // namespace Kokkos

#endif  // KOKKOS_COMPLEX_HPP
