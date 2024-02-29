//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_HALF_NUMERIC_TRAITS_HPP_
#define KOKKOS_HALF_NUMERIC_TRAITS_HPP_

#include <Kokkos_NumericTraits.hpp>

////////////// BEGIN HALF_T (binary16) limits //////////////
// clang-format off
// '\brief:' below are from the libc definitions for float and double:
// https://www.gnu.org/software/libc/manual/html_node/Floating-Point-Parameters.html
//
// The arithmetic encoding and equations below are derived from:
// Ref1: https://en.wikipedia.org/wiki/Single-precision_floating-point_format
// Ref2: https://en.wikipedia.org/wiki/Exponent_bias
// Ref3; https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
//
// Some background on the magic numbers 2**10=1024 and 2**15=32768 used below:
//
// IMPORTANT: For IEEE754 encodings, see Ref1.
//
// For binary16, we have B = 2 and p = 16 with 2**16 possible significands.
// The binary16 format is: [s  e  e  e  e  e  f f f f f f f f f f]
//              bit index:  15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
// s: signed bit (1 bit)
// e: exponent bits (5 bits)
// f: fractional bits (10 bits)
//
// E_bias      = 2**(n_exponent_bits - 1) - 1 = 2**(5 - 1) - 1 = 15
// E_subnormal = 00000 (base2)
// E_infinity  = 11111 (base2)
// E_min       = 1 - E_bias = 1 - 15
// E_max       = 2**5 - 1 - E_bias = 2**5 - 1 - 15 = 16
//
// 2**10=1024 is the smallest denominator that is representable in binary16:
// [s  e  e  e  e  e  f f f f f f f f f f]
// [0  0  0  0  0  0  0 0 0 0 0 0 0 0 0 1]
// which is: 1 / 2**-10
//
//
// 2**15 is the largest exponent factor representable in binary16, for example the
// largest integer value representable in binary16 is:
// [s  e  e  e  e  e  f f f f f f f f f f]
// [0  1  1  1  1  0  1 1 1 1 1 1 1 1 1 1]
// which is: 2**(2**4 + 2**3 + 2**2 + 2**1 - 15) * (1 + 2**-10 + 2**-9 + 2**-8 + 2**-7 + 2**-6 + 2**-5 + 2**-4 + 2**-3 + 2**-2 + 2**-1)) =
//           2**15 * (1 + 0.9990234375) =
//           65504.0
//
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
/// \brief: Infinity
///
/// Binary16 encoding:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  1  1  1  1  1  0 0 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
template <>
struct Kokkos::Experimental::Impl::infinity_helper<Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'11111'0000000000};
};

/// \brief: Minimum normalized number
///
/// Stdc defines this as the smallest number (representable in binary16).
///
/// Binary16 encoding:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [1  1  1  1  1  0  1 1 1 1 1 1 1 1 1 1]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
/// and in base10: -1 * 2**(2**4 + 2**3 + 2**2 + 2**1 - 15) * (1 + 2**-10 + 2**-9 + 2**-8 + 2**-7 + 2**-6 + 2**-5 + 2**-4 + 2**-3 + 2**-2 + 2**-1)
///              = -2**15 * (1 + (2**10 - 1) / 2**10)
template <>
struct Kokkos::Experimental::Impl::finite_min_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b1'11110'1111111111}; // -65504
};

/// \brief: Maximum normalized number
///
/// Stdc defines this as the maximum number (representable in binary16).
///
/// Binary16 encoding:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  1  1  1  1  0  1 1 1 1 1 1 1 1 1 1]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
/// and in base10: 1 * 2**(2**4 + 2**3 + 2**2 + 2**1 - 15) * (1 + 2**-10 + 2**-9 + 2**-8 + 2**-7 + 2**-6 + 2**-5 + 2**-4 + 2**-3 + 2**-2 + 2**-1)
///              = 2**15 * (1 + (2**10 - 1) / 2**10)
template <>
struct Kokkos::Experimental::Impl::finite_max_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'11110'1111111111}; // +65504
};

/// \brief: This is the difference between 1 and the smallest floating point
///         number of type binary16 that is greater than 1
///
/// Smallest number in binary16 that is greater than 1 encoding:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  0  1  1  1  1  0 0 0 0 0 0 0 0 0 1]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
/// and in base10: 1 * 2**(2**3 + 2**2 + 2**1 + 2**0 - 15) * (1 + 2**-10)
///                = 2**0 * (1 + 2**-10)
///                = 1.0009765625
///
/// Lastly, 1 - 1.0009765625 = 0.0009765625.
template <>
struct Kokkos::Experimental::Impl::epsilon_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'00101'0000000000}; // 0.0009765625
};

/// @brief: The largest possible rounding error in ULPs
///
/// This simply uses the maximum rounding error.
///
/// Reference: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html#689
template <>
struct Kokkos::Experimental::Impl::round_error_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'01110'0000000000}; // 0.5
};

/// \brief: Minimum normalized positive half precision number
///
/// Stdc defines this as the minimum normalized positive floating
/// point number that is representable in type binary16
///
/// Smallest number in binary16 that is greater than 1 encoding:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  0  0  0  0  1  0 0 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
/// and in base10: 1 * 2**(2**0 - 15) * (1)
///                = 2**-14
template <>
struct Kokkos::Experimental::Impl::norm_min_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'00001'0000000000}; // 0.00006103515625
};

/// \brief: Quiet not a half precision number
///
/// IEEE 754 defines this as all exponent bits and the first fraction bit high.
///
/// Quiet NaN in binary16:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  1  1  1  1  1  1 0 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
template <>
struct Kokkos::Experimental::Impl::quiet_NaN_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'11111'1000000000};
};

/// \brief: Signaling not a half precision number
///
/// IEEE 754 defines this as all exponent bits and the second fraction bit high.
///
/// Quiet NaN in binary16:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  1  1  1  1  1  0 1 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
template <>
struct Kokkos::Experimental::Impl::signaling_NaN_helper<
    Kokkos::Experimental::half_t> {
  static constexpr Kokkos::Experimental::half_t::bit_comparison_type value{0b0'11111'0100000000};
};

/// \brief: Number of digits in the matissa that can be represented
///         without losing precision.
///
/// Stdc defines this as the number of base-RADIX digits in the floating point mantissa for the binary16 data type.
///
/// In binary16, we have 10 fractional bits plus the implicit leading 1.
template <>
struct Kokkos::Experimental::Impl::digits_helper<Kokkos::Experimental::half_t> {
  static constexpr int value = 11;
};

/// \brief: "The number of base-10 digits that can be represented by the type T without change"
/// Reference: https://en.cppreference.com/w/cpp/types/numeric_limits/digits10.
///
/// "For base-radix types, it is the value of digits() (digits - 1 for floating-point types) multiplied by log10(radix) and rounded down."
/// Reference: https://en.cppreference.com/w/cpp/types/numeric_limits/digits10.
///
/// This is: floor(11 - 1 * log10(2))
template <>
struct Kokkos::Experimental::Impl::digits10_helper<
    Kokkos::Experimental::half_t> {
  static constexpr int value = 3;
};

/// \brief: Value of the base of the exponent representation.
///
/// Stdc defined this as the value of the base, or radix, of the exponent representation.
template <>
struct Kokkos::Experimental::Impl::radix_helper<Kokkos::Experimental::half_t> {
  static constexpr int value = 2;
};

/// \brief: This is the smallest possible exponent value
///
/// Stdc defines this as the smallest possible exponent value for type binary16. 
/// More precisely, it is the minimum negative integer such that the value min_exponent_helper
/// raised to this power minus 1 can be represented as a normalized floating point number of type float.
///
/// In binary16:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  0  0  0  0  1  0 0 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
/// 
/// and in base10: 1 * 2**(2**0 - 15) * (1 + 0)
///                = 2**-14
/// 
/// with a bias of one from (C11 5.2.4.2.2), gives -13;
template <>
struct Kokkos::Experimental::Impl::min_exponent_helper<
    Kokkos::Experimental::half_t> {
  static constexpr int value = -13;
};

/// \brief: This is the largest possible exponent value
///
/// In binary16:
///             [s  e  e  e  e  e  f f f f f f f f f f]
///             [0  1  1  1  1  0  0 0 0 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
/// 
/// and in base10: 1 * 2**(2**4 + 2**3 + 2**2 + 2**1 - 15) * (1 + 0)
///                = 2**(30 - 15)
///                = 2**15
/// 
/// with a bias of one from (C11 5.2.4.2.2), gives 16;
template <>
struct Kokkos::Experimental::Impl::max_exponent_helper<
    Kokkos::Experimental::half_t> {
  static constexpr int value = 16;
};
#endif
////////////// END HALF_T (binary16) limits //////////////

////////////// BEGIN BHALF_T (bfloat16) limits //////////////
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
/// \brief: Infinity
///
/// Bfloat16 encoding:
///             [s  e  e  e  e  e  e e e f f f f f f f]
///             [0  1  1  1  1  1  1 1 1 0 0 0 0 0 0 0]
/// bit index:   15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
///
template <>
struct Kokkos::Experimental::Impl::infinity_helper<Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'11111111'0000000};
};

// Minimum normalized number
template <>
struct Kokkos::Experimental::Impl::finite_min_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b1'11111110'1111111}; // -3.38953139e38
};
// Maximum normalized number
template <>
struct Kokkos::Experimental::Impl::finite_max_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'11111110'1111111}; // +3.38953139e3
};
// 1/2^7
template <>
struct Kokkos::Experimental::Impl::epsilon_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'01111000'0000000}; // 0.0078125
};
template <>
struct Kokkos::Experimental::Impl::round_error_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'01111110'0000000}; // 0.5
};
// Minimum normalized positive bhalf number
template <>
struct Kokkos::Experimental::Impl::norm_min_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'00000001'0000000}; // 1.175494351e-38
};
// Quiet not a bhalf number
template <>
struct Kokkos::Experimental::Impl::quiet_NaN_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'11111111'1000000};
};
// Signaling not a bhalf number
template <>
struct Kokkos::Experimental::Impl::signaling_NaN_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr Kokkos::Experimental::bhalf_t::bit_comparison_type value{0b0'11111111'0100000};
};
// Number of digits in the matissa that can be represented
// without losing precision.
template <>
struct Kokkos::Experimental::Impl::digits_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr int value = 2;
};
// 7 - 1 * log10(2)
template <>
struct Kokkos::Experimental::Impl::digits10_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr int value = 1;
};
// Value of the base of the exponent representation.
template <>
struct Kokkos::Experimental::Impl::radix_helper<Kokkos::Experimental::bhalf_t> {
  static constexpr int value = 2;
};
// This is the smallest possible exponent value
// with a bias of one (C11 5.2.4.2.2).
template <>
struct Kokkos::Experimental::Impl::min_exponent_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr int value = -125;
};
// This is the largest possible exponent value
// with a bias of one (C11 5.2.4.2.2).
template <>
struct Kokkos::Experimental::Impl::max_exponent_helper<
    Kokkos::Experimental::bhalf_t> {
  static constexpr int value = 128;
};
#endif
////////////// END BHALF_T (bfloat16) limits //////////

#endif  // KOKKOS_HALF_NUMERIC_TRAITS_HPP_
