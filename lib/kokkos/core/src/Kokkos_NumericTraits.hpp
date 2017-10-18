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

#ifndef KOKKOS_NUMERICTRAITS_HPP
#define KOKKOS_NUMERICTRAITS_HPP

#include<climits>
#include<cfloat>

namespace Kokkos {

template<class T>
struct reduction_identity; /*{
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T sum() { return T(); }  // 0
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T prod()  // 1
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom prod reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T max()   // minimum value
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom max reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T min()   // maximum value
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom min reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T bor()   // 0, only for integer type
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom bor reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T band()  // !0, only for integer type
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom band reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T lor()   // 0, only for integer type
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom lor reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T land()  // !0, only for integer type
    { static_assert( false, "Missing specialization of Kokkos::reduction_identity for custom land reduction type"); return T(); }
};*/

template<>
struct reduction_identity<signed char> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char sum()  {return static_cast<signed char>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char prod() {return static_cast<signed char>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char max()  {return SCHAR_MIN;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char min()  {return SCHAR_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char bor()  {return static_cast<signed char>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char band() {return ~static_cast<signed char>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char lor()  {return static_cast<signed char>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char land() {return static_cast<signed char>(1);}
};

template<>
struct reduction_identity<short> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short sum()  {return static_cast<short>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short prod() {return static_cast<short>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short max()  {return SHRT_MIN;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short min()  {return SHRT_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short bor()  {return static_cast<short>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short band() {return ~static_cast<short>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short lor()  {return static_cast<short>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short land() {return static_cast<short>(1);}
};

template<>
struct reduction_identity<int> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int sum()  {return static_cast<int>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int prod() {return static_cast<int>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int max()  {return INT_MIN;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int min()  {return INT_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int bor()  {return static_cast<int>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int band() {return ~static_cast<int>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int lor()  {return static_cast<int>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int land() {return static_cast<int>(1);}
};

template<>
struct reduction_identity<long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long sum()  {return static_cast<long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long prod() {return static_cast<long>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long max()  {return LLONG_MIN;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long min()  {return LLONG_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long bor()  {return static_cast<long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long band() {return ~static_cast<long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long lor()  {return static_cast<long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long land() {return static_cast<long>(1);}
};

template<>
struct reduction_identity<long long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long sum()  {return static_cast<long long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long prod() {return static_cast<long long>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long max()  {return LLONG_MIN;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long min()  {return LLONG_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long bor()  {return static_cast<long long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long band() {return ~static_cast<long long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long lor()  {return static_cast<long long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long land() {return static_cast<long long>(1);}
};

template<>
struct reduction_identity<unsigned char> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char sum()  {return static_cast<unsigned char>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char prod() {return static_cast<unsigned char>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char max()  {return static_cast<unsigned char>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char min()  {return UCHAR_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char bor()  {return static_cast<unsigned char>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char band() {return ~static_cast<unsigned char>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char lor()  {return static_cast<unsigned char>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char land() {return static_cast<unsigned char>(1);}
};

template<>
struct reduction_identity<unsigned short> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short sum()  {return static_cast<unsigned short>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short prod() {return static_cast<unsigned short>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short max()  {return static_cast<unsigned short>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short min()  {return USHRT_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short bor()  {return static_cast<unsigned short>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short band() {return ~static_cast<unsigned short>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short lor()  {return static_cast<unsigned short>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short land() {return static_cast<unsigned short>(1);}
};

template<>
struct reduction_identity<unsigned int> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int sum()  {return static_cast<unsigned int>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int prod() {return static_cast<unsigned int>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int max()  {return static_cast<unsigned int>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int min()  {return UINT_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int bor()  {return static_cast<unsigned int>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int band() {return ~static_cast<unsigned int>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int lor()  {return static_cast<unsigned int>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int land() {return static_cast<unsigned int>(1);}
};

template<>
struct reduction_identity<unsigned long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long sum()  {return static_cast<unsigned long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long prod() {return static_cast<unsigned long>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long max()  {return static_cast<unsigned long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long min()  {return ULONG_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long bor()  {return static_cast<unsigned long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long band() {return ~static_cast<unsigned long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long lor()  {return static_cast<unsigned long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long land() {return static_cast<unsigned long>(1);}
};

template<>
struct reduction_identity<unsigned long long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long sum()  {return static_cast<unsigned long long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long prod() {return static_cast<unsigned long long>(1);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long max()  {return static_cast<unsigned long long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long min()  {return ULLONG_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long bor()  {return static_cast<unsigned long long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long band() {return ~static_cast<unsigned long long>(0x0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long lor()  {return static_cast<unsigned long long>(0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long land() {return static_cast<unsigned long long>(1);}
};

template<>
struct reduction_identity<float> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float sum()  {return static_cast<float>(0.0f);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float prod() {return static_cast<float>(1.0f);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float max()  {return -FLT_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float min()  {return FLT_MAX;}
};

template<>
struct reduction_identity<double> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double sum()  {return static_cast<double>(0.0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double prod() {return static_cast<double>(1.0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double max()  {return -DBL_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double min()  {return DBL_MAX;}
};

template<>
struct reduction_identity<long double> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long double sum()  {return static_cast<long double>(0.0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long double prod() {return static_cast<long double>(1.0);}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long double max()  {return -LDBL_MAX;}
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long double min()  {return LDBL_MAX;}
};

}

#endif
