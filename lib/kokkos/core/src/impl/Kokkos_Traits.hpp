/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

#ifndef KOKKOSTRAITS_HPP
#define KOKKOSTRAITS_HPP

#include <stddef.h>
#include <Kokkos_Macros.hpp>
#include <stdint.h>

namespace Kokkos {
namespace Impl {

/* C++11 conformal compile-time type traits utilities.
 * Prefer to use C++11 when portably available.
 */
//----------------------------------------------------------------------------
// C++11 Helpers:

template < class T , T v >
struct integral_constant
{
  // Declaration of 'static const' causes an unresolved linker symbol in debug
  // static const T value = v ;
  enum { value = T(v) };
  typedef T value_type;
  typedef integral_constant<T,v> type;
  KOKKOS_INLINE_FUNCTION operator T() { return v ; }
};

typedef integral_constant<bool,false> false_type ;
typedef integral_constant<bool,true>  true_type ;

//----------------------------------------------------------------------------
// C++11 Type relationships:

template< class X , class Y > struct is_same : public false_type {};
template< class X >           struct is_same<X,X> : public true_type {};

//----------------------------------------------------------------------------
// C++11 Type properties:

template <typename T> struct is_const : public false_type {};
template <typename T> struct is_const<const T> : public true_type {};
template <typename T> struct is_const<const T & > : public true_type {};

//----------------------------------------------------------------------------
// C++11 Type transformations:

template <typename T> struct remove_const { typedef T type; };
template <typename T> struct remove_const<const T> { typedef T type; };
template <typename T> struct remove_const<const T & > { typedef T & type; };

template <typename T> struct add_const { typedef const T type; };
template <typename T> struct add_const<T & > { typedef const T & type; };
template <typename T> struct add_const<const T> { typedef const T type; };
template <typename T> struct add_const<const T & > { typedef const T & type; };

template<typename T> struct remove_reference { typedef T type ; };
template<typename T> struct remove_reference< T & > { typedef T type ; };
template<typename T> struct remove_reference< const T & > { typedef const T type ; };

//----------------------------------------------------------------------------
// C++11 Other type generators:

template< bool , class T , class F >
struct condition { typedef F type ; };

template< class T , class F >
struct condition<true,T,F> { typedef T type ; };

template< bool , class = void >
struct enable_if ;

template< class T >
struct enable_if< true , T > { typedef T type ; };

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Other traits

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< class , class T = void >
struct enable_if_type { typedef T type ; };

//----------------------------------------------------------------------------

template< bool B >
struct bool_ : public integral_constant<bool,B> {};

template< unsigned I >
struct unsigned_ : public integral_constant<unsigned,I> {};

template< int I >
struct int_ : public integral_constant<int,I> {};

typedef bool_<true> true_;
typedef bool_<false> false_;
//----------------------------------------------------------------------------
// if_

template < bool Cond , typename TrueType , typename FalseType>
struct if_c
{
  enum { value = Cond };

  typedef FalseType type;


  typedef typename remove_const<
          typename remove_reference<type>::type >::type value_type ;

  typedef typename add_const<value_type>::type const_value_type ;

  static KOKKOS_INLINE_FUNCTION
  const_value_type & select( const_value_type & v ) { return v ; }

  static KOKKOS_INLINE_FUNCTION
  value_type & select( value_type & v ) { return v ; }

  template< class T >
  static KOKKOS_INLINE_FUNCTION
  value_type & select( const T & ) { value_type * ptr(0); return *ptr ; }


  template< class T >
  static KOKKOS_INLINE_FUNCTION
  const_value_type & select( const T & , const_value_type & v ) { return v ; }

  template< class T >
  static KOKKOS_INLINE_FUNCTION
  value_type & select( const T & , value_type & v ) { return v ; }
};

template <typename TrueType, typename FalseType>
struct if_c< true , TrueType , FalseType >
{
  enum { value = true };

  typedef TrueType type;


  typedef typename remove_const<
          typename remove_reference<type>::type >::type value_type ;

  typedef typename add_const<value_type>::type const_value_type ;

  static KOKKOS_INLINE_FUNCTION
  const_value_type & select( const_value_type & v ) { return v ; }

  static KOKKOS_INLINE_FUNCTION
  value_type & select( value_type & v ) { return v ; }

  template< class T >
  static KOKKOS_INLINE_FUNCTION
  value_type & select( const T & ) { value_type * ptr(0); return *ptr ; }


  template< class F >
  static KOKKOS_INLINE_FUNCTION
  const_value_type & select( const_value_type & v , const F & ) { return v ; }

  template< class F >
  static KOKKOS_INLINE_FUNCTION
  value_type & select( value_type & v , const F & ) { return v ; }
};

template< typename TrueType >
struct if_c< false , TrueType , void >
{
  enum { value = false };

  typedef void type ;
  typedef void value_type ;
};

template< typename FalseType >
struct if_c< true , void , FalseType >
{
  enum { value = true };

  typedef void type ;
  typedef void value_type ;
};

template <typename Cond, typename TrueType, typename FalseType>
struct if_ : public if_c<Cond::value, TrueType, FalseType> {};

//----------------------------------------------------------------------------

template < size_t N >
struct is_power_of_two
{
  enum type { value = (N > 0) && !(N & (N-1)) };
};

template < size_t N , bool OK = is_power_of_two<N>::value >
struct power_of_two ;

template < size_t N >
struct power_of_two<N,true>
{
  enum type { value = 1+ power_of_two<(N>>1),true>::value };
};

template <>
struct power_of_two<2,true>
{
  enum type { value = 1 };
};

template <>
struct power_of_two<1,true>
{
  enum type { value = 0 };
};

/** \brief  If power of two then return power,
 *          otherwise return ~0u.
 */
static KOKKOS_FORCEINLINE_FUNCTION
unsigned power_of_two_if_valid( const unsigned N )
{
  unsigned p = ~0u ;
  if ( N && ! ( N & ( N - 1 ) ) ) {
#if defined( __CUDA_ARCH__ )
    p = __ffs(N) - 1 ;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
    p = __builtin_ffs(N) - 1 ;
#elif defined( __INTEL_COMPILER )
    p = _bit_scan_forward(N);
#else
    p = 0 ;
    for ( unsigned j = 1 ; ! ( N & j ) ; j <<= 1 ) { ++p ; }
#endif
  }
  return p ;
}

//----------------------------------------------------------------------------

template< typename T , T v , bool NonZero = ( v != T(0) ) >
struct integral_nonzero_constant
{
  // Declaration of 'static const' causes an unresolved linker symbol in debug
  // static const T value = v ;
  enum { value = T(v) };
  typedef T value_type ;
  typedef integral_nonzero_constant<T,v> type ;
  KOKKOS_INLINE_FUNCTION integral_nonzero_constant( const T & ) {}
};

template< typename T , T zero >
struct integral_nonzero_constant<T,zero,false>
{
  const T value ;
  typedef T value_type ;
  typedef integral_nonzero_constant<T,0> type ;
  KOKKOS_INLINE_FUNCTION integral_nonzero_constant( const T & v ) : value(v) {}
};

//----------------------------------------------------------------------------

template <typename T> struct is_integral : public false_ {};

template <> struct is_integral<int8_t>  : public true_ {};
template <> struct is_integral<int16_t> : public true_ {};
template <> struct is_integral<int32_t> : public true_ {};
template <> struct is_integral<int64_t> : public true_ {};

template <> struct is_integral<uint8_t>  : public true_ {};
template <> struct is_integral<uint16_t> : public true_ {};
template <> struct is_integral<uint32_t> : public true_ {};
template <> struct is_integral<uint64_t> : public true_ {};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSTRAITS_HPP */

