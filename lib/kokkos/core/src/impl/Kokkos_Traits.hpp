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

#ifndef KOKKOSTRAITS_HPP
#define KOKKOSTRAITS_HPP

#include <cstddef>
#include <cstdint>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_BitOps.hpp>
#include <string>
#include <type_traits>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// Help with C++11 variadic argument packs

template< unsigned I , typename ... Pack >
struct get_type { typedef void type ; };

template< typename T , typename ... Pack >
struct get_type< 0 , T , Pack ... >
{ typedef T type ; };

template< unsigned I , typename T , typename ... Pack >
struct get_type< I , T , Pack ... >
{ typedef typename get_type< I - 1 , Pack ... >::type type ; };


template< typename T , typename ... Pack >
struct has_type { enum { value = false }; };

template< typename T , typename S , typename ... Pack >
struct has_type<T,S,Pack...>
{
private:

  enum { self_value = std::is_same<T,S>::value };

  typedef has_type<T,Pack...> next ;

  static_assert( ! ( self_value && next::value )
               , "Error: more than one member of the argument pack matches the type" );

public:

  enum { value = self_value || next::value };

};


template< typename DefaultType
        , template< typename > class Condition
        , typename ... Pack >
struct has_condition
{
  enum { value = false };
  typedef DefaultType type ;
};

template< typename DefaultType
        , template< typename > class Condition
        , typename S
        , typename ... Pack >
struct has_condition< DefaultType , Condition , S , Pack... >
{
private:

  enum { self_value = Condition<S>::value };

  typedef has_condition< DefaultType , Condition , Pack... > next ;

  static_assert( ! ( self_value && next::value )
               , "Error: more than one member of the argument pack satisfies condition" );

public:

  enum { value = self_value || next::value };

  typedef typename
    std::conditional< self_value , S , typename next::type >::type
      type ;
};


template< class ... Args >
struct are_integral { enum { value = true }; };

template< typename T , class ... Args >
struct are_integral<T,Args...> {
  enum { value =
    // Accept std::is_integral OR std::is_enum as an integral value
    // since a simple enum value is automically convertable to an
    // integral value.
    ( std::is_integral<T>::value || std::is_enum<T>::value )
    &&
    are_integral<Args...>::value };
};

//----------------------------------------------------------------------------
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

template <typename T> struct is_array : public false_type {};
template <typename T> struct is_array< T[] > : public true_type {};
template <typename T, unsigned N > struct is_array< T[N] > : public true_type {};

//----------------------------------------------------------------------------
// C++11 Type transformations:

template <typename T> struct remove_const { typedef T type; };
template <typename T> struct remove_const<const T> { typedef T type; };
template <typename T> struct remove_const<const T & > { typedef T & type; };

template <typename T> struct add_const { typedef const T type; };
template <typename T> struct add_const<T & > { typedef const T & type; };
template <typename T> struct add_const<const T> { typedef const T type; };
template <typename T> struct add_const<const T & > { typedef const T & type; };

template <typename T> struct remove_reference { typedef T type ; };
template <typename T> struct remove_reference< T & > { typedef T type ; };
template <typename T> struct remove_reference< const T & > { typedef const T type ; };

template <typename T> struct remove_extent { typedef T type ; };
template <typename T> struct remove_extent<T[]> { typedef T type ; };
template <typename T, unsigned N > struct remove_extent<T[N]> { typedef T type ; };

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

// Allows aliased types:
template< typename T >
struct is_integral : public integral_constant< bool ,
  (
    std::is_same< T ,          char >::value ||
    std::is_same< T , unsigned char >::value ||
    std::is_same< T ,          short int >::value ||
    std::is_same< T , unsigned short int >::value ||
    std::is_same< T ,          int >::value ||
    std::is_same< T , unsigned int >::value ||
    std::is_same< T ,          long int >::value ||
    std::is_same< T , unsigned long int >::value ||
    std::is_same< T ,          long long int >::value ||
    std::is_same< T , unsigned long long int >::value ||

    std::is_same< T , int8_t   >::value ||
    std::is_same< T , int16_t  >::value ||
    std::is_same< T , int32_t  >::value ||
    std::is_same< T , int64_t  >::value ||
    std::is_same< T , uint8_t  >::value ||
    std::is_same< T , uint16_t >::value ||
    std::is_same< T , uint32_t >::value ||
    std::is_same< T , uint64_t >::value
  )>
{};
//----------------------------------------------------------------------------

template<typename T>
struct is_label : public false_type {};

template<>
struct is_label<const char*> : public true_type {};

template<>
struct is_label<char*> : public true_type {};


template<int N>
struct is_label<const char[N]> : public true_type {};

template<int N>
struct is_label<char[N]> : public true_type {};


template<>
struct is_label<const std::string> : public true_type {};

template<>
struct is_label<std::string> : public true_type {};

// These 'constexpr'functions can be used as
// both regular functions and meta-function.

/**\brief  There exists integral 'k' such that N = 2^k */
KOKKOS_INLINE_FUNCTION
constexpr bool is_integral_power_of_two( const size_t N )
{ return ( 0 < N ) && ( 0 == ( N & ( N - 1 ) ) ); }

/**\brief  Return integral 'k' such that N = 2^k, assuming valid.  */
KOKKOS_INLINE_FUNCTION
constexpr unsigned integral_power_of_two_assume_valid( const size_t N )
{ return N == 1 ? 0 : 1 + integral_power_of_two_assume_valid( N >> 1 ); }

/**\brief  Return integral 'k' such that N = 2^k, if exists.
 *         If does not exist return ~0u.
 */
KOKKOS_INLINE_FUNCTION
constexpr unsigned integral_power_of_two( const size_t N )
{ return is_integral_power_of_two(N) ? integral_power_of_two_assume_valid(N) : ~0u ; }

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
KOKKOS_FORCEINLINE_FUNCTION
unsigned power_of_two_if_valid( const unsigned N )
{
  unsigned p = ~0u ;
  if ( is_integral_power_of_two ( N ) ) {
    p = bit_scan_forward ( N ) ;
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

template < class C > struct is_integral_constant : public false_
{
  typedef void integral_type ;
  enum { integral_value = 0 };
};

template < typename T , T v >
struct is_integral_constant< integral_constant<T,v> > : public true_
{
  typedef T integral_type ;
  enum { integral_value = v };
};

//----------------------------------------------------------------------------

template <class...>
class TypeList;

//----------------------------------------------------------------------------

template <class>
struct ReverseTypeList;

template <class Head, class... Tail>
struct ReverseTypeList<TypeList<Head, Tail...>> {
  template <class... ReversedTail>
  struct impl {
    using type = typename ReverseTypeList<TypeList<Tail...>>::template impl<Head, ReversedTail...>::type;
  };
  using type = typename impl<>::type;
};

template <>
struct ReverseTypeList<TypeList<>> {
  template <class... ReversedTail>
  struct impl {
    using type = TypeList<ReversedTail...>;
  };
  using type = TypeList<>;
};

//----------------------------------------------------------------------------

template <class T>
struct make_all_extents_into_pointers
{
  using type = T;
};

template <class T, unsigned N>
struct make_all_extents_into_pointers<T[N]>
{
  using type = typename make_all_extents_into_pointers<T>::type*;
};

template <class T>
struct make_all_extents_into_pointers<T*>
{
  using type = typename make_all_extents_into_pointers<T>::type*;
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSTRAITS_HPP */

