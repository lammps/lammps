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

#ifndef KOKKOS_ARRAY_HPP
#define KOKKOS_ARRAY_HPP

#include <type_traits>
#include <algorithm>
#include <limits>
#include <cstddef>

namespace Kokkos {

/**\brief  Derived from the C++17 'std::array'.
 *         Dropping the iterator interface.
 */
template< class T      = void
        , size_t N     = ~size_t(0)
        , class Proxy  = void
        >
struct Array {
public:
  /**
   * The elements of this C array shall not be accessed directly. The data
   * member has to be declared public to enable aggregate initialization as for
   * std::array. We mark it as private in the documentation.
   * @private
   */
  T m_internal_implementation_private_member_data[N];
public:

  typedef T &                                 reference ;
  typedef typename std::add_const<T>::type &  const_reference ;
  typedef size_t                              size_type ;
  typedef ptrdiff_t                           difference_type ;
  typedef T                                   value_type ;
  typedef T *                                 pointer ;
  typedef typename std::add_const<T>::type *  const_pointer ;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size() { return N ; }
  KOKKOS_INLINE_FUNCTION static constexpr bool      empty(){ return false ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference operator[]( const iType & i )
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_internal_implementation_private_member_data[i];
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  const_reference operator[]( const iType & i ) const
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_internal_implementation_private_member_data[i];
    }

  KOKKOS_INLINE_FUNCTION pointer       data()
    {
      return & m_internal_implementation_private_member_data[0];
    }
  KOKKOS_INLINE_FUNCTION const_pointer data() const
    {
      return & m_internal_implementation_private_member_data[0];
    }

  // Do not default unless move and move-assignment are also defined
  // ~Array() = default ;
  // Array() = default ;
  // Array( const Array & ) = default ;
  // Array & operator = ( const Array & ) = default ;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && ) = default ;
  // Array & operator = ( Array && ) = default ;
};


template< class T , class Proxy >
struct Array<T,0,Proxy> {
public:

  typedef typename std::add_const<T>::type &  reference ;
  typedef typename std::add_const<T>::type &  const_reference ;
  typedef size_t                              size_type ;
  typedef ptrdiff_t                           difference_type ;
  typedef typename std::add_const<T>::type    value_type ;
  typedef typename std::add_const<T>::type *  pointer ;
  typedef typename std::add_const<T>::type *  const_pointer ;

  KOKKOS_INLINE_FUNCTION static constexpr size_type size()  { return 0 ; }
  KOKKOS_INLINE_FUNCTION static constexpr bool      empty() { return true ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  value_type operator[]( const iType & )
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integer argument" );
      return value_type();
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  value_type operator[]( const iType & ) const
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integer argument" );
      return value_type();
    }

  KOKKOS_INLINE_FUNCTION pointer       data()       { return pointer(0) ; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return const_pointer(0); }

  KOKKOS_FUNCTION_DEFAULTED ~Array() = default ;
  KOKKOS_FUNCTION_DEFAULTED Array() = default ;
  KOKKOS_FUNCTION_DEFAULTED Array( const Array & ) = default ;
  KOKKOS_FUNCTION_DEFAULTED Array & operator = ( const Array & ) = default ;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && ) = default ;
  // Array & operator = ( Array && ) = default ;
};


template<>
struct Array<void,~size_t(0),void>
{
  struct contiguous {};
  struct strided {};
};

template< class T >
struct Array< T , ~size_t(0) , Array<>::contiguous >
{
private:
  T *    m_elem ;
  size_t m_size ;
public:

  typedef T &                                 reference ;
  typedef typename std::add_const<T>::type &  const_reference ;
  typedef size_t                              size_type ;
  typedef ptrdiff_t                           difference_type ;
  typedef T                                   value_type ;
  typedef T *                                 pointer ;
  typedef typename std::add_const<T>::type *  const_pointer ;

  KOKKOS_INLINE_FUNCTION constexpr size_type size()  const { return m_size ; }
  KOKKOS_INLINE_FUNCTION constexpr bool      empty() const { return 0 != m_size ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference operator[]( const iType & i )
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_elem[i];
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  const_reference operator[]( const iType & i ) const
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_elem[i];
    }

  KOKKOS_INLINE_FUNCTION pointer       data()       { return m_elem ; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return m_elem ; }

  KOKKOS_FUNCTION_DEFAULTED ~Array() = default ;
  Array() = delete ;
  Array( const Array & rhs ) = delete ;

  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && rhs ) = default ;
  // Array & operator = ( Array && rhs ) = delete ;

  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
    {
      const size_t n = std::min( m_size , rhs.size() );
      for ( size_t i = 0 ; i < n ; ++i ) m_elem[i] = rhs[i] ;
      return *this ;
    }

  template< size_t N , class P >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,P> & rhs )
    {
      const size_t n = std::min( m_size , rhs.size() );
      for ( size_t i = 0 ; i < n ; ++i ) m_elem[i] = rhs[i] ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION constexpr Array( pointer arg_ptr , size_type arg_size , size_type = 0 )
    : m_elem(arg_ptr), m_size(arg_size) {}
};

template< class T >
struct Array< T , ~size_t(0) , Array<>::strided >
{
private:
  T *    m_elem ;
  size_t m_size ;
  size_t m_stride ;
public:

  typedef T &                                 reference ;
  typedef typename std::add_const<T>::type &  const_reference ;
  typedef size_t                              size_type ;
  typedef ptrdiff_t                           difference_type ;
  typedef T                                   value_type ;
  typedef T *                                 pointer ;
  typedef typename std::add_const<T>::type *  const_pointer ;

  KOKKOS_INLINE_FUNCTION constexpr size_type size()  const { return m_size ; }
  KOKKOS_INLINE_FUNCTION constexpr bool      empty() const { return 0 != m_size ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  reference operator[]( const iType & i )
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_elem[i*m_stride];
    }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  const_reference operator[]( const iType & i ) const
    {
      static_assert( ( std::is_integral<iType>::value || std::is_enum<iType>::value ) , "Must be integral argument" );
      return m_elem[i*m_stride];
    }

  KOKKOS_INLINE_FUNCTION pointer       data()       { return m_elem ; }
  KOKKOS_INLINE_FUNCTION const_pointer data() const { return m_elem ; }

  KOKKOS_FUNCTION_DEFAULTED ~Array() = default ;
  Array()  = delete ;
  Array( const Array & ) = delete ;


  // Some supported compilers are not sufficiently C++11 compliant
  // for default move constructor and move assignment operator.
  // Array( Array && rhs ) = default ;
  // Array & operator = ( Array && rhs ) = delete ;

  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array & rhs )
    {
      const size_t n = std::min( m_size , rhs.size() );
      for ( size_t i = 0 ; i < n ; ++i ) m_elem[i] = rhs[i] ;
      return *this ;
    }

  template< size_t N , class P >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,P> & rhs )
    {
      const size_t n = std::min( m_size , rhs.size() );
      for ( size_t i = 0 ; i < n ; ++i ) m_elem[i] = rhs[i] ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION constexpr Array( pointer arg_ptr , size_type arg_size , size_type arg_stride )
    : m_elem(arg_ptr), m_size(arg_size), m_stride(arg_stride) {}
};

} // namespace Kokkos

#endif /* #ifndef KOKKOS_ARRAY_HPP */

