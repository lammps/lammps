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

#ifndef KOKKOS_ANALYZESHAPE_HPP
#define KOKKOS_ANALYZESHAPE_HPP

#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

/** \brief  Analyze the array shape defined by a Kokkos::View data type.
 *
 *  It is presumed that the data type can be mapped down to a multidimensional
 *  array of an intrinsic scalar numerical type (double, float, int, ... ).
 *  The 'value_type' of an array may be an embedded aggregate type such
 *  as a fixed length array 'Array<T,N>'.
 *  In this case the 'array_intrinsic_type' represents the
 *  underlying array of intrinsic scalar numerical type.
 *
 *  The embedded aggregate type must have an AnalyzeShape specialization
 *  to map it down to a shape and intrinsic scalar numerical type.
 */
template< class T >
struct AnalyzeShape : public Shape< sizeof(T) , 0 >
{
  typedef void specialize ;

  typedef Shape< sizeof(T), 0 >  shape ;

  typedef       T  array_intrinsic_type ;
  typedef       T  value_type ;
  typedef       T  type ;

  typedef const T  const_array_intrinsic_type ;
  typedef const T  const_value_type ;
  typedef const T  const_type ;

  typedef       T  non_const_array_intrinsic_type ;
  typedef       T  non_const_value_type ;
  typedef       T  non_const_type ;
};

template<>
struct AnalyzeShape<void> : public Shape< 0 , 0 >
{
  typedef void specialize ;

  typedef Shape< 0 , 0 >  shape ;

  typedef       void  array_intrinsic_type ;
  typedef       void  value_type ;
  typedef       void  type ;
  typedef const void  const_array_intrinsic_type ;
  typedef const void  const_value_type ;
  typedef const void  const_type ;
  typedef       void  non_const_array_intrinsic_type ;
  typedef       void  non_const_value_type ;
  typedef       void  non_const_type ;
};

template< class T >
struct AnalyzeShape< const T > : public AnalyzeShape<T>::shape
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename nested::shape shape ;

  typedef typename nested::const_array_intrinsic_type  array_intrinsic_type ;
  typedef typename nested::const_value_type            value_type ;
  typedef typename nested::const_type                  type ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type ;
};

template< class T >
struct AnalyzeShape< T * >
  : public ShapeInsert< typename AnalyzeShape<T>::shape , 0 >::type
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type * array_intrinsic_type ;
  typedef typename nested::value_type             value_type ;
  typedef typename nested::type                 * type ;

  typedef typename nested::const_array_intrinsic_type * const_array_intrinsic_type ;
  typedef typename nested::const_value_type             const_value_type ;
  typedef typename nested::const_type                 * const_type ;

  typedef typename nested::non_const_array_intrinsic_type * non_const_array_intrinsic_type ;
  typedef typename nested::non_const_value_type             non_const_value_type ;
  typedef typename nested::non_const_type                 * non_const_type ;
};

template< class T >
struct AnalyzeShape< T[] >
  : public ShapeInsert< typename AnalyzeShape<T>::shape , 0 >::type
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [] ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type [] ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type [] ;
};

template< class T >
struct AnalyzeShape< const T[] >
  : public ShapeInsert< typename AnalyzeShape< const T >::shape , 0 >::type
{
private:
  typedef AnalyzeShape< const T > nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [] ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type [] ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type [] ;
};

template< class T , unsigned N >
struct AnalyzeShape< T[N] >
  : public ShapeInsert< typename AnalyzeShape<T>::shape , N >::type
{
private:
  typedef AnalyzeShape<T> nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [N] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [N] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [N] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [N] ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type [N] ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type [N] ;
};

template< class T , unsigned N >
struct AnalyzeShape< const T[N] >
  : public ShapeInsert< typename AnalyzeShape< const T >::shape , N >::type
{
private:
  typedef AnalyzeShape< const T > nested ;
public:

  typedef typename nested::specialize specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_intrinsic_type  array_intrinsic_type [N] ;
  typedef typename nested::value_type            value_type ;
  typedef typename nested::type                  type [N] ;

  typedef typename nested::const_array_intrinsic_type  const_array_intrinsic_type [N] ;
  typedef typename nested::const_value_type            const_value_type ;
  typedef typename nested::const_type                  const_type [N] ;

  typedef typename nested::non_const_array_intrinsic_type  non_const_array_intrinsic_type [N] ;
  typedef typename nested::non_const_value_type            non_const_value_type ;
  typedef typename nested::non_const_type                  non_const_type [N] ;
};

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_ANALYZESHAPE_HPP */

