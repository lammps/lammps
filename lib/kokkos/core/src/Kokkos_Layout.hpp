/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Layout.hpp
/// \brief Declaration of various \c MemoryLayout options.

#ifndef KOKKOS_LAYOUT_HPP
#define KOKKOS_LAYOUT_HPP

#include <stddef.h>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/// \struct LayoutLeft
/// \brief Memory layout tag indicating left-to-right (Fortran scheme)
///   striding of multi-indices.
///
/// This is an example of a \c MemoryLayout template parameter of
/// View.  The memory layout describes how View maps from a
/// multi-index (i0, i1, ..., ik) to a memory location.  
///
/// "Layout left" indicates a mapping where the leftmost index i0
/// refers to contiguous access, and strides increase for dimensions
/// going right from there (i1, i2, ...).  This layout imitates how
/// Fortran stores multi-dimensional arrays.  For the special case of
/// a two-dimensional array, "layout left" is also called "column
/// major."
struct LayoutLeft {
  //! The tag (what type of kokkos_object is this).
  typedef Impl::LayoutTag       kokkos_tag ;
  typedef LayoutLeft array_layout ;
};

//----------------------------------------------------------------------------
/// \struct LayoutRight
/// \brief Memory layout tag indicating right-to-left (C or
///   lexigraphical scheme) striding of multi-indices.
///
/// This is an example of a \c MemoryLayout template parameter of
/// View.  The memory layout describes how View maps from a
/// multi-index (i0, i1, ..., ik) to a memory location.  
///
/// "Right layout" indicates a mapping where the rightmost index ik
/// refers to contiguous access, and strides increase for dimensions
/// going left from there.  This layout imitates how C stores
/// multi-dimensional arrays.  For the special case of a
/// two-dimensional array, "layout right" is also called "row major."
struct LayoutRight {
  //! The tag (what type of kokkos_object is this).
  typedef Impl::LayoutTag       kokkos_tag ;
  typedef LayoutRight array_layout ;
};

//----------------------------------------------------------------------------
/// \struct LayoutStride
/// \brief  Memory layout tag indicated arbitrarily strided
///         multi-index mapping into contiguous memory.
struct LayoutStride {

  //! The tag (what type of kokkos_object is this).
  typedef Impl::LayoutTag       kokkos_tag ;

  typedef LayoutStride array_layout ;

  enum { MAX_RANK = 8 };

  size_t dimension[ MAX_RANK ] ;
  size_t stride[ MAX_RANK ] ; 

  /** \brief  Compute strides from ordered dimensions.
   *
   *  Values of order uniquely form the set [0..rank)
   *  and specify ordering of the dimensions.
   *  Order = {0,1,2,...} is LayoutLeft
   *  Order = {...,2,1,0} is LayoutRight
   */
  template< typename iTypeOrder , typename iTypeDimen >
  KOKKOS_INLINE_FUNCTION static
  LayoutStride order_dimensions( int const rank
                               , iTypeOrder const * const order
                               , iTypeDimen const * const dimen )
    {
      LayoutStride tmp ;
      // Verify valid rank order:
      int check_input = MAX_RANK < rank ? 0 : int( 1 << rank ) - 1 ;
      for ( int r = 0 ; r < MAX_RANK ; ++r ) {
        tmp.dimension[r] = 0 ;
        tmp.stride[r]    = 0 ;
        check_input &= ~int( 1 << order[r] );
      }
      if ( 0 == check_input ) {
        size_t n = 1 ;
        for ( int r = 0 ; r < rank ; ++r ) {
          tmp.stride[ order[r] ] = n ;
          n *= ( tmp.dimension[r] = dimen[r] );
        }
      }
      return tmp ;
    }
};

//----------------------------------------------------------------------------
/// \struct LayoutTileLeft
/// \brief Memory layout tag indicating left-to-right (Fortran scheme)
///   striding of multi-indices by tiles.
///
/// This is an example of a \c MemoryLayout template parameter of
/// View.  The memory layout describes how View maps from a
/// multi-index (i0, i1, ..., ik) to a memory location.  
///
/// "Tiled layout" indicates a mapping to contiguously stored
/// <tt>ArgN0</tt> by <tt>ArgN1</tt> tiles for the rightmost two
/// dimensions.  Indices are LayoutLeft within each tile, and the
/// tiles themselves are arranged using LayoutLeft.  Note that the
/// dimensions <tt>ArgN0</tt> and <tt>ArgN1</tt> of the tiles must be
/// compile-time constants.  This speeds up index calculations.  If
/// both tile dimensions are powers of two, Kokkos can optimize
/// further.
template < unsigned ArgN0 , unsigned ArgN1 ,
           bool IsPowerOfTwo = ( Impl::is_power_of_two<ArgN0>::value &&
                                 Impl::is_power_of_two<ArgN1>::value )
         >
struct LayoutTileLeft {
  //! The tag (what type of kokkos_object is this).
  typedef Impl::LayoutTag       kokkos_tag ;

  typedef LayoutTileLeft<ArgN0,ArgN1,IsPowerOfTwo> array_layout ;
  enum { N0 = ArgN0 };
  enum { N1 = ArgN1 };
};

} // namespace Kokkos

#endif // #ifndef KOKKOS_LAYOUT_HPP

