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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP

#include <type_traits>
#include <initializer_list>

#include <Kokkos_Pair.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Atomic_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ExecPolicy > class ParallelFor ;

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

template< long sN0 = -1
        , long sN1 = -1
        , long sN2 = -1
        , long sN3 = -1
        , long sN4 = -1
        , long sN5 = -1
        , long sN6 = -1
        , long sN7 = -1
        >
struct ViewDimension {

  enum { arg_N0 = sN0 };
  enum { arg_N1 = sN1 };
  enum { arg_N2 = sN2 };
  enum { arg_N3 = sN3 };
  enum { arg_N4 = sN4 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN0 < 0 ? 0 :
                ( sN1 < 0 ? 1 :
                ( sN2 < 0 ? 2 :
                ( sN3 < 0 ? 3 :
                ( sN4 < 0 ? 4 :
                ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 )))))))) };
  enum { rank_dynamic = 0 };

  enum { N0 = 0 < sN0 ? sN0 : 1 };
  enum { N1 = 0 < sN1 ? sN1 : 1 };
  enum { N2 = 0 < sN2 ? sN2 : 1 };
  enum { N3 = 0 < sN3 ? sN3 : 1 };
  enum { N4 = 0 < sN4 ? sN4 : 1 };
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   , unsigned , unsigned , unsigned 
                         , unsigned , unsigned , unsigned , unsigned ) {}
};

template< long sN1
        , long sN2
        , long sN3
        , long sN4
        , long sN5
        , long sN6
        , long sN7
        >
struct ViewDimension< 0, sN1, sN2, sN3, sN4, sN5, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = sN1 };
  enum { arg_N2 = sN2 };
  enum { arg_N3 = sN3 };
  enum { arg_N4 = sN4 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN1 < 0 ? 1 :
                ( sN2 < 0 ? 2 :
                ( sN3 < 0 ? 3 :
                ( sN4 < 0 ? 4 :
                ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 ))))))) };
  enum { rank_dynamic = 1 };

  size_t N0 ; /* When 1 == rank_dynamic allow N0 >= 2^32 */
  enum { N1 = 0 < sN1 ? sN1 : 1 };
  enum { N2 = 0 < sN2 ? sN2 : 1 };
  enum { N3 = 0 < sN3 ? sN3 : 1 };
  enum { N4 = 0 < sN4 ? sN4 : 1 };
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned , unsigned , unsigned
                         , unsigned , unsigned , unsigned , unsigned )
    : N0( aN0 ) {}
};

template< long sN2
        , long sN3
        , long sN4
        , long sN5
        , long sN6
        , long sN7
        >
struct ViewDimension< 0, 0, sN2, sN3, sN4, sN5, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = sN2 };
  enum { arg_N3 = sN3 };
  enum { arg_N4 = sN4 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN2 < 0 ? 2 :
                ( sN3 < 0 ? 3 :
                ( sN4 < 0 ? 4 :
                ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 )))))) };
  enum { rank_dynamic = 2 };

  size_t N0 ; /* When 2 == rank_dynamic allow N0 >= 2^32 */
  size_t N1 ; /* When 2 == rank_dynamic allow N1 >= 2^32 */
  enum { N2 = 0 < sN2 ? sN2 : 1 };
  enum { N3 = 0 < sN3 ? sN3 : 1 };
  enum { N4 = 0 < sN4 ? sN4 : 1 };
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned , unsigned
                         , unsigned , unsigned , unsigned , unsigned )
    : N0( aN0 ) , N1( aN1 ) {}
};

template< long sN3
        , long sN4
        , long sN5
        , long sN6
        , long sN7
        >
struct ViewDimension< 0, 0, 0, sN3, sN4, sN5, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = sN3 };
  enum { arg_N4 = sN4 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN3 < 0 ? 3 :
                ( sN4 < 0 ? 4 :
                ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 ))))) };
  enum { rank_dynamic = 3 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  enum { N3 = 0 < sN3 ? sN3 : 1 };
  enum { N4 = 0 < sN4 ? sN4 : 1 };
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned
                         , unsigned , unsigned , unsigned , unsigned )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) {}
};

template< long sN4
        , long sN5
        , long sN6
        , long sN7
        >
struct ViewDimension< 0, 0, 0, 0, sN4, sN5, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = 0 };
  enum { arg_N4 = sN4 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN4 < 0 ? 4 :
                ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 )))) };
  enum { rank_dynamic = 4 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  enum { N4 = 0 < sN4 ? sN4 : 1 };
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned aN3
                         , unsigned , unsigned , unsigned , unsigned )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) , N3( aN3 ) {}
};

template< long sN5
        , long sN6
        , long sN7
        >
struct ViewDimension< 0, 0, 0, 0, 0, sN5, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = 0 };
  enum { arg_N4 = 0 };
  enum { arg_N5 = sN5 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN5 < 0 ? 5 :
                ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 ))) };
  enum { rank_dynamic = 5 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  enum { N5 = 0 < sN5 ? sN5 : 1 };
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned aN3
                         , unsigned aN4 , unsigned , unsigned , unsigned )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) , N3( aN3 ) , N4( aN4 ) {}
};

template< long sN6
        , long sN7
        >
struct ViewDimension< 0, 0, 0, 0, 0, 0, sN6, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = 0 };
  enum { arg_N4 = 0 };
  enum { arg_N5 = 0 };
  enum { arg_N6 = sN6 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN6 < 0 ? 6 :
                ( sN7 < 0 ? 7 : 8 )) };
  enum { rank_dynamic = 6 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  enum { N6 = 0 < sN6 ? sN6 : 1 };
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned aN3
                         , unsigned aN4 , unsigned aN5 , unsigned , unsigned )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) , N3( aN3 ) , N4( aN4 ) , N5( aN5 ) {}
};

template< long sN7 >
struct ViewDimension< 0, 0, 0, 0, 0, 0, 0, sN7 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = 0 };
  enum { arg_N4 = 0 };
  enum { arg_N5 = 0 };
  enum { arg_N6 = 0 };
  enum { arg_N7 = sN7 };

  enum { rank = ( sN7 < 0 ? 7 : 8 ) };
  enum { rank_dynamic = 7 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  enum { N7 = 0 < sN7 ? sN7 : 1 };

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned aN3
                         , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) , N3( aN3 ) , N4( aN4 ) , N5( aN5 ) , N6( aN6 ) {}
};

template<>
struct ViewDimension< 0, 0, 0, 0, 0, 0, 0, 0 > {

  enum { arg_N0 = 0 };
  enum { arg_N1 = 0 };
  enum { arg_N2 = 0 };
  enum { arg_N3 = 0 };
  enum { arg_N4 = 0 };
  enum { arg_N5 = 0 };
  enum { arg_N6 = 0 };
  enum { arg_N7 = 0 };

  enum { rank = 8 };
  enum { rank_dynamic = 8 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;

  ViewDimension() = default ;
  ViewDimension( const ViewDimension & ) = default ;
  ViewDimension & operator = ( const ViewDimension & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension( size_t   aN0 , unsigned aN1 , unsigned aN2 , unsigned aN3
                         , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned aN7 )
    : N0( aN0 ) , N1( aN1 ) , N2( aN2 ) , N3( aN3 ) , N4( aN4 ) , N5( aN5 ) , N6( aN6 ) , N7( aN7 ) {}
};

//----------------------------------------------------------------------------

template< class DstDim , class SrcDim >
struct ViewDimensionAssignable ;

template< long dN0 , long dN1 , long dN2 , long dN3 , long dN4 , long dN5 , long dN6 , long dN7 
        , long sN0 , long sN1 , long sN2 , long sN3 , long sN4 , long sN5 , long sN6 , long sN7 >
struct ViewDimensionAssignable< ViewDimension<dN0,dN1,dN2,dN3,dN4,dN5,dN6,dN7>
                              , ViewDimension<sN0,sN1,sN2,sN3,sN4,sN5,sN6,sN7> >
{
  typedef ViewDimension<dN0,dN1,dN2,dN3,dN4,dN5,dN6,dN7>  dst ;
  typedef ViewDimension<sN0,sN1,sN2,sN3,sN4,sN5,sN6,sN7>  src ;

  enum { value = dst::rank == src::rank &&
                 dst::rank_dynamic >= src::rank_dynamic &&
                 ( 0 < dst::rank_dynamic || dN0 == sN0 ) &&
                 ( 1 < dst::rank_dynamic || dN1 == sN1 ) &&
                 ( 2 < dst::rank_dynamic || dN2 == sN2 ) &&
                 ( 3 < dst::rank_dynamic || dN3 == sN3 ) &&
                 ( 4 < dst::rank_dynamic || dN4 == sN4 ) &&
                 ( 5 < dst::rank_dynamic || dN5 == sN5 ) &&
                 ( 6 < dst::rank_dynamic || dN6 == sN6 ) &&
                 ( 7 < dst::rank_dynamic || dN7 == sN7 ) };
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

/** \brief  Given a value type and dimension generate the View data type */
template< class T , class Dim /* ViewDimension */ >
struct ViewDataType {
  enum { R  = Dim::rank };
  enum { RD = Dim::rank_dynamic };

  // Unused static dimensions are set to 1 (instead of 0 or -1L) to avoid compile errors
  // in the 'false' clauses of the std::conditional.

  enum { N0 = 0 < Dim::arg_N0 ? Dim::arg_N0 : 1 };
  enum { N1 = 0 < Dim::arg_N1 ? Dim::arg_N1 : 1 };
  enum { N2 = 0 < Dim::arg_N2 ? Dim::arg_N2 : 1 };
  enum { N3 = 0 < Dim::arg_N3 ? Dim::arg_N3 : 1 };
  enum { N4 = 0 < Dim::arg_N4 ? Dim::arg_N4 : 1 };
  enum { N5 = 0 < Dim::arg_N5 ? Dim::arg_N5 : 1 };
  enum { N6 = 0 < Dim::arg_N6 ? Dim::arg_N6 : 1 };
  enum { N7 = 0 < Dim::arg_N7 ? Dim::arg_N7 : 1 };

  typedef typename std::conditional< R == 0 , T ,
          typename std::conditional< R == 1 ,
            typename std::conditional< RD == 0 , T[N0] , T * >::type ,

          typename std::conditional< R == 2 ,
            typename std::conditional< RD == 0 , T[N0][N1] ,
            typename std::conditional< RD == 1 , T*   [N1] ,
                                                 T**
            >::type >::type ,

          typename std::conditional< R == 3 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2] ,
            typename std::conditional< RD == 1 , T*   [N1][N2] ,
            typename std::conditional< RD == 2 , T**      [N2] ,
                                                 T***
            >::type >::type >::type ,

          typename std::conditional< R == 4 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2][N3] ,
            typename std::conditional< RD == 1 , T*   [N1][N2][N3] ,
            typename std::conditional< RD == 2 , T**      [N2][N3] ,
            typename std::conditional< RD == 3 , T***         [N3] ,
                                                 T****
            >::type >::type >::type >::type ,

          typename std::conditional< R == 5 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2][N3][N4] ,
            typename std::conditional< RD == 1 , T*   [N1][N2][N3][N4] ,
            typename std::conditional< RD == 2 , T**      [N2][N3][N4] ,
            typename std::conditional< RD == 3 , T***         [N3][N4] ,
            typename std::conditional< RD == 4 , T****            [N4] ,
                                                 T*****
            >::type >::type >::type >::type >::type ,

          typename std::conditional< R == 6 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2][N3][N4][N5] ,
            typename std::conditional< RD == 1 , T*   [N1][N2][N3][N4][N5] ,
            typename std::conditional< RD == 2 , T**      [N2][N3][N4][N5] ,
            typename std::conditional< RD == 3 , T***         [N3][N4][N5] ,
            typename std::conditional< RD == 4 , T****            [N4][N5] ,
            typename std::conditional< RD == 5 , T*****               [N5] ,
                                                 T******
            >::type >::type >::type >::type >::type >::type ,

          typename std::conditional< R == 7 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2][N3][N4][N5][N6] ,
            typename std::conditional< RD == 1 , T*   [N1][N2][N3][N4][N5][N6] ,
            typename std::conditional< RD == 2 , T**      [N2][N3][N4][N5][N6] ,
            typename std::conditional< RD == 3 , T***         [N3][N4][N5][N6] ,
            typename std::conditional< RD == 4 , T****            [N4][N5][N6] ,
            typename std::conditional< RD == 5 , T*****               [N5][N6] ,
            typename std::conditional< RD == 6 , T******                  [N6] ,
                                                 T*******
            >::type >::type >::type >::type >::type >::type >::type ,

          typename std::conditional< R == 8 ,
            typename std::conditional< RD == 0 , T[N0][N1][N2][N3][N4][N5][N6][N7] ,
            typename std::conditional< RD == 1 , T*   [N1][N2][N3][N4][N5][N6][N7] ,
            typename std::conditional< RD == 2 , T**      [N2][N3][N4][N5][N6][N7] ,
            typename std::conditional< RD == 3 , T***         [N3][N4][N5][N6][N7] ,
            typename std::conditional< RD == 4 , T****            [N4][N5][N6][N7] ,
            typename std::conditional< RD == 5 , T*****               [N5][N6][N7] ,
            typename std::conditional< RD == 6 , T******                  [N6][N7] ,
            typename std::conditional< RD == 7 , T*******                     [N7] ,
                                                 T******** 
            >::type >::type >::type >::type >::type >::type >::type >::type ,

          void >::type >::type >::type >::type >::type >::type >::type >::type >::type
    type ;
};

/**\brief  Analysis of View data type.
 *
 *  Data type conforms to one of the following patterns :
 *    {const} value_type [][#][#][#]
 *    {const} value_type ***[#][#][#]
 *  Where the sum of counts of '*' and '[#]' is at most ten.
 *
 *  Provide typedef for the ViewDimension<...> and value_type.
 */
template< class T >
struct ViewArrayAnalysis
{
private:
  // std::rank<T>, std::extent<T,i>, and std::remove_all_extents<T>
  // consider "const value_type***" to be the type.

  // Strip away pointers and count them
  typedef typename std::remove_all_extents< T >::type  t_0 ; // brackets removed
  typedef typename std::remove_pointer< t_0   >::type  t_1 ;
  typedef typename std::remove_pointer< t_1   >::type  t_2 ;
  typedef typename std::remove_pointer< t_2   >::type  t_3 ;
  typedef typename std::remove_pointer< t_3   >::type  t_4 ;
  typedef typename std::remove_pointer< t_4   >::type  t_5 ;
  typedef typename std::remove_pointer< t_5   >::type  t_6 ;
  typedef typename std::remove_pointer< t_6   >::type  t_7 ;
  typedef typename std::remove_pointer< t_7   >::type  t_8 ;
  typedef typename std::remove_pointer< t_8   >::type  t_9 ;
  typedef typename std::remove_pointer< t_9   >::type  t_10 ;

  enum { rank_pointer =
    ( ! std::is_pointer< t_0 >::value ? 0 :
    ( ! std::is_pointer< t_1 >::value ? 1 :
    ( ! std::is_pointer< t_2 >::value ? 2 :
    ( ! std::is_pointer< t_3 >::value ? 3 :
    ( ! std::is_pointer< t_4 >::value ? 4 :
    ( ! std::is_pointer< t_5 >::value ? 5 :
    ( ! std::is_pointer< t_6 >::value ? 6 :
    ( ! std::is_pointer< t_7 >::value ? 7 :
    ( ! std::is_pointer< t_8 >::value ? 8 :
    ( ! std::is_pointer< t_9 >::value ? 9 :
    ( ! std::is_pointer< t_10 >::value ? 10 : 0x7fffffff ))))))))))) };

  // The pointer-stripped type t_10 may have been an array typedef of the form 'type[#][#]...'
  // Append those dimensions.

  enum { rank_bracket        = std::rank< T >::value };
  enum { rank_bracket_nested = std::rank< t_10 >::value };
  enum { rank_base           = rank_pointer + rank_bracket };
  enum { rank                = rank_pointer + rank_bracket + rank_bracket_nested };

  static_assert( rank <= 10 , "Maximum ten dimensional array" );

  enum { extent_0 = 0 < rank_base ? std::extent< T , rank_pointer <= 0 ? 0 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 0 ? 0 - rank_base : 10 >::value };

  enum { extent_1 = 1 < rank_base ? std::extent< T , rank_pointer <= 1 ? 1 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 1 ? 1 - rank_base : 10 >::value };

  enum { extent_2 = 2 < rank_base ? std::extent< T , rank_pointer <= 2 ? 2 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 2 ? 2 - rank_base : 10 >::value };

  enum { extent_3 = 3 < rank_base ? std::extent< T , rank_pointer <= 3 ? 3 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 3 ? 3 - rank_base : 10 >::value };

  enum { extent_4 = 4 < rank_base ? std::extent< T , rank_pointer <= 4 ? 4 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 4 ? 4 - rank_base : 10 >::value };

  enum { extent_5 = 5 < rank_base ? std::extent< T , rank_pointer <= 5 ? 5 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 5 ? 5 - rank_base : 10 >::value };

  enum { extent_6 = 6 < rank_base ? std::extent< T , rank_pointer <= 6 ? 6 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 6 ? 6 - rank_base : 10 >::value };

  enum { extent_7 = 7 < rank_base ? std::extent< T , rank_pointer <= 7 ? 7 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 7 ? 7 - rank_base : 10 >::value };

  enum { extent_8 = 8 < rank_base ? std::extent< T , rank_pointer <= 8 ? 8 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 8 ? 8 - rank_base : 10 >::value };

  enum { extent_9 = 9 < rank_base ? std::extent< T , rank_pointer <= 9 ? 9 - rank_pointer : 10 >::value
                                  : std::extent< t_10 , rank_base <= 9 ? 9 - rank_base : 10 >::value };

  typedef typename std::remove_all_extents< t_10 >::type base_type ;

  enum { rank_dynamic = rank_pointer ? rank_pointer : ( ( rank_bracket && extent_0 == 0 ) ? 1 : 0 ) };

public:

  typedef ViewDimension< ( rank <= 0 ? -1L : extent_0 )
                       , ( rank <= 1 ? -1L : extent_1 )
                       , ( rank <= 2 ? -1L : extent_2 )
                       , ( rank <= 3 ? -1L : extent_3 )
                       , ( rank <= 4 ? -1L : extent_4 )
                       , ( rank <= 5 ? -1L : extent_5 )
                       , ( rank <= 6 ? -1L : extent_6 )
                       , ( rank <= 7 ? -1L : extent_7 )
                       > dimension ;

  typedef base_type                                     value_type ;
  typedef typename std::add_const<    base_type >::type const_value_type ;
  typedef typename std::remove_const< base_type >::type non_const_value_type ;

  static_assert( unsigned(dimension::rank) == unsigned(rank) , "" );
  static_assert( unsigned(dimension::rank_dynamic) == unsigned(rank_dynamic) , "" );
};

template< class DataType , class ValueType , class ArrayLayout >
struct ViewDataAnalysis
{
private:

  typedef ViewArrayAnalysis< DataType > array_analysis ;

  // ValueType is opportunity for partial specialization.
  // Must match array analysis when this default template is used.
  static_assert( std::is_same< ValueType , typename array_analysis::non_const_value_type >::value , "" );

public:

  typedef void specialize ; // No specialization

  typedef typename array_analysis::dimension             dimension ;
  typedef typename array_analysis::value_type            value_type ;
  typedef typename array_analysis::const_value_type      const_value_type ;
  typedef typename array_analysis::non_const_value_type  non_const_value_type ;

  // Generate analogous multidimensional array specification type.
  typedef typename ViewDataType<           value_type , dimension >::type  type ;
  typedef typename ViewDataType<     const_value_type , dimension >::type  const_type ;
  typedef typename ViewDataType< non_const_value_type , dimension >::type  non_const_type ;

  // Generate "flattened" multidimensional array specification type.
  typedef type            array_scalar_type ;
  typedef const_type      const_array_scalar_type ;
  typedef non_const_type  non_const_array_scalar_type ;
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

template < class Dimension , class Layout , typename Enable = void >
struct ViewOffset {
  using is_mapping_plugin = std::false_type ;
};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 >= rank OR 0 == rank_dynamic ) : no padding / striding
template < class Dimension >
struct ViewOffset< Dimension , Kokkos::LayoutLeft
                 , typename std::enable_if<( 1 >= Dimension::rank
                                             ||
                                             0 == Dimension::rank_dynamic
                                           )>::type >
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::true_type ;

  typedef size_t             size_type ;
  typedef Dimension          dimension_type ;
  typedef Kokkos::LayoutLeft array_layout ;

  dimension_type m_dim ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    { return i0 + m_dim.N0 * i1 ; }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
  {
    return i0 + m_dim.N0 * ( i1 + m_dim.N1 * i2 );
  }

  //rank 4
  template < typename I0, typename I1, typename I2, typename I3 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3 ) const
  {
    return i0 + m_dim.N0 * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * i3 ));
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4 ) const
  {
    return i0 + m_dim.N0 * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * i4 )));
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5 ) const
  {
    return i0 + m_dim.N0 * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * i5 ))));
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6 ) const
  {
    return i0 + m_dim.N0 * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * (
           i5 + m_dim.N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6, I7 const & i7 ) const
  {
    return i0 + m_dim.N0 * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * (
           i5 + m_dim.N5 * (
           i6 + m_dim.N6 * i7 ))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return true ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_dim.N0 * m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      s[0] = 1 ;
      if ( 0 < dimension_type::rank ) { s[1] = m_dim.N0 ; }
      if ( 1 < dimension_type::rank ) { s[2] = s[1] * m_dim.N1 ; }
      if ( 2 < dimension_type::rank ) { s[3] = s[2] * m_dim.N2 ; }
      if ( 3 < dimension_type::rank ) { s[4] = s[3] * m_dim.N3 ; }
      if ( 4 < dimension_type::rank ) { s[5] = s[4] * m_dim.N4 ; }
      if ( 5 < dimension_type::rank ) { s[6] = s[5] * m_dim.N5 ; }
      if ( 6 < dimension_type::rank ) { s[7] = s[6] * m_dim.N6 ; }
      if ( 7 < dimension_type::rank ) { s[8] = s[7] * m_dim.N7 ; }
    }

  //----------------------------------------

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const &
                      , size_t aN0   , unsigned aN1 , unsigned aN2 , unsigned aN3
                      , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned aN7 )
    : m_dim( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    } 

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutRight , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutLeft and LayoutRight are only compatible when rank == 1" );
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutStride , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutLeft and LayoutStride are only compatible when rank == 1" );
      if ( rhs.m_stride.S0 != 1 ) {
        Kokkos::abort("Kokkos::Experimental::ViewOffset assignment of LayoutLeft from LayoutStride  requires stride == 1" );
      }
    }

  //----------------------------------------
  // Subview construction

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs
                      , const size_t n0
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      )
    : m_dim( n0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( ( 0 == dimension_type::rank ) ||
                     ( 1 == dimension_type::rank && 1 == dimension_type::rank_dynamic && 1 <= DimRHS::rank )
                   , "ViewOffset subview construction requires compatible rank" );
    }
};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template < class Dimension >
struct ViewOffset< Dimension , Kokkos::LayoutLeft
                 , typename std::enable_if<( 1 < Dimension::rank
                                             &&
                                             0 < Dimension::rank_dynamic
                                           )>::type >
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::true_type ;

  typedef size_t             size_type ;
  typedef Dimension          dimension_type ;
  typedef Kokkos::LayoutLeft array_layout ;

  dimension_type m_dim ;
  size_type      m_stride ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    { return i0 + m_stride * i1 ; }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
  {
    return i0 + m_stride * ( i1 + m_dim.N1 * i2 );
  }

  //rank 4
  template < typename I0, typename I1, typename I2, typename I3 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3 ) const
  {
    return i0 + m_stride * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * i3 ));
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4 ) const
  {
    return i0 + m_stride * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * i4 )));
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5 ) const
  {
    return i0 + m_stride * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * i5 ))));
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6 ) const
  {
    return i0 + m_stride * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * (
           i5 + m_dim.N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6, I7 const & i7 ) const
  {
    return i0 + m_stride * (
           i1 + m_dim.N1 * (
           i2 + m_dim.N2 * (
           i3 + m_dim.N3 * (
           i4 + m_dim.N4 * (
           i5 + m_dim.N5 * (
           i6 + m_dim.N6 * i7 ))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return m_stride == m_dim.N0 ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_stride ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_stride * m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_stride * m_dim.N1 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      s[0] = 1 ;
      if ( 0 < dimension_type::rank ) { s[1] = m_stride ; }
      if ( 1 < dimension_type::rank ) { s[2] = s[1] * m_dim.N1 ; }
      if ( 2 < dimension_type::rank ) { s[3] = s[2] * m_dim.N2 ; }
      if ( 3 < dimension_type::rank ) { s[4] = s[3] * m_dim.N3 ; }
      if ( 4 < dimension_type::rank ) { s[5] = s[4] * m_dim.N4 ; }
      if ( 5 < dimension_type::rank ) { s[6] = s[5] * m_dim.N5 ; }
      if ( 6 < dimension_type::rank ) { s[7] = s[6] * m_dim.N6 ; }
      if ( 7 < dimension_type::rank ) { s[8] = s[7] * m_dim.N7 ; }
    }

  //----------------------------------------

private:

  template< unsigned TrivialScalarSize >
  struct Padding {
    enum { div = TrivialScalarSize == 0 ? 0 : Kokkos::Impl::MEMORY_ALIGNMENT / ( TrivialScalarSize ? TrivialScalarSize : 1 ) };
    enum { mod = TrivialScalarSize == 0 ? 0 : Kokkos::Impl::MEMORY_ALIGNMENT % ( TrivialScalarSize ? TrivialScalarSize : 1 ) };

    // If memory alignment is a multiple of the trivial scalar size then attempt to align.
    enum { align = 0 != TrivialScalarSize && 0 == mod ? div : 0 };
    enum { div_ok = div ? div : 1 }; // To valid modulo zero in constexpr

    KOKKOS_INLINE_FUNCTION
    static constexpr size_t stride( size_t const N )
      {
        return ( align && ( Kokkos::Impl::MEMORY_ALIGNMENT_THRESHOLD * align < N ) && ( N % div_ok ) )
               ? N + align - ( N % div_ok ) : N ;
      }
  };

public:

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  /* Enable padding for trivial scalar types with non-zero trivial scalar size */
  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const & padding_type_size
                      , size_t aN0   , unsigned aN1 , unsigned aN2 , unsigned aN3
                      , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned aN7 )
    : m_dim( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
    , m_stride( Padding<TrivialScalarSize>::stride( aN0 ) )
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    , m_stride( rhs.stride_1() )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    } 

  //----------------------------------------
  // Subview construction

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs
                      , const size_t aN0
                      , const size_t aN1
                      , const size_t aN2
                      , const size_t aN3
                      , const size_t aN4
                      , const size_t aN5
                      , const size_t aN6
                      , const size_t aN7
                      )
    : m_dim( aN0
           , ( 1 < DimRHS::rank && aN1 ? aN1 :
             ( 2 < DimRHS::rank && aN2 ? aN2 :
             ( 3 < DimRHS::rank && aN3 ? aN3 :
             ( 4 < DimRHS::rank && aN4 ? aN4 :
             ( 5 < DimRHS::rank && aN5 ? aN5 :
             ( 6 < DimRHS::rank && aN6 ? aN6 :
             ( 7 < DimRHS::rank && aN7 ? aN7 : 0 )))))))
           , 0, 0, 0, 0, 0, 0 )
    , m_stride( ( 1 < DimRHS::rank && aN1 ? rhs.stride_1() :
                ( 2 < DimRHS::rank && aN2 ? rhs.stride_2() :
                ( 3 < DimRHS::rank && aN3 ? rhs.stride_3() :
                ( 4 < DimRHS::rank && aN4 ? rhs.stride_4() :
                ( 5 < DimRHS::rank && aN5 ? rhs.stride_5() :
                ( 6 < DimRHS::rank && aN6 ? rhs.stride_6() :
                ( 7 < DimRHS::rank && aN7 ? rhs.stride_7() : 0 ))))))) )
    {
      // This subview must be 2 == rank and 2 == rank_dynamic
      // due to only having stride #0.
      // The source dimension #0 must be non-zero for stride-one leading dimension.
      // At most subsequent dimension can be non-zero.

      static_assert( ( 2 == dimension_type::rank ) &&
                     ( 2 == dimension_type::rank_dynamic ) &&
                     ( 2 <= DimRHS::rank )
                   , "ViewOffset subview construction requires compatible rank" );
    }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 >= rank OR 0 == rank_dynamic ) : no padding / striding
template < class Dimension >
struct ViewOffset< Dimension , Kokkos::LayoutRight
                 , typename std::enable_if<( 1 >= Dimension::rank
                                             ||
                                             0 == Dimension::rank_dynamic
                                           )>::type >
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::true_type ;

  typedef size_t              size_type ;
  typedef Dimension           dimension_type ;
  typedef Kokkos::LayoutRight array_layout ;

  dimension_type m_dim ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    { return i1 + m_dim.N1 * i0 ; }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
  {
    return i2 + m_dim.N2 * ( i1 + m_dim.N1 * ( i0 ));
  }

  //rank 4
  template < typename I0, typename I1, typename I2, typename I3 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3 ) const
  {
    return i3 + m_dim.N3 * (
           i2 + m_dim.N2 * (
           i1 + m_dim.N1 * ( i0 )));
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4 ) const
  {
    return i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * (
           i1 + m_dim.N1 * ( i0 ))));
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5 ) const
  {
    return i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * (
           i1 + m_dim.N1 * ( i0 )))));
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6 ) const
  {
    return i6 + m_dim.N6 * (
           i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * (
           i1 + m_dim.N1 * ( i0 ))))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6, I7 const & i7 ) const
  {
    return i7 + m_dim.N7 * (
           i6 + m_dim.N6 * (
           i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * (
           i1 + m_dim.N1 * ( i0 )))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return true ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_dim.N7 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_dim.N7 * m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2 * m_dim.N1 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < dimension_type::rank ) { s[7] = n ; n *= m_dim.N7 ; }
      if ( 6 < dimension_type::rank ) { s[6] = n ; n *= m_dim.N6 ; }
      if ( 5 < dimension_type::rank ) { s[5] = n ; n *= m_dim.N5 ; }
      if ( 4 < dimension_type::rank ) { s[4] = n ; n *= m_dim.N4 ; }
      if ( 3 < dimension_type::rank ) { s[3] = n ; n *= m_dim.N3 ; }
      if ( 2 < dimension_type::rank ) { s[2] = n ; n *= m_dim.N2 ; }
      if ( 1 < dimension_type::rank ) { s[1] = n ; n *= m_dim.N1 ; }
      if ( 0 < dimension_type::rank ) { s[0] = n ; }
      s[dimension_type::rank] = n * m_dim.N0 ;
    }

  //----------------------------------------

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const &
                      , size_t aN0   , unsigned aN1 , unsigned aN2 , unsigned aN3
                      , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned aN7 )
    : m_dim( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutRight , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    } 

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutRight and LayoutLeft are only compatible when rank == 1" );
    }

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutStride , void > & rhs )
    : m_dim( rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( DimRHS::rank == 1 && dimension_type::rank == 1 && dimension_type::rank_dynamic == 1
                   , "ViewOffset LayoutLeft and LayoutStride are only compatible when rank == 1" );
      if ( rhs.m_stride.S0 != 1 ) {
        Kokkos::abort("Kokkos::Experimental::ViewOffset assignment of LayoutRight from LayoutStride  requires stride == 1" );
      }
    }

  //----------------------------------------
  // Subview construction

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutRight , void > & rhs
                      , const size_t n0
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      , const size_t
                      )
    : m_dim( n0, 0, 0, 0, 0, 0, 0, 0 )
    {
      static_assert( ( 0 == dimension_type::rank ) ||
                     ( 1 == dimension_type::rank && 1 == dimension_type::rank_dynamic && 1 <= DimRHS::rank )
                   , "ViewOffset subview construction requires compatible rank" );
    }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template < class Dimension >
struct ViewOffset< Dimension , Kokkos::LayoutRight
                 , typename std::enable_if<( 1 < Dimension::rank
                                             &&
                                             0 < Dimension::rank_dynamic
                                           )>::type >
{
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::true_type ;

  typedef size_t               size_type ;
  typedef Dimension            dimension_type ;
  typedef Kokkos::LayoutRight  array_layout ;

  dimension_type m_dim ;
  size_type      m_stride ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
  { return i1 + i0 * m_stride ; }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
  { return i2 + m_dim.N2 * ( i1 ) + i0 * m_stride ; }

  //rank 4
  template < typename I0, typename I1, typename I2, typename I3 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3 ) const
  {
    return i3 + m_dim.N3 * (
           i2 + m_dim.N2 * ( i1 )) +
           i0 * m_stride ;
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4 ) const
  {
    return i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * ( i1 ))) +
           i0 * m_stride ;
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5 ) const
  {
    return i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * ( i1 )))) +
           i0 * m_stride ;
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6 ) const
  {
    return i6 + m_dim.N6 * (
           i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * ( i1 ))))) +
           i0 * m_stride ;
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6, I7 const & i7 ) const
  {
    return i7 + m_dim.N7 * (
           i6 + m_dim.N6 * (
           i5 + m_dim.N5 * (
           i4 + m_dim.N4 * (
           i3 + m_dim.N3 * (
           i2 + m_dim.N2 * ( i1 )))))) +
           i0 * m_stride ;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    { return m_dim.N0 * m_stride ; }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_stride == m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2 * m_dim.N1 ; }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_dim.N7 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_dim.N7 * m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return m_stride ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < dimension_type::rank ) { s[7] = n ; n *= m_dim.N7 ; }
      if ( 6 < dimension_type::rank ) { s[6] = n ; n *= m_dim.N6 ; }
      if ( 5 < dimension_type::rank ) { s[5] = n ; n *= m_dim.N5 ; }
      if ( 4 < dimension_type::rank ) { s[4] = n ; n *= m_dim.N4 ; }
      if ( 3 < dimension_type::rank ) { s[3] = n ; n *= m_dim.N3 ; }
      if ( 2 < dimension_type::rank ) { s[2] = n ; n *= m_dim.N2 ; }
      if ( 1 < dimension_type::rank ) { s[1] = n ; }
      if ( 0 < dimension_type::rank ) { s[0] = m_stride ; }
      s[dimension_type::rank] = m_stride * m_dim.N0 ;
    }

  //----------------------------------------

private:

  template< unsigned TrivialScalarSize >
  struct Padding {
    enum { div = TrivialScalarSize == 0 ? 0 : Kokkos::Impl::MEMORY_ALIGNMENT / ( TrivialScalarSize ? TrivialScalarSize : 1 ) };
    enum { mod = TrivialScalarSize == 0 ? 0 : Kokkos::Impl::MEMORY_ALIGNMENT % ( TrivialScalarSize ? TrivialScalarSize : 1 ) };

    // If memory alignment is a multiple of the trivial scalar size then attempt to align.
    enum { align = 0 != TrivialScalarSize && 0 == mod ? div : 0 };
    enum { div_ok = div ? div : 1 }; // To valid modulo zero in constexpr

    KOKKOS_INLINE_FUNCTION
    static constexpr size_t stride( size_t const N )
    {
      return ( align && ( Kokkos::Impl::MEMORY_ALIGNMENT_THRESHOLD * align < N ) && ( N % div_ok ) )
             ? N + align - ( N % div_ok ) : N ;
    }
  };

public:

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  /* Enable padding for trivial scalar types with non-zero trivial scalar size.  */
  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const & padding_type_size
                      , size_t aN0   , unsigned aN1 , unsigned aN2 , unsigned aN3
                      , unsigned aN4 , unsigned aN5 , unsigned aN6 , unsigned aN7 )
    : m_dim( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
    , m_stride( Padding<TrivialScalarSize>::
                  stride( /* 2 <= rank */
                          m_dim.N1 * ( dimension_type::rank == 2 ? 1 :
                          m_dim.N2 * ( dimension_type::rank == 3 ? 1 :
                          m_dim.N3 * ( dimension_type::rank == 4 ? 1 :
                          m_dim.N4 * ( dimension_type::rank == 5 ? 1 :
                          m_dim.N5 * ( dimension_type::rank == 6 ? 1 :
                          m_dim.N6 * ( dimension_type::rank == 7 ? 1 : m_dim.N7 )))))) ))
    {}

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutLeft , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    , m_stride( rhs.stride_0() )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    } 

  //----------------------------------------
  // Subview construction
  // Last dimension must be non-zero

  template< class DimRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , Kokkos::LayoutRight , void > & rhs
                      , const size_t aN0
                      , const size_t aN1
                      , const size_t aN2
                      , const size_t aN3
                      , const size_t aN4
                      , const size_t aN5
                      , const size_t aN6
                      , const size_t aN7
                      )
    : m_dim( // N0 == First non-zero dimension before the last dimension.
             ( 1 < DimRHS::rank && aN0 ? aN0 :
             ( 2 < DimRHS::rank && aN1 ? aN1 :
             ( 3 < DimRHS::rank && aN2 ? aN2 :
             ( 4 < DimRHS::rank && aN3 ? aN3 :
             ( 5 < DimRHS::rank && aN4 ? aN4 :
             ( 6 < DimRHS::rank && aN5 ? aN5 :
             ( 7 < DimRHS::rank && aN6 ? aN6 : 0 )))))))
           , // N1 == Last dimension.
             ( 2 == DimRHS::rank ? aN1 :
             ( 3 == DimRHS::rank ? aN2 :
             ( 4 == DimRHS::rank ? aN3 :
             ( 5 == DimRHS::rank ? aN4 :
             ( 6 == DimRHS::rank ? aN5 :
             ( 7 == DimRHS::rank ? aN6 : aN7 ))))))
           , 0, 0, 0, 0, 0, 0 )
    , m_stride( ( 1 < DimRHS::rank && aN0 ? rhs.stride_0() :
                ( 2 < DimRHS::rank && aN1 ? rhs.stride_1() :
                ( 3 < DimRHS::rank && aN2 ? rhs.stride_2() :
                ( 4 < DimRHS::rank && aN3 ? rhs.stride_3() :
                ( 5 < DimRHS::rank && aN4 ? rhs.stride_4() :
                ( 6 < DimRHS::rank && aN5 ? rhs.stride_5() :
                ( 7 < DimRHS::rank && aN6 ? rhs.stride_6() : 0 ))))))) )
    {
      // This subview must be 2 == rank and 2 == rank_dynamic
      // due to only having stride #0.
      // The source dimension #0 must be non-zero for stride-one leading dimension.
      // At most subsequent dimension can be non-zero.

      static_assert( ( 2 == dimension_type::rank ) &&
                     ( 2 == dimension_type::rank_dynamic ) &&
                     ( 2 <= DimRHS::rank )
                   , "ViewOffset subview construction requires compatible rank" );
    }
};

//----------------------------------------------------------------------------
/* Strided array layout only makes sense for 0 < rank */

template< unsigned Rank >
struct ViewStride ;

template<>
struct ViewStride<1> {
  size_t S0 ;
  enum { S1 = 0 , S2 = 0 , S3 = 0 , S4 = 0 , S5 = 0 , S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t , size_t , size_t
                      , size_t , size_t , size_t , size_t )
    : S0( aS0 )
    {}
};

template<>
struct ViewStride<2> {
  size_t S0 , S1 ;
  enum { S2 = 0 , S3 = 0 , S4 = 0 , S5 = 0 , S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t , size_t
                      , size_t , size_t , size_t , size_t )
    : S0( aS0 ) , S1( aS1 )
    {}
};

template<>
struct ViewStride<3> {
  size_t S0 , S1 , S2 ;
  enum { S3 = 0 , S4 = 0 , S5 = 0 , S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t
                      , size_t , size_t , size_t , size_t )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 )
    {}
};

template<>
struct ViewStride<4> {
  size_t S0 , S1 , S2 , S3 ;
  enum { S4 = 0 , S5 = 0 , S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t aS3
                      , size_t , size_t , size_t , size_t )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 ) , S3( aS3 )
    {}
};

template<>
struct ViewStride<5> {
  size_t S0 , S1 , S2 , S3 , S4 ;
  enum { S5 = 0 , S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t aS3
                      , size_t aS4 , size_t , size_t , size_t )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 ) , S3( aS3 )
    , S4( aS4 )
    {}
};

template<>
struct ViewStride<6> {
  size_t S0 , S1 , S2 , S3 , S4 , S5 ;
  enum { S6 = 0 , S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t aS3
                      , size_t aS4 , size_t aS5 , size_t , size_t )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 ) , S3( aS3 )
    , S4( aS4 ) , S5( aS5 )
    {}
};

template<>
struct ViewStride<7> {
  size_t S0 , S1 , S2 , S3 , S4 , S5 , S6 ;
  enum { S7 = 0 };

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t aS3
                      , size_t aS4 , size_t aS5 , size_t aS6 , size_t )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 ) , S3( aS3 )
    , S4( aS4 ) , S5( aS5 ) , S6( aS6 )
    {}
};

template<>
struct ViewStride<8> {
  size_t S0 , S1 , S2 , S3 , S4 , S5 , S6 , S7 ;

  ViewStride() = default ;
  ViewStride( const ViewStride & ) = default ;
  ViewStride & operator = ( const ViewStride & ) = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride( size_t aS0 , size_t aS1 , size_t aS2 , size_t aS3
                      , size_t aS4 , size_t aS5 , size_t aS6 , size_t aS7 )
    : S0( aS0 ) , S1( aS1 ) , S2( aS2 ) , S3( aS3 )
    , S4( aS4 ) , S5( aS5 ) , S6( aS6 ) , S7( aS7 )
    {}
};

template < class Dimension >
struct ViewOffset< Dimension , Kokkos::LayoutStride
                 , typename std::enable_if<( 0 < Dimension::rank )>::type >
{
private:
  typedef ViewStride< Dimension::rank >  stride_type ;
public:

  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::true_type ;

  typedef size_t                size_type ;
  typedef Dimension             dimension_type ;
  typedef Kokkos::LayoutStride  array_layout ;

  dimension_type  m_dim ;
  stride_type     m_stride ;

  //----------------------------------------

  // rank 1
  template< typename I0 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 ) const
  {
    return i0 * m_stride.S0 ;
  }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 ;
  }

  //rank 3
  template < typename I0, typename I1, typename I2 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 ;
  }

  //rank 4
  template < typename I0, typename I1, typename I2, typename I3 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 +
           i3 * m_stride.S3 ;
  }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 +
           i3 * m_stride.S3 +
           i4 * m_stride.S4 ;
  }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 +
           i3 * m_stride.S3 +
           i4 * m_stride.S4 +
           i5 * m_stride.S5 ;
  }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 +
           i3 * m_stride.S3 +
           i4 * m_stride.S4 +
           i5 * m_stride.S5 +
           i6 * m_stride.S6 ;
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
           , typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0, I1 const & i1, I2 const & i2, I3 const & i3
                      , I4 const & i4, I5 const & i5, I6 const & i6, I7 const & i7 ) const
  {
    return i0 * m_stride.S0 +
           i1 * m_stride.S1 +
           i2 * m_stride.S2 +
           i3 * m_stride.S3 +
           i4 * m_stride.S4 +
           i5 * m_stride.S5 +
           i6 * m_stride.S6 +
           i7 * m_stride.S7 ;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return m_dim.N2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return m_dim.N3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return m_dim.N4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return m_dim.N5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return m_dim.N6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return m_dim.N7 ; }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const
    { return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7 ; }

private:

  KOKKOS_INLINE_FUNCTION
  static constexpr size_type Max( size_type lhs , size_type rhs )
    { return lhs < rhs ? rhs : lhs ; }

public:

  /* Span of the range space, largest stride * dimension */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const
    {
      return Max( m_dim.N0 * m_stride.S0 ,
             Max( m_dim.N1 * m_stride.S1 ,
             Max( m_dim.N2 * m_stride.S2 ,
             Max( m_dim.N3 * m_stride.S3 ,
             Max( m_dim.N4 * m_stride.S4 ,
             Max( m_dim.N5 * m_stride.S5 ,
             Max( m_dim.N6 * m_stride.S6 ,
                  m_dim.N7 * m_stride.S7 )))))));
    }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return span() == size(); }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return m_stride.S0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return m_stride.S1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return m_stride.S2 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return m_stride.S3 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return m_stride.S4 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return m_stride.S5 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return m_stride.S6 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return m_stride.S7 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      if ( 0 < dimension_type::rank ) { s[0] = m_stride.S0 ; }
      if ( 1 < dimension_type::rank ) { s[1] = m_stride.S1 ; }
      if ( 2 < dimension_type::rank ) { s[2] = m_stride.S2 ; }
      if ( 3 < dimension_type::rank ) { s[3] = m_stride.S3 ; }
      if ( 4 < dimension_type::rank ) { s[4] = m_stride.S4 ; }
      if ( 5 < dimension_type::rank ) { s[5] = m_stride.S5 ; }
      if ( 6 < dimension_type::rank ) { s[6] = m_stride.S6 ; }
      if ( 7 < dimension_type::rank ) { s[7] = m_stride.S7 ; }
      s[dimension_type::rank] = span();
    }

  //----------------------------------------

  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  KOKKOS_INLINE_FUNCTION
  ViewOffset( const Kokkos::LayoutStride & rhs )
    : m_dim( rhs.dimension[0] , rhs.dimension[1] , rhs.dimension[2] , rhs.dimension[3]
           , rhs.dimension[4] , rhs.dimension[5] , rhs.dimension[6] , rhs.dimension[7] )
    , m_stride( rhs.stride[0] , rhs.stride[1] , rhs.stride[2] , rhs.stride[3]
              , rhs.stride[4] , rhs.stride[5] , rhs.stride[6] , rhs.stride[7] )
    {}

  template< class DimRHS , class LayoutRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , LayoutRHS , void > & rhs )
    : m_dim( rhs.m_dim.N0 , rhs.m_dim.N1 , rhs.m_dim.N2 , rhs.m_dim.N3 
           , rhs.m_dim.N4 , rhs.m_dim.N5 , rhs.m_dim.N6 , rhs.m_dim.N7 )
    , m_stride( rhs.stride_0() , rhs.stride_1() , rhs.stride_2() , rhs.stride_3()
              , rhs.stride_4() , rhs.stride_5() , rhs.stride_6() , rhs.stride_7() )
    {
      static_assert( int(DimRHS::rank) == int(dimension_type::rank) , "ViewOffset assignment requires equal rank" );
      // Also requires equal static dimensions ...
    }

  //----------------------------------------
  // Subview construction

private:

  KOKKOS_INLINE_FUNCTION
  static constexpr unsigned
    count_non_zero( const size_t aN0 = 0 
                  , const size_t aN1 = 0 
                  , const size_t aN2 = 0 
                  , const size_t aN3 = 0 
                  , const size_t aN4 = 0 
                  , const size_t aN5 = 0 
                  , const size_t aN6 = 0 
                  , const size_t aN7 = 0 
                  )
    {
      return ( aN0 ? 1 : 0 ) +
             ( aN1 ? 1 : 0 ) +
             ( aN2 ? 1 : 0 ) +
             ( aN3 ? 1 : 0 ) +
             ( aN4 ? 1 : 0 ) +
             ( aN5 ? 1 : 0 ) +
             ( aN6 ? 1 : 0 ) +
             ( aN7 ? 1 : 0 );
    }

  template< unsigned Rank , unsigned I >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t
    get_non_zero( const size_t aN0
                , const size_t aN1
                , const size_t aN2
                , const size_t aN3
                , const size_t aN4
                , const size_t aN5
                , const size_t aN6
                , const size_t aN7
                )
    {
      return ( 0 < Rank && I < 1                                                     && aN0 ? aN0 :
             ( 1 < Rank && I < 2 && I == count_non_zero(aN0)                         && aN1 ? aN1 :
             ( 2 < Rank && I < 3 && I == count_non_zero(aN0,aN1)                     && aN2 ? aN2 :
             ( 3 < Rank && I < 4 && I == count_non_zero(aN0,aN1,aN2)                 && aN3 ? aN3 :
             ( 4 < Rank && I < 5 && I == count_non_zero(aN0,aN1,aN2,aN3)             && aN4 ? aN4 :
             ( 5 < Rank && I < 6 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4)         && aN5 ? aN5 :
             ( 6 < Rank && I < 7 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4,aN5)     && aN6 ? aN6 :
             ( 7 < Rank && I < 8 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4,aN5,aN6) && aN7 ? aN7 : 0 ))))))));
    }
  
  template< unsigned Rank , unsigned I , class DimRHS , class LayoutRHS >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t
    get_non_zero( const size_t aN0 , const size_t aN1 , const size_t aN2 , const size_t aN3
                , const size_t aN4 , const size_t aN5 , const size_t aN6 , const size_t aN7
                , const ViewOffset< DimRHS , LayoutRHS , void > & rhs )
    {
      return ( 0 < Rank && I < 1                                                     && aN0 ? rhs.stride_0() :
             ( 1 < Rank && I < 2 && I == count_non_zero(aN0)                         && aN1 ? rhs.stride_1() :
             ( 2 < Rank && I < 3 && I == count_non_zero(aN0,aN1)                     && aN2 ? rhs.stride_2() :
             ( 3 < Rank && I < 4 && I == count_non_zero(aN0,aN1,aN2)                 && aN3 ? rhs.stride_3() :
             ( 4 < Rank && I < 5 && I == count_non_zero(aN0,aN1,aN2,aN3)             && aN4 ? rhs.stride_4() :
             ( 5 < Rank && I < 6 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4)         && aN5 ? rhs.stride_5() :
             ( 6 < Rank && I < 7 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4,aN5)     && aN6 ? rhs.stride_6() :
             ( 7 < Rank && I < 8 && I == count_non_zero(aN0,aN1,aN2,aN3,aN4,aN5,aN6) && aN7 ? rhs.stride_7() : 0 ))))))));
    }
  

public:

  template< class DimRHS , class LayoutRHS >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( const ViewOffset< DimRHS , LayoutRHS , void > & rhs
                      , const size_t aN0
                      , const size_t aN1
                      , const size_t aN2
                      , const size_t aN3
                      , const size_t aN4
                      , const size_t aN5
                      , const size_t aN6
                      , const size_t aN7
                      )
    // Contract the non-zero dimensions
    : m_dim( ViewOffset::template get_non_zero<DimRHS::rank,0>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,1>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,2>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,3>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,4>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,5>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,6>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           , ViewOffset::template get_non_zero<DimRHS::rank,7>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7 )
           )
    , m_stride( ViewOffset::template get_non_zero<DimRHS::rank,0>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,1>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,2>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,3>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,4>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,5>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,6>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              , ViewOffset::template get_non_zero<DimRHS::rank,7>( aN0, aN1, aN2, aN3, aN4, aN5, aN6, aN7, rhs )
              )
    {
    }

  //----------------------------------------
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ALL_t {
  KOKKOS_INLINE_FUNCTION
  constexpr const ALL_t & operator()() const { return *this ; }
};

template< class T >
struct ViewOffsetRange {

  static_assert( std::is_integral<T>::value , "Non-range must be an integral type" );

  enum { is_range = false };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const , T const & ) { return 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( T const & i ) { return size_t(i) ; }
};

template<>
struct ViewOffsetRange<void> {
  enum { is_range = false };
};

template<>
struct ViewOffsetRange< Kokkos::Experimental::Impl::ALL_t > {
  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , Experimental::Impl::ALL_t const & ) { return n ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( Experimental::Impl::ALL_t const & ) { return 0 ; }
};

template< typename iType >
struct ViewOffsetRange< std::pair<iType,iType> > {

  static_assert( std::is_integral<iType>::value , "Range bounds must be an integral type" );

  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , std::pair<iType,iType> const & r )
    { return ( size_t(r.first) < size_t(r.second) && size_t(r.second) <= n ) ? size_t(r.second) - size_t(r.first) : 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( std::pair<iType,iType> const & r ) { return size_t(r.first) ; }
};

template< typename iType >
struct ViewOffsetRange< Kokkos::pair<iType,iType> > {

  static_assert( std::is_integral<iType>::value , "Range bounds must be an integral type" );

  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , Kokkos::pair<iType,iType> const & r )
    { return ( size_t(r.first) < size_t(r.second) && size_t(r.second) <= n ) ? size_t(r.second) - size_t(r.first) : 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( Kokkos::pair<iType,iType> const & r ) { return size_t(r.first) ; }
};

template< typename iType >
struct ViewOffsetRange< std::initializer_list< iType > > {

  static_assert( std::is_integral<iType>::value , "Range bounds must be an integral type" );

  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , std::initializer_list< iType > const & r )
    {
      return ( size_t(r.begin()[0]) < size_t(r.begin()[1]) && size_t(r.begin()[1]) <= n )
             ? size_t(r.begin()[1]) - size_t(r.begin()[0]) : 0 ;
    }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( std::initializer_list< iType > const & r ) { return size_t(r.begin()[0]) ; }
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

/** \brief  ViewDataHandle provides the type of the 'data handle' which the view
 *          uses to access data with the [] operator. It also provides
 *          an allocate function and a function to extract a raw ptr from the
 *          data handle. ViewDataHandle also defines an enum ReferenceAble which
 *          specifies whether references/pointers to elements can be taken and a
 *          'return_type' which is what the view operators will give back.
 *          Specialisation of this object allows three things depending
 *          on ViewTraits and compiler options:
 *          (i)   Use special allocator (e.g. huge pages/small pages and pinned memory)
 *          (ii)  Use special data handle type (e.g. add Cuda Texture Object)
 *          (iii) Use special access intrinsics (e.g. texture fetch and non-caching loads)
 */
template< class Traits , class Enable = void >
struct ViewDataHandle {

  typedef typename Traits::value_type   value_type  ;
  typedef typename Traits::value_type * handle_type ;
  typedef typename Traits::value_type & return_type ;
  typedef Kokkos::Experimental::Impl::SharedAllocationTracker  track_type  ;

  KOKKOS_INLINE_FUNCTION
  static handle_type assign( value_type * arg_data_ptr
                           , track_type const & /*arg_tracker*/ )
  {
    return handle_type( arg_data_ptr );
  }
};

template< class Traits >
struct ViewDataHandle< Traits ,
  typename std::enable_if<( std::is_same< typename Traits::non_const_value_type
                                        , typename Traits::value_type >::value
                            &&
                            std::is_same< typename Traits::specialize , void >::value
                            &&
                            Traits::memory_traits::Atomic
                          )>::type >
{
  typedef typename Traits::value_type  value_type ;
  typedef typename Kokkos::Impl::AtomicViewDataHandle< Traits >  handle_type ;
  typedef typename Kokkos::Impl::AtomicDataElement< Traits >     return_type ;
  typedef Kokkos::Experimental::Impl::SharedAllocationTracker    track_type  ;

  KOKKOS_INLINE_FUNCTION
  static handle_type assign( value_type * arg_data_ptr
                           , track_type const & /*arg_tracker*/ )
  {
    return handle_type( arg_data_ptr );
  }
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

template< class Traits
        , bool R0 = false
        , bool R1 = false
        , bool R2 = false
        , bool R3 = false
        , bool R4 = false
        , bool R5 = false
        , bool R6 = false
        , bool R7 = false
        , typename Enable = void >
struct SubviewMapping ;

/** \brief  View mapping for non-specialized data type and standard layout */
template< class Traits >
class ViewMapping< Traits , void ,
  typename std::enable_if<(
    std::is_same< typename Traits::specialize , void >::value
    &&
    ViewOffset< typename Traits::dimension
              , typename Traits::array_layout
              , void >::is_mapping_plugin::value
  )>::type >
{
private:

  template< class , class , typename > friend class ViewMapping ;
  template< class , bool , bool , bool , bool , bool , bool , bool , bool , class > friend struct SubviewMapping ;
  template< class , class , class , class > friend class Kokkos::Experimental::View ;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  typedef typename ViewDataHandle< Traits >::handle_type  handle_type ;

  handle_type  m_handle ;
  offset_type  m_offset ;

  KOKKOS_INLINE_FUNCTION
  ViewMapping( const handle_type & arg_handle , const offset_type & arg_offset )
    : m_handle( arg_handle )
    , m_offset( arg_offset )
    {}

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const { return m_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return m_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return m_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return m_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return m_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return m_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return m_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return m_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  using is_regular = typename offset_type::is_regular ;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const { return m_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const { return m_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const { return m_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const { return m_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const { return m_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const { return m_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const { return m_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const { return m_offset.stride_7(); }

  //----------------------------------------
  // Range span

  /** \brief  Span of the mapped range */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return m_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const { return m_offset.span_is_contiguous(); }

  typedef typename ViewDataHandle< Traits >::return_type  reference_type ;

  /** \brief  If data references are lvalue_reference than can query pointer to memory */
  KOKKOS_INLINE_FUNCTION constexpr typename Traits::value_type * data() const
    {
      typedef typename Traits::value_type * ptr_type ;

      return std::is_lvalue_reference< reference_type >::value
             ? (ptr_type) m_handle
             : (ptr_type) 0 ;
    }

  //----------------------------------------
  // The View class performs all rank and bounds checking before
  // calling these element reference methods.

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const { return m_handle[0]; }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    ! std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const { return m_handle[i0]; }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename
    std::enable_if< std::is_integral<I0>::value &&
                    std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
                  , reference_type >::type
  reference( const I0 & i0 ) const { return m_handle[ m_offset(i0) ]; }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 ) const
    { return m_handle[ m_offset(i0,i1) ]; }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { return m_handle[ m_offset(i0,i1,i2) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { return m_handle[ m_offset(i0,i1,i2,i3) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 ) const
    { return m_handle[ m_offset(i0,i1,i2,i3,i4) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 ) const
    { return m_handle[ m_offset(i0,i1,i2,i3,i4,i5) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { return m_handle[ m_offset(i0,i1,i2,i3,i4,i5,i6) ]; }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 , typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
                          , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 ) const
    { return m_handle[ m_offset(i0,i1,i2,i3,i4,i5,i6,i7) ]; }

  //----------------------------------------

private:

  enum { MemorySpanMask = 8 - 1 /* Force alignment on 8 byte boundary */ };
  enum { MemorySpanSize = sizeof(typename Traits::value_type) };

public:

  /** \brief  Span, in bytes, of the referenced memory */
  KOKKOS_INLINE_FUNCTION constexpr size_t memory_span() const
    {
      return ( m_offset.span() * sizeof(typename Traits::value_type) + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  /** \brief  Span, in bytes, of the required memory */
  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span( const std::integral_constant<bool,AllowPadding> &
                                     , const size_t N0 , const size_t N1 , const size_t N2 , const size_t N3
                                      , const size_t N4 , const size_t N5 , const size_t N6 , const size_t N7 )
    {
      typedef std::integral_constant< unsigned , AllowPadding ? MemorySpanSize : 0 >  padding ;
      return ( offset_type( padding(), N0, N1, N2, N3, N4, N5, N6, N7 ).span() * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  /** \brief  Span, in bytes, of the required memory */
  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span( const std::integral_constant<bool,AllowPadding> &
                                       , const typename Traits::array_layout & layout )
    {
      return ( offset_type( layout ).span() * MemorySpanSize + MemorySpanMask ) & ~size_t(MemorySpanMask);
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION ~ViewMapping() {}
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_handle(), m_offset() {}
  KOKKOS_INLINE_FUNCTION ViewMapping( const ViewMapping & rhs )
    : m_handle( rhs.m_handle ), m_offset( rhs.m_offset ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( const ViewMapping & rhs )
    { m_handle = rhs.m_handle ; m_offset = rhs.m_offset ; return *this ; }

  KOKKOS_INLINE_FUNCTION ViewMapping( ViewMapping && rhs )
    : m_handle( rhs.m_handle ), m_offset( rhs.m_offset ) {}
  KOKKOS_INLINE_FUNCTION ViewMapping & operator = ( ViewMapping && rhs )
    { m_handle = rhs.m_handle ; m_offset = rhs.m_offset ; return *this ; }

  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  ViewMapping( void * ptr
             , const std::integral_constant<bool,AllowPadding> &
             , const size_t N0 , const size_t N1 , const size_t N2 , const size_t N3
             , const size_t N4 , const size_t N5 , const size_t N6 , const size_t N7 )
    : m_handle( reinterpret_cast< handle_type >( ptr ) )
    , m_offset( std::integral_constant< unsigned , AllowPadding ? sizeof(typename Traits::value_type) : 0 >()
              , N0, N1, N2, N3, N4, N5, N6, N7 )
    {}

  template< bool AllowPadding >
  KOKKOS_INLINE_FUNCTION
  ViewMapping( void * ptr
             , const std::integral_constant<bool,AllowPadding> &
             , const typename Traits::array_layout & layout )
    : m_handle( reinterpret_cast< handle_type >( ptr ) )
    , m_offset( layout )
    {}

  //----------------------------------------
  // If the View is to construct or destroy the elements.

  struct FunctorTagConstructScalar {};
  struct FunctorTagConstructNonScalar {};
  struct FunctorTagDestructNonScalar {};

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const FunctorTagConstructScalar & , const size_t i ) const
    { m_handle[i] = 0 ; }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const FunctorTagConstructNonScalar & , const size_t i ) const
    { 
      typedef typename Traits::value_type  value_type ;
      new( & m_handle[i] ) value_type();
    }

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const FunctorTagDestructNonScalar & , const size_t i ) const
    { 
      typedef typename Traits::value_type  value_type ;
      ( & (m_handle[i]) )->~value_type();
    }

  template< class ExecSpace >
  typename std::enable_if< Kokkos::Impl::is_execution_space<ExecSpace>::value &&
                           std::is_scalar< typename Traits::value_type >::value >::type
  construct( const ExecSpace & space ) const
    {
      typedef Kokkos::RangePolicy< ExecSpace , FunctorTagConstructScalar , size_t > Policy ;

      (void) Kokkos::Impl::ParallelFor< ViewMapping , Policy >( *this , Policy( 0 , m_offset.span() ) );
      ExecSpace::fence();
    }

  template< class ExecSpace >
  typename std::enable_if< Kokkos::Impl::is_execution_space<ExecSpace>::value &&
                           ! std::is_scalar< typename Traits::value_type >::value >::type
  construct( const ExecSpace & space ) const
    {
      typedef Kokkos::RangePolicy< ExecSpace , FunctorTagConstructNonScalar , size_t > Policy ;

      (void) Kokkos::Impl::ParallelFor< ViewMapping , Policy >( *this , Policy( 0 , m_offset.span() ) );
      ExecSpace::fence();
    }

  template< class ExecSpace >
  typename std::enable_if< Kokkos::Impl::is_execution_space<ExecSpace>::value &&
                           std::is_scalar< typename Traits::value_type >::value >::type
  destroy( const ExecSpace & ) const {}

  template< class ExecSpace >
  typename std::enable_if< Kokkos::Impl::is_execution_space<ExecSpace>::value &&
                           ! std::is_scalar< typename Traits::value_type >::value >::type
  destroy( const ExecSpace & space ) const
    {
      typedef Kokkos::RangePolicy< ExecSpace , FunctorTagDestructNonScalar , size_t > Policy ;

      (void) Kokkos::Impl::ParallelFor< ViewMapping , Policy >( *this , Policy( 0 , m_offset.span() ) );
      ExecSpace::fence();
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Assign compatible default mappings */

template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space , typename SrcTraits::memory_space >::value
    &&
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    (
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename DstTraits::array_layout , Kokkos::LayoutStride >::value
    )
    &&
    std::is_same< typename SrcTraits::specialize , void >::value
    &&
    (
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutStride >::value
    )
  )>::type >
{
public:

  enum { is_assignable = true };

  typedef Kokkos::Experimental::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , void , void >  DstType ;
  typedef ViewMapping< SrcTraits , void , void >  SrcType ;

  KOKKOS_INLINE_FUNCTION
  static void assign( DstType & dst , const SrcType & src , const TrackType & src_track )
    {
      static_assert( std::is_same< typename DstTraits::value_type , typename SrcTraits::value_type >::value ||
                     std::is_same< typename DstTraits::value_type , typename SrcTraits::const_value_type >::value
                   , "View assignment must have same value type or const = non-const" );

      static_assert( ViewDimensionAssignable< typename DstTraits::dimension , typename SrcTraits::dimension >::value
                   , "View assignment must have compatible dimensions" );

      static_assert( std::is_same< typename DstTraits::array_layout , typename SrcTraits::array_layout >::value ||
                     std::is_same< typename DstTraits::array_layout , Kokkos::LayoutStride >::value ||
                     ( DstTraits::dimension::rank == 0 ) ||
                     ( DstTraits::dimension::rank == 1 && DstTraits::dimension::rank_dynamic == 1 )
                   , "View assignment must have compatible layout or have rank <= 1" );

      typedef typename DstType::offset_type  dst_offset_type ;

      dst.m_offset = dst_offset_type( src.m_offset );
      dst.m_handle = Kokkos::Experimental::Impl::ViewDataHandle< DstTraits >::assign( src.m_handle , src_track );
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief  View mapping for non-specialized data type and standard layout */
template< class Traits , bool R0 , bool R1 , bool R2 , bool R3 , bool R4 , bool R5 , bool R6 , bool R7 >
struct SubviewMapping< Traits, R0, R1, R2, R3, R4, R5, R6, R7 ,
  typename std::enable_if<(
    std::is_same< typename Traits::specialize , void >::value
    &&
    (
      std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value ||
      std::is_same< typename Traits::array_layout , Kokkos::LayoutStride >::value
    )
  )>::type >
{
private:

  // Subview's rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7) };

  // Whether right-most rank is a range.
  enum { R0_rev = 0 == Traits::rank ? false : (
                  1 == Traits::rank ? R0 : (
                  2 == Traits::rank ? R1 : (
                  3 == Traits::rank ? R2 : (
                  4 == Traits::rank ? R3 : (
                  5 == Traits::rank ? R4 : (
                  6 == Traits::rank ? R5 : (
                  7 == Traits::rank ? R6 : R7 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1 or 2, InputLayout Left, Interval 0
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0 && std::is_same< typename Traits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
        // because single stride one or second index has a stride.
        ( rank <= 2 && R0_rev && std::is_same< typename Traits::array_layout , Kokkos::LayoutRight >::value )
      ), typename Traits::array_layout , Kokkos::LayoutStride
      >::type array_layout ;

  typedef typename Traits::value_type  value_type ;

  typedef typename std::conditional< rank == 0 , value_type ,
          typename std::conditional< rank == 1 , value_type * ,
          typename std::conditional< rank == 2 , value_type ** ,
          typename std::conditional< rank == 3 , value_type *** ,
          typename std::conditional< rank == 4 , value_type **** ,
          typename std::conditional< rank == 5 , value_type ***** ,
          typename std::conditional< rank == 6 , value_type ****** ,
          typename std::conditional< rank == 7 , value_type ******* ,
                                                 value_type ********
          >::type >::type >::type >::type >::type >::type >::type >::type
     data_type ;

public:

  typedef 
    Kokkos::Experimental::ViewTraits< data_type , array_layout
                                    , typename Traits::device_type
                                    , typename Traits::memory_traits > traits_type ;

  typedef Kokkos::Experimental::View< data_type
                                    , array_layout
                                    , typename Traits::device_type
                                    , typename Traits::memory_traits > type ;

  template< class T0 , class T1 , class T2 , class T3
          , class T4 , class T5 , class T6 , class T7 >
  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< traits_type , void , void > & dst
                    , ViewMapping< Traits , void , void > const & src
                    , T0 const & arg0
                    , T1 const & arg1
                    , T2 const & arg2
                    , T3 const & arg3
                    , T4 const & arg4
                    , T5 const & arg5
                    , T6 const & arg6
                    , T7 const & arg7
                    )
    {
      typedef ViewMapping< traits_type , void , void >  DstType ;

      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T0>  V0 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T1>  V1 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T2>  V2 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T3>  V3 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T4>  V4 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T5>  V5 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T6>  V6 ;
      typedef Kokkos::Experimental::Impl::ViewOffsetRange<T7>  V7 ;

      dst.m_offset = dst_offset_type
        ( src.m_offset
        , V0::dimension( src.m_offset.dimension_0() , arg0 )
        , V1::dimension( src.m_offset.dimension_1() , arg1 )
        , V2::dimension( src.m_offset.dimension_2() , arg2 )
        , V3::dimension( src.m_offset.dimension_3() , arg3 )
        , V4::dimension( src.m_offset.dimension_4() , arg4 )
        , V5::dimension( src.m_offset.dimension_5() , arg5 )
        , V6::dimension( src.m_offset.dimension_6() , arg6 )
        , V7::dimension( src.m_offset.dimension_7() , arg7 )
        );

      dst.m_handle = dst_handle_type( src.m_handle +
                                      src.m_offset( V0::begin( arg0 )
                                                  , V1::begin( arg1 )
                                                  , V2::begin( arg2 )
                                                  , V3::begin( arg3 )
                                                  , V4::begin( arg4 )
                                                  , V5::begin( arg5 )
                                                  , V6::begin( arg6 )
                                                  , V7::begin( arg7 )
                                                  ) );
    }
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

template< class V
        , bool R0 = false , bool R1 = false , bool R2 = false , bool R3 = false
        , bool R4 = false , bool R5 = false , bool R6 = false , bool R7 = false >
struct SubviewType ;

template< class D , class A1, class A2, class A3
        , bool R0 , bool R1 , bool R2 , bool R3
        , bool R4 , bool R5 , bool R6 , bool R7 >
struct SubviewType< Kokkos::Experimental::View< D , A1, A2, A3 > , R0 , R1 , R2 , R3 , R4 , R5 , R6 , R7 >
{
private:
  typedef Kokkos::Experimental::ViewTraits< D , A1 , A2 , A3 >  traits ;
  typedef Kokkos::Experimental::Impl::SubviewMapping< traits , R0 , R1 , R2 , R3 , R4 , R5 , R6 , R7 >  mapping ;
public:
  typedef typename mapping::type  type ;
};

}}} // namespace Kokkos::Experimental::Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

class Error_view_scalar_reference_to_non_scalar_view ;

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

#if defined( KOKKOS_EXPRESSION_CHECK )

#define KOKKOS_ASSERT_VIEW_MAPPING_ACCESS( SPACE , MAP , RANK , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , SPACE >::verify( MAP.data() ); \
  /* array bounds checking */

#else

#define KOKKOS_ASSERT_VIEW_MAPPING_ACCESS( SPACE , MAP , RANK , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , SPACE >::verify( MAP.data() )

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP */

