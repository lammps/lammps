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

#include <utility>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <BoxElemPart.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

void box_partition( const size_t global_size ,
                    const size_t global_rank ,
                    const size_t global_box[][2] ,
                          size_t box[][2] )
{
  box[0][0] = global_box[0][0] ; box[0][1] = global_box[0][1] ;
  box[1][0] = global_box[1][0] ; box[1][1] = global_box[1][1] ;
  box[2][0] = global_box[2][0] ; box[2][1] = global_box[2][1] ;

  size_t ip = 0 ;
  size_t np = global_size ;

  while ( 1 < np ) {

    // P = [ ip + j * portion , ip + ( j + 1 ) * portion )

    size_t jip , jup ;

    {
      const size_t part = ( 0 == ( np % 5 ) ) ? 5 : (
                          ( 0 == ( np % 3 ) ) ? 3 : 2 );

      const size_t portion = np / part ;

      if ( 2 < part || global_rank < ip + portion ) {
        jip = portion * size_t( double( global_rank - ip ) / double(portion) );
        jup = jip + portion ;
      }
      else {
        jip = portion ;
        jup = np ;
      }
    }

    // Choose axis with largest count:

    const size_t nb[3] = {
      box[0][1] - box[0][0] ,
      box[1][1] - box[1][0] ,
      box[2][1] - box[2][0] };

    const int axis = nb[2] > nb[1] ? ( nb[2] > nb[0] ? 2 : 0 )
                                        : ( nb[1] > nb[0] ? 1 : 0 );

    box[ axis ][1] = box[ axis ][0] + size_t( double(nb[axis]) * ( double(jup) / double(np) ));
    box[ axis ][0] = box[ axis ][0] + size_t( double(nb[axis]) * ( double(jip) / double(np) ));

    np = jup - jip ;
    ip = ip + jip ;
  }
}

} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

void BoxElemPart::local( const size_t  rank ,
                               size_t  uses_elem[][2] ,
                               size_t  owns_node[][2] ,
                               size_t  uses_node[][2] ) const
{
  if ( BoxElemPart::DecomposeElem == m_decompose ) {

    Kokkos::Example::box_partition( m_global_size , rank , m_global_elem_box , uses_elem );

    for ( int i = 0 ; i < 3 ; ++i ) {
      owns_node[i][0] = uses_elem[i][0] ;
      owns_node[i][1] = uses_elem[i][1] + ( m_global_elem_box[i][1] == uses_elem[i][1] ? 1 : 0 );
    }
  }
  else {

    const size_t global_vert[3][2] =
      { { 0 , m_global_elem_box[0][1] + 1 },
        { 0 , m_global_elem_box[1][1] + 1 },
        { 0 , m_global_elem_box[2][1] + 1 } };

    Kokkos::Example::box_partition( m_global_size , rank , global_vert , owns_node );

    for ( int i = 0 ; i < 3 ; ++i ) {
      uses_elem[i][0] = global_vert[i][0] == owns_node[i][0] ? owns_node[i][0] : owns_node[i][0] - 1 ;
      uses_elem[i][1] = global_vert[i][1] == owns_node[i][1] ? owns_node[i][1] - 1 : owns_node[i][1] ;
    }
  }

  for ( int i = 0 ; i < 3 ; ++i ) {
    uses_node[i][0] = uses_elem[i][0] ;
    uses_node[i][1] = uses_elem[i][1] + 1 ;
  }

  if ( BoxElemPart::ElemQuadratic == m_elem_order ) {
    for ( int i = 0 ; i < 3 ; ++i ) {
      owns_node[i][0] = 2 * owns_node[i][0] ;
      uses_node[i][0] = 2 * uses_node[i][0] ;
      owns_node[i][1] = 2 * owns_node[i][1] - 1 ;
      uses_node[i][1] = 2 * uses_node[i][1] - 1 ;
    }
  }
}

BoxElemPart::BoxElemPart(
  const BoxElemPart::ElemOrder elem_order ,
  const BoxElemPart::Decompose decompose ,
  const size_t global_size ,
  const size_t global_rank ,
  const size_t elem_nx ,
  const size_t elem_ny ,
  const size_t elem_nz )
{
  m_global_size = global_size ;
  m_global_rank = global_rank ;
  m_decompose   = decompose ;
  m_elem_order  = elem_order ;

  m_global_elem_box[0][0] = 0 ; m_global_elem_box[0][1] = elem_nx ;
  m_global_elem_box[1][0] = 0 ; m_global_elem_box[1][1] = elem_ny ;
  m_global_elem_box[2][0] = 0 ; m_global_elem_box[2][1] = elem_nz ;

  m_global_node_box[0][0] = 0 ; m_global_node_box[0][1] = 0 ;
  m_global_node_box[1][0] = 0 ; m_global_node_box[1][1] = 0 ;
  m_global_node_box[2][0] = 0 ; m_global_node_box[2][1] = 0 ;

  m_owns_node_count = 0 ;
  m_send_node_count = 0 ;

  m_ok = true ;

  //----------------------------------------

  if ( ElemLinear == elem_order ) {
    m_global_node_box[0][1] = elem_nx + 1 ;
    m_global_node_box[1][1] = elem_ny + 1 ;
    m_global_node_box[2][1] = elem_nz + 1 ;
  }
  else if ( ElemQuadratic == elem_order ) {
    m_global_node_box[0][1] = 2 * elem_nx + 1 ;
    m_global_node_box[1][1] = 2 * elem_ny + 1 ;
    m_global_node_box[2][1] = 2 * elem_nz + 1 ;
  }

  //----------------------------------------

  local( m_global_rank , m_uses_elem_box , m_owns_node_box[0] , m_uses_node_box );

  const size_t global_node_count_ = Kokkos::Example::box_count( m_global_node_box );
  const size_t global_elem_count_ = Kokkos::Example::box_count( m_global_elem_box );

  //----------------------------------------

  size_t elem_count = Kokkos::Example::box_count( m_uses_elem_box );
  size_t node_count = Kokkos::Example::box_count( m_owns_node_box[0] );

  m_owns_node[0][0] = global_rank ;
  m_owns_node[0][1] = node_count ;
  m_owns_node_count = 1 ;
  m_send_node_count = 0 ;

  for ( size_t rr = 1 ; rr < m_global_size && m_ok ; ++rr ) {

    const size_t rank = ( m_global_rank + rr ) % m_global_size ;

    size_t elem_box[3][2] , o_node_box[3][2] , u_node_box[3][2] ;

    // Boxes for process 'rank'
    local( rank , elem_box , o_node_box , u_node_box );

    // Box that this process uses but is owned by process 'rank'
    Kokkos::Example::box_intersect( m_owns_node_box[ m_owns_node_count ] , m_uses_node_box , o_node_box );

    m_owns_node[ m_owns_node_count ][1] = Kokkos::Example::box_count( m_owns_node_box[ m_owns_node_count ] );

    if ( m_owns_node[ m_owns_node_count ][1] ) {

      if ( ( PROC_NEIGH_MAX - 1 ) <= m_owns_node_count ) {
        std::cout << "BoxElemPart exceeded maximum neighbor count" << std::endl ;
        m_ok = false ;
        break ;
      }

      m_owns_node[ m_owns_node_count ][0] = rank ;

      ++m_owns_node_count ;
    }

    // Box that this process owns and is used by process 'rank'
    Kokkos::Example::box_intersect( m_send_node_box[ m_send_node_count ] , m_owns_node_box[0] , u_node_box );

    m_send_node[ m_send_node_count ][1] = Kokkos::Example::box_count( m_send_node_box[ m_send_node_count ] );

    if ( m_send_node[ m_send_node_count ][1] ) {

      if ( ( PROC_NEIGH_MAX - 1 ) <= m_send_node_count ) {
        std::cout << "BoxElemPart exceeded maximum neighbor count" << std::endl ;
        m_ok = false ;
        break ;
      }

      m_send_node[ m_send_node_count ][0] = rank ;
      ++m_send_node_count ;
    }

    // Error checking:

    size_t test_box[3][2] ;

    elem_count += Kokkos::Example::box_count( elem_box );
    node_count += Kokkos::Example::box_count( o_node_box );

    {
      Kokkos::Example::box_intersect( test_box , m_owns_node_box[0] , o_node_box );

      if ( Kokkos::Example::box_count( test_box ) ) {
        std::cout << "Box partitioning error" << std::endl ;
        std::cout << "owns_node[" << m_global_rank << "]{"
                  << " [" << m_owns_node_box[0][0][0] << "," << m_owns_node_box[0][0][1] << ")"
                  << " [" << m_owns_node_box[0][1][0] << "," << m_owns_node_box[0][1][1] << ")"
                  << " [" << m_owns_node_box[0][2][0] << "," << m_owns_node_box[0][2][1] << ")"
                  << "} intersects"
                  << " owns_node[" << rank << "]{"
                  << " [" << o_node_box[0][0] << "," << o_node_box[0][1] << ")"
                  << " [" << o_node_box[1][0] << "," << o_node_box[1][1] << ")"
                  << " [" << o_node_box[2][0] << "," << o_node_box[2][1] << ")"
                  << "}" << std::endl ;
        m_ok = false ;
        break ;
      }
    }

    if ( DecomposeElem == decompose ) {

      Kokkos::Example::box_intersect( test_box , m_uses_elem_box , elem_box );

      if ( Kokkos::Example::box_count( test_box ) ) {
        std::cout << "Box partitioning error" << std::endl ;
        std::cout << "ElemBox[" << m_global_rank << "]{"
                  << " [" << m_uses_elem_box[0][0] << "," << m_uses_elem_box[0][1] << ")"
                  << " [" << m_uses_elem_box[1][0] << "," << m_uses_elem_box[1][1] << ")"
                  << " [" << m_uses_elem_box[2][0] << "," << m_uses_elem_box[2][1] << ")"
                  << "} intersects"
                  << " ElemBox[" << rank << "]{"
                  << " [" << elem_box[0][0] << "," << elem_box[0][1] << ")"
                  << " [" << elem_box[1][0] << "," << elem_box[1][1] << ")"
                  << " [" << elem_box[2][0] << "," << elem_box[2][1] << ")"
                  << "}" << std::endl ;
        m_ok = false ;
        break ;
      }
    }
  }

  // Sentinal values at the end of the owns and send lists:

  m_owns_node[ m_owns_node_count ][0] = ~0u ;
  m_owns_node[ m_owns_node_count ][1] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][0][0] = 0u ; m_owns_node_box[ m_owns_node_count ][0][0] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][1][0] = 0u ; m_owns_node_box[ m_owns_node_count ][1][0] = ~0u ;
  m_owns_node_box[ m_owns_node_count ][2][0] = 0u ; m_owns_node_box[ m_owns_node_count ][2][0] = ~0u ;

  m_send_node[ m_send_node_count ][0] = ~0u ;
  m_send_node[ m_send_node_count ][1] = ~0u ;
  m_send_node_box[ m_send_node_count ][0][0] = 0u ; m_send_node_box[ m_send_node_count ][0][0] = ~0u ;
  m_send_node_box[ m_send_node_count ][1][0] = 0u ; m_send_node_box[ m_send_node_count ][1][0] = ~0u ;
  m_send_node_box[ m_send_node_count ][2][0] = 0u ; m_send_node_box[ m_send_node_count ][2][0] = ~0u ;

  {
    size_t count = 0 ;
    for ( size_t i = 0 ; i < m_owns_node_count ; ++i ) {
      count += m_owns_node[i][1] ;
    }
    if ( count != Kokkos::Example::box_count( m_uses_node_box ) ) {
      std::cout << "Node uses count = " << Kokkos::Example::box_count( m_uses_node_box )
                << " error count = " << count << std::endl ;
      m_ok = false ;
    }
  }

  if ( global_node_count_ != node_count ) {
    std::cout << "Node count = " << global_node_count_ << " overlap error count = " << node_count << std::endl ;
    m_ok = false ;
  }

  if ( DecomposeElem == decompose && global_elem_count_ != elem_count ) {
    std::cout << "Elem count = " << global_elem_count_ << " overlap error count = " << elem_count << std::endl ;
    m_ok = false ;
  }

  if ( ! m_ok ) {
    for ( int i = 0 ; i < 3 ; ++i ) { for ( int j = 0 ; j < 2 ; ++j ) {
      m_global_elem_box[i][j] = 0 ;
      m_global_node_box[i][j] = 0 ;
      m_uses_elem_box[i][j] = 0 ;
      m_uses_node_box[i][j] = 0 ;
    }}
    m_owns_node_count = 0 ;
    m_send_node_count = 0 ;
  }
}

void BoxElemPart::print( std::ostream & s ) const
{
  s << "BoxElemPart P[" << m_global_rank << ":" << m_global_size << "]"
    << std::endl
    << "  elem_box {"
    << " [" << m_uses_elem_box[0][0] << "," << m_uses_elem_box[0][1] << ")"
    << " [" << m_uses_elem_box[1][0] << "," << m_uses_elem_box[1][1] << ")"
    << " [" << m_uses_elem_box[2][0] << "," << m_uses_elem_box[2][1] << ")"
    << " } / {"
    << " [" << m_global_elem_box[0][0] << "," << m_global_elem_box[0][1] << ")"
    << " [" << m_global_elem_box[1][0] << "," << m_global_elem_box[1][1] << ")"
    << " [" << m_global_elem_box[2][0] << "," << m_global_elem_box[2][1] << ")"
    << " }"
    << std::endl
    << "  node_box {"
    << " [" << m_owns_node_box[0][0][0] << "," << m_owns_node_box[0][0][1] << ")"
    << " [" << m_owns_node_box[0][1][0] << "," << m_owns_node_box[0][1][1] << ")"
    << " [" << m_owns_node_box[0][2][0] << "," << m_owns_node_box[0][2][1] << ")"
    << " } / {"
    << " [" << m_uses_node_box[0][0] << "," << m_uses_node_box[0][1] << ")"
    << " [" << m_uses_node_box[1][0] << "," << m_uses_node_box[1][1] << ")"
    << " [" << m_uses_node_box[2][0] << "," << m_uses_node_box[2][1] << ")"
    << " } / {"
    << " [" << m_global_node_box[0][0] << "," << m_global_node_box[0][1] << ")"
    << " [" << m_global_node_box[1][0] << "," << m_global_node_box[1][1] << ")"
    << " [" << m_global_node_box[2][0] << "," << m_global_node_box[2][1] << ")"
    << " }"
    << std::endl ;

  for ( size_t i = 1 ; i < m_owns_node_count ; ++i ) {
    s << "  P[" << m_owns_node[i][0] << "]"
      << " recv node_box {"
      << " [" << m_owns_node_box[i][0][0] << "," << m_owns_node_box[i][0][1] << ")"
      << " [" << m_owns_node_box[i][1][0] << "," << m_owns_node_box[i][1][1] << ")"
      << " [" << m_owns_node_box[i][2][0] << "," << m_owns_node_box[i][2][1] << ")"
      << " }"
      << std::endl ;
  }

  for ( size_t i = 0 ; i < m_send_node_count ; ++i ) {
    s << "  P[" << m_send_node[i][0] << "]"
      << " send node_box {"
      << " [" << m_send_node_box[i][0][0] << "," << m_send_node_box[i][0][1] << ")"
      << " [" << m_send_node_box[i][1][0] << "," << m_send_node_box[i][1][1] << ")"
      << " [" << m_send_node_box[i][2][0] << "," << m_send_node_box[i][2][1] << ")"
      << " }"
      << std::endl ;
  }
}

} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------


