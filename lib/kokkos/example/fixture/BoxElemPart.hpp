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

#ifndef KOKKOS_BOXELEMPART_HPP
#define KOKKOS_BOXELEMPART_HPP

#include <utility>
#include <ostream>
#include <Kokkos_Macros.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

KOKKOS_INLINE_FUNCTION
void box_intersect( size_t box[][2] ,
                    const size_t boxA[][2] ,
                    const size_t boxB[][2] )
{
  for ( int i = 0 ; i < 3 ; ++i ) {
    box[i][0] = boxA[i][0] > boxB[i][0] ? boxA[i][0] : boxB[i][0] ;
    box[i][1] = boxA[i][1] < boxB[i][1] ? boxA[i][1] : boxB[i][1] ;
    if ( box[i][0] > box[i][1] ) box[i][1] = box[i][0] ;
  }
}

KOKKOS_INLINE_FUNCTION
size_t box_count( const size_t box[][2] )
{
  return size_t( box[0][1] - box[0][0] ) *
         size_t( box[1][1] - box[1][0] ) *
         size_t( box[2][1] - box[2][0] );
}

KOKKOS_INLINE_FUNCTION
void box_ghost_layer( const size_t global_box[][2] ,
                      const size_t local_box[][2] ,
                      const size_t ghost_layer ,
                            size_t ghost_box[][2] )
{
  for ( int i = 0 ; i < 3 ; ++i ) {
    ghost_box[i][0] = global_box[i][0] + ghost_layer > local_box[i][0] ? global_box[i][0] : local_box[i][0] - ghost_layer ;
    ghost_box[i][1] = global_box[i][1] < local_box[i][1] + ghost_layer ? global_box[i][1] : local_box[i][1] + ghost_layer ;
  }
}

void box_partition( const size_t global_size ,
                    const size_t global_rank ,
                    const size_t global_box[][2] ,
                          size_t box[][2] );

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

/** \brief Partition a box of hexahedral elements among subdomains.
 *
 *  Nodes are ordered locally as follows:
 *    { owned_by[ this_process ] ,
 *      owned_by[ neighbor_process[0] ] ,
 *      owned_by[ neighbor_process[1] ] ,
 *      owned_by[ neighbor_process[2] ] ,
 *      ... };
 */
class BoxElemPart {
public:

  enum Decompose { DecomposeNode , DecomposeElem };
  enum ElemOrder { ElemLinear , ElemQuadratic };

  bool ok() const { return m_ok ; }

  BoxElemPart( const ElemOrder elem_order ,
               const Decompose decompose ,
               const size_t global_size ,
               const size_t global_rank ,
               const size_t elem_nx ,
               const size_t elem_ny ,
               const size_t elem_nz );

  KOKKOS_INLINE_FUNCTION
  size_t global_elem_count() const
    { return Kokkos::Example::box_count( m_global_elem_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t global_node_count() const
    { return Kokkos::Example::box_count( m_global_node_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t uses_elem_count() const
    { return Kokkos::Example::box_count( m_uses_elem_box ); }

  KOKKOS_INLINE_FUNCTION
  size_t owns_node_count() const
    { return Kokkos::Example::box_count( m_owns_node_box[0] ); }

  KOKKOS_INLINE_FUNCTION
  size_t uses_node_count() const
    { return Kokkos::Example::box_count( m_uses_node_box ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  size_t uses_elem_offset( const size_t ix ,
                           const size_t iy ,
                           const size_t iz ) const
  {
    return size_t( ix - m_uses_elem_box[0][0] ) + size_t( m_uses_elem_box[0][1] - m_uses_elem_box[0][0] ) * (
           size_t( iy - m_uses_elem_box[1][0] ) + size_t( m_uses_elem_box[1][1] - m_uses_elem_box[1][0] ) * (
           size_t( iz - m_uses_elem_box[2][0] ) ) );
  }

  KOKKOS_INLINE_FUNCTION
  void uses_elem_coord( size_t lid , size_t c[] ) const
  {
    const size_t nx = m_uses_elem_box[0][1] - m_uses_elem_box[0][0] ;
    const size_t ny = m_uses_elem_box[1][1] - m_uses_elem_box[1][0] ;

    c[0] = m_uses_elem_box[0][0] + lid % nx ; lid /= nx ;
    c[1] = m_uses_elem_box[1][0] + lid % ny ; lid /= ny ;
    c[2] = m_uses_elem_box[2][0] + lid ;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  size_t global_coord_max( size_t axis ) const
  { return m_global_node_box[axis][1] - 1 ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  void local_node_coord( size_t lid , size_t coord[] ) const
  {
    // Local id within an 'owns' block (has sentinal)
    size_t j = 0 ;
    while ( m_owns_node[j][1] <= lid ) { lid -= m_owns_node[j][1] ; ++j ; }

    // Map to global coordinates:
    const size_t nx = m_owns_node_box[j][0][1] - m_owns_node_box[j][0][0] ;
    const size_t ny = m_owns_node_box[j][1][1] - m_owns_node_box[j][1][0] ;

    coord[0] = m_owns_node_box[j][0][0] + lid % nx ; lid /= nx ;
    coord[1] = m_owns_node_box[j][1][0] + lid % ny ; lid /= ny ;
    coord[2] = m_owns_node_box[j][2][0] + lid ;
  }

  KOKKOS_INLINE_FUNCTION
  size_t local_node_id( const size_t c[] ) const
  {
    // Find which 'owns' block and accumulate the offset of this block:
    size_t lid = 0 ;
    size_t j = 0 ;
    while ( ! ( m_owns_node_box[j][0][0] <= c[0] && c[0] < m_owns_node_box[j][0][1] &&
                m_owns_node_box[j][1][0] <= c[1] && c[1] < m_owns_node_box[j][1][1] &&
                m_owns_node_box[j][2][0] <= c[2] && c[2] < m_owns_node_box[j][2][1] ) ) {
      
      lid += m_owns_node[j][1] ;
      ++j ;
    }

    // Map offset to the block plus offset within the block:
    return lid +
           size_t( c[0] - m_owns_node_box[j][0][0] ) + size_t( m_owns_node_box[j][0][1] - m_owns_node_box[j][0][0] ) * (
           size_t( c[1] - m_owns_node_box[j][1][0] ) + size_t( m_owns_node_box[j][1][1] - m_owns_node_box[j][1][0] ) * (
           size_t( c[2] - m_owns_node_box[j][2][0] ) ) );
  }

  KOKKOS_INLINE_FUNCTION
  size_t global_node_id( const size_t c[] ) const
  {
    return size_t( c[0] - m_global_node_box[0][0] ) + size_t( m_global_node_box[0][1] - m_global_node_box[0][0] ) * (
           size_t( c[1] - m_global_node_box[1][0] ) + size_t( m_global_node_box[1][1] - m_global_node_box[1][0] ) * (
           size_t( c[2] - m_global_node_box[2][0] ) ) );
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  size_t recv_node_msg_count() const { return m_owns_node_count - 1 ; }

  KOKKOS_INLINE_FUNCTION
  size_t recv_node_rank(  size_t msg ) const { return m_owns_node[msg+1][0] ; }

  KOKKOS_INLINE_FUNCTION
  size_t recv_node_count( size_t msg ) const { return m_owns_node[msg+1][1] ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  size_t send_node_msg_count() const { return m_send_node_count ; }

  KOKKOS_INLINE_FUNCTION
  size_t send_node_rank(  size_t msg ) const { return m_send_node[msg][0] ; }

  KOKKOS_INLINE_FUNCTION
  size_t send_node_count( size_t msg ) const { return m_send_node[msg][1] ; }

  KOKKOS_INLINE_FUNCTION
  size_t send_node_id_count() const
  {
    size_t count = 0 ;
    for ( size_t i = 0 ; i < m_send_node_count ; ++i ) {
      count += m_send_node[i][1] ;
    }
    return count ;
  }

  KOKKOS_INLINE_FUNCTION
  size_t send_node_id( size_t item ) const
  {
    // Find which send list this send item is in:
    size_t j = 0 ;
    while ( m_send_node[j][1] <= item ) { item -= m_send_node[j][1] ; ++j ; }

    // Map to global coordinate:
    const size_t nx = m_send_node_box[j][0][1] - m_send_node_box[j][0][0] ;
    const size_t ny = m_send_node_box[j][1][1] - m_send_node_box[j][1][0] ;

    size_t c[3] ;

    c[0] = m_send_node_box[j][0][0] + item % nx ; item /= nx ;
    c[1] = m_send_node_box[j][1][0] + item % ny ; item /= ny ;
    c[2] = m_send_node_box[j][2][0] + item ;

    // Map to local id:
    return size_t( c[0] - m_owns_node_box[0][0][0] ) + size_t( m_owns_node_box[0][0][1] - m_owns_node_box[0][0][0] ) * (
           size_t( c[1] - m_owns_node_box[0][1][0] ) + size_t( m_owns_node_box[0][1][1] - m_owns_node_box[0][1][0] ) * (
           size_t( c[2] - m_owns_node_box[0][2][0] ) ) );
  }

  //----------------------------------------

  void print( std::ostream & s ) const ;

private:

  // Maximum number of processes in a neighborhood, including this process
  enum { PROC_NEIGH_MAX = 64 };

  void local( const size_t  rank ,
                    size_t  uses_elem[][2] ,
                    size_t  owns_node[][2] ,
                    size_t  uses_node[][2] ) const ;

  size_t  m_global_size ;
  size_t  m_global_rank ;

  Decompose m_decompose ;
  ElemOrder m_elem_order ;

  size_t m_global_elem_box[3][2] ;
  size_t m_global_node_box[3][2] ;
  size_t m_uses_elem_box[3][2] ;
  size_t m_uses_node_box[3][2] ;

  // [ processor rank , count ]
  size_t m_owns_node_box[ PROC_NEIGH_MAX ][3][2] ;
  size_t m_owns_node[     PROC_NEIGH_MAX ][2] ;
  size_t m_owns_node_count ;

  size_t m_send_node_box[ PROC_NEIGH_MAX ][3][2] ;
  size_t m_send_node[     PROC_NEIGH_MAX ][2] ;
  size_t m_send_node_count ;

  bool   m_ok ;
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_BOXELEMPART_HPP */

