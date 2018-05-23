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

#ifndef KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP
#define KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP

#include <cstdio>
#include <utility>

#include <Kokkos_Core.hpp>

#include <HexElement.hpp>
#include <BoxElemPart.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

/** \brief  Map a grid onto a unit cube with smooth nonlinear grading
 *          of the map.
 */
struct MapGridUnitCube {

  const float m_a ;
  const float m_b ;
  const float m_c ;
  const size_t m_max_x ;
  const size_t m_max_y ;
  const size_t m_max_z ;

  MapGridUnitCube( const size_t grid_max_x ,
                   const size_t grid_max_y ,
                   const size_t grid_max_z ,
                   const float bubble_x ,
                   const float bubble_y ,
                   const float bubble_z )
    : m_a( bubble_x )
    , m_b( bubble_y )
    , m_c( bubble_z )
    , m_max_x( grid_max_x )
    , m_max_y( grid_max_y )
    , m_max_z( grid_max_z )
    {}

  template< typename Scalar >
  KOKKOS_INLINE_FUNCTION
  void operator()( int grid_x ,
                   int grid_y ,
                   int grid_z ,
                   Scalar & coord_x ,
                   Scalar & coord_y ,
                   Scalar & coord_z ) const
    {
      // Map to a unit cube [0,1]^3

      const double x = double(grid_x) / double(m_max_x);
      const double y = double(grid_y) / double(m_max_y);
      const double z = double(grid_z) / double(m_max_z);

      coord_x = x + x * x * ( x - 1 ) * ( x - 1 ) * m_a ;
      coord_y = y + y * y * ( y - 1 ) * ( y - 1 ) * m_b ;
      coord_z = z + z * z * ( z - 1 ) * ( z - 1 ) * m_c ;
    }
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

/** \brief  Generate a distributed unstructured finite element mesh
 *          from a partitioned NX*NY*NZ box of elements.
 *
 *  Order owned nodes first followed by off-process nodes
 *  grouped by owning process.
 */
template< class Device ,
          BoxElemPart::ElemOrder Order ,
          class CoordinateMap = MapGridUnitCube >
class BoxElemFixture {
public:

  typedef Device execution_space ;

  enum { SpaceDim = 3 };
  enum { ElemNode = Order == BoxElemPart::ElemLinear ? 8 :
                    Order == BoxElemPart::ElemQuadratic ? 27 : 0 };

private:

  typedef Kokkos::Example::HexElement_TensorData< ElemNode > hex_data ;

  Kokkos::Example::BoxElemPart m_box_part ;
  CoordinateMap                m_coord_map ;

  Kokkos::View< double *[SpaceDim] , Device > m_node_coord ;
  Kokkos::View< size_t *[SpaceDim] , Device > m_node_grid ;
  Kokkos::View< size_t *[ElemNode] , Device > m_elem_node ;
  Kokkos::View< size_t *[2] ,        Device > m_recv_node ;
  Kokkos::View< size_t *[2] ,        Device > m_send_node ;
  Kokkos::View< size_t * ,           Device > m_send_node_id ;

  unsigned char m_elem_node_local[ ElemNode ][4] ;

public:

  typedef Kokkos::View< const size_t  * [ElemNode], Device > elem_node_type ;
  typedef Kokkos::View< const double  * [SpaceDim], Device > node_coord_type ;
  typedef Kokkos::View< const size_t  * [SpaceDim], Device > node_grid_type ;
  typedef Kokkos::View< const size_t  * [2] , Device > comm_list_type ;
  typedef Kokkos::View< const size_t  *     , Device > send_nodeid_type ;

  inline bool ok() const { return m_box_part.ok(); }

  KOKKOS_INLINE_FUNCTION
  size_t node_count() const { return m_node_grid.extent(0); }

  KOKKOS_INLINE_FUNCTION
  size_t node_count_owned() const { return m_box_part.owns_node_count(); }

  KOKKOS_INLINE_FUNCTION
  size_t node_count_global() const { return m_box_part.global_node_count(); }

  KOKKOS_INLINE_FUNCTION
  size_t elem_count() const { return m_elem_node.extent(0); }

  KOKKOS_INLINE_FUNCTION
  size_t elem_count_global() const { return m_box_part.global_elem_count(); }

  KOKKOS_INLINE_FUNCTION
  size_t elem_node_local( size_t inode , int k ) const
    { return m_elem_node_local[inode][k] ; }

  KOKKOS_INLINE_FUNCTION
  size_t node_grid( size_t inode , int iaxis ) const
    { return m_node_grid(inode,iaxis); }

  KOKKOS_INLINE_FUNCTION
  size_t node_global_index( size_t local ) const
    {
      const size_t tmp_node_grid[SpaceDim] =
        { m_node_grid(local,0) , m_node_grid(local,1) , m_node_grid(local,2) };
      return m_box_part.global_node_id( tmp_node_grid );
    }

  KOKKOS_INLINE_FUNCTION
  double node_coord( size_t inode , int iaxis ) const
    { return m_node_coord(inode,iaxis); }

  KOKKOS_INLINE_FUNCTION
  size_t node_grid_max( int iaxis ) const
    { return m_box_part.global_coord_max(iaxis); }

  KOKKOS_INLINE_FUNCTION
  size_t elem_node( size_t ielem , size_t inode ) const
    { return m_elem_node(ielem,inode); }

  elem_node_type   elem_node()   const { return m_elem_node ; }
  node_coord_type  node_coord()  const { return m_node_coord ; }
  node_grid_type   node_grid()   const { return m_node_grid ; }
  comm_list_type   recv_node()   const { return m_recv_node ; }
  comm_list_type   send_node()   const { return m_send_node ; }
  send_nodeid_type send_nodeid() const { return m_send_node_id ; }

  KOKKOS_INLINE_FUNCTION
  BoxElemFixture( const BoxElemFixture & rhs )
    : m_box_part(   rhs.m_box_part )
    , m_coord_map(  rhs.m_coord_map )
    , m_node_coord( rhs.m_node_coord )
    , m_node_grid(  rhs.m_node_grid )
    , m_elem_node(  rhs.m_elem_node )
    , m_recv_node(  rhs.m_recv_node )
    , m_send_node(  rhs.m_send_node )
    , m_send_node_id( rhs.m_send_node_id )
    {
      for ( int i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = rhs.m_elem_node_local[i][0] ;
        m_elem_node_local[i][1] = rhs.m_elem_node_local[i][1] ;
        m_elem_node_local[i][2] = rhs.m_elem_node_local[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
    }

  BoxElemFixture & operator = ( const BoxElemFixture & rhs )
    {
      m_box_part      = rhs.m_box_part ;
      m_coord_map     = rhs.m_coord_map ;
      m_node_coord    = rhs.m_node_coord ;
      m_node_grid     = rhs.m_node_grid ;
      m_elem_node     = rhs.m_elem_node ;
      m_recv_node     = rhs.m_recv_node ;
      m_send_node     = rhs.m_send_node ;
      m_send_node_id  = rhs.m_send_node_id ;

      for ( int i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = rhs.m_elem_node_local[i][0] ;
        m_elem_node_local[i][1] = rhs.m_elem_node_local[i][1] ;
        m_elem_node_local[i][2] = rhs.m_elem_node_local[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
      return *this ;
    }

  BoxElemFixture( const BoxElemPart::Decompose decompose ,
                  const size_t global_size ,
                  const size_t global_rank ,
                  const size_t elem_nx ,
                  const size_t elem_ny ,
                  const size_t elem_nz ,
                  const float bubble_x = 1.1f ,
                  const float bubble_y = 1.2f ,
                  const float bubble_z = 1.3f )
  : m_box_part( Order , decompose , global_size , global_rank , elem_nx , elem_ny , elem_nz )
  , m_coord_map( m_box_part.global_coord_max(0) ,
                 m_box_part.global_coord_max(1) ,
                 m_box_part.global_coord_max(2) ,
                 bubble_x ,
                 bubble_y ,
                 bubble_z )
  , m_node_coord( "fixture_node_coord" , m_box_part.uses_node_count() )
  , m_node_grid(  "fixture_node_grid" , m_box_part.uses_node_count() )
  , m_elem_node(  "fixture_elem_node" , m_box_part.uses_elem_count() )
  , m_recv_node(  "fixture_recv_node" , m_box_part.recv_node_msg_count() )
  , m_send_node(  "fixture_send_node" , m_box_part.send_node_msg_count() )
  , m_send_node_id( "fixture_send_node_id" , m_box_part.send_node_id_count() )
  {
    {
      const hex_data elem_data ;

      for ( int i = 0 ; i < ElemNode ; ++i ) {
        m_elem_node_local[i][0] = elem_data.eval_map[i][0] ;
        m_elem_node_local[i][1] = elem_data.eval_map[i][1] ;
        m_elem_node_local[i][2] = elem_data.eval_map[i][2] ;
        m_elem_node_local[i][3] = 0 ;
      }
    }

    const size_t nwork =
      std::max( m_recv_node.extent(0) ,
      std::max( m_send_node.extent(0) ,
      std::max( m_send_node_id.extent(0) ,
      std::max( m_node_grid.extent(0) ,
                m_elem_node.extent(0) * m_elem_node.extent(1) ))));

    Kokkos::parallel_for( nwork , *this );
  }


  // Initialization:

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t i ) const
  {
    if ( i < m_elem_node.extent(0) * m_elem_node.extent(1) ) {

      const size_t ielem = i / ElemNode ;
      const size_t inode = i % ElemNode ;

      size_t elem_grid[SpaceDim] ;
      size_t tmp_node_grid[SpaceDim] ;

      m_box_part.uses_elem_coord( ielem , elem_grid );

      enum { elem_node_scale = Order == BoxElemPart::ElemLinear ? 1 :
                               Order == BoxElemPart::ElemQuadratic ? 2 : 0 };

      tmp_node_grid[0] = elem_node_scale * elem_grid[0] + m_elem_node_local[inode][0] ;
      tmp_node_grid[1] = elem_node_scale * elem_grid[1] + m_elem_node_local[inode][1] ;
      tmp_node_grid[2] = elem_node_scale * elem_grid[2] + m_elem_node_local[inode][2] ;

      m_elem_node(ielem,inode) = m_box_part.local_node_id( tmp_node_grid );
    }

    if ( i < m_node_grid.extent(0) ) {
      size_t tmp_node_grid[SpaceDim] ;
      m_box_part.local_node_coord( i , tmp_node_grid );
      m_node_grid(i,0) = tmp_node_grid[0] ;
      m_node_grid(i,1) = tmp_node_grid[1] ;
      m_node_grid(i,2) = tmp_node_grid[2] ;

      m_coord_map( tmp_node_grid[0] ,
                   tmp_node_grid[1] ,
                   tmp_node_grid[2] ,
                   m_node_coord(i,0) ,
                   m_node_coord(i,1) ,
                   m_node_coord(i,2) );
    }

    if ( i < m_recv_node.extent(0) ) {
      m_recv_node(i,0) = m_box_part.recv_node_rank(i);
      m_recv_node(i,1) = m_box_part.recv_node_count(i);
    }

    if ( i < m_send_node.extent(0) ) {
      m_send_node(i,0) = m_box_part.send_node_rank(i);
      m_send_node(i,1) = m_box_part.send_node_count(i);
    }

    if ( i < m_send_node_id.extent(0) ) {
      m_send_node_id(i) = m_box_part.send_node_id(i);
    }
  }
};

} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_BOXELEMFIXTURE_HPP */

