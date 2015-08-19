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

#ifndef KOKKOS_BOXMESHFIXTURE_HPP
#define KOKKOS_BOXMESHFIXTURE_HPP

#include <cmath>
#include <stdexcept>
#include <sstream>

#include <Kokkos_Core.hpp>
#include <BoxMeshPartition.hpp>
#include <FEMesh.hpp>
#include <HexElement.hpp>

//----------------------------------------------------------------------------

struct FixtureElementHex8 {

  static const unsigned element_node_count = 8 ;

  HybridFEM::HexElement_TensorData< element_node_count > elem_data ;
  BoxBoundsLinear box_bounds ;

  FixtureElementHex8() : elem_data(), box_bounds() {}

  static void create_node_boxes_from_vertex_boxes(
    const BoxType                & vertex_box_global ,
    const std::vector< BoxType > & vertex_box_parts ,
          BoxType                & node_box_global ,
          std::vector< BoxType > & node_box_parts )
  {
    node_box_global = vertex_box_global ;
    node_box_parts  = vertex_box_parts  ;
  }

  void elem_to_node( const unsigned node_local , unsigned coord[] ) const
  {
    coord[0] += elem_data.eval_map[ node_local ][0] ;
    coord[1] += elem_data.eval_map[ node_local ][1] ;
    coord[2] += elem_data.eval_map[ node_local ][2] ;
  }
};

struct FixtureElementHex27 {
  static const unsigned element_node_count = 27 ;

  HybridFEM::HexElement_TensorData< element_node_count > elem_data ;
  BoxBoundsQuadratic box_bounds ;

  FixtureElementHex27() : elem_data(), box_bounds() {}

  static void create_node_boxes_from_vertex_boxes(
    const BoxType                & vertex_box_global ,
    const std::vector< BoxType > & vertex_box_parts ,
          BoxType                & node_box_global ,
          std::vector< BoxType > & node_box_parts )
  {
    node_box_global = vertex_box_global ;
    node_box_parts  = vertex_box_parts  ;

    node_box_global[0][1] = 2 * node_box_global[0][1] - 1 ;
    node_box_global[1][1] = 2 * node_box_global[1][1] - 1 ;
    node_box_global[2][1] = 2 * node_box_global[2][1] - 1 ;

    for ( unsigned i = 0 ; i < vertex_box_parts.size() ; ++i ) {
      node_box_parts[i][0][0] = 2 * node_box_parts[i][0][0] ;
      node_box_parts[i][1][0] = 2 * node_box_parts[i][1][0] ;
      node_box_parts[i][2][0] = 2 * node_box_parts[i][2][0] ;

      node_box_parts[i][0][1] =
        std::min( node_box_global[0][1] , 2 * node_box_parts[i][0][1] );
      node_box_parts[i][1][1] =
        std::min( node_box_global[1][1] , 2 * node_box_parts[i][1][1] );
      node_box_parts[i][2][1] =
        std::min( node_box_global[2][1] , 2 * node_box_parts[i][2][1] );
    }
  }

  void elem_to_node( const unsigned node_local , unsigned coord[] ) const
  {
    coord[0] = 2 * coord[0] + elem_data.eval_map[ node_local ][0] ;
    coord[1] = 2 * coord[1] + elem_data.eval_map[ node_local ][1] ;
    coord[2] = 2 * coord[2] + elem_data.eval_map[ node_local ][2] ;
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , class Device , class ElementSpec >
struct BoxMeshFixture {

  typedef Scalar  coordinate_scalar_type ;
  typedef Device  execution_space ;

  static const unsigned element_node_count = ElementSpec::element_node_count ;

  typedef HybridFEM::FEMesh< coordinate_scalar_type ,
                             element_node_count ,
                             execution_space > FEMeshType ;

  typedef typename FEMeshType::node_coords_type    node_coords_type ;
  typedef typename FEMeshType::elem_node_ids_type  elem_node_ids_type ;
  typedef typename FEMeshType::node_elem_ids_type  node_elem_ids_type ;


  static void verify(
    const typename FEMeshType::node_coords_type::HostMirror   & node_coords ,
    const typename FEMeshType::elem_node_ids_type::HostMirror & elem_node_ids ,
    const typename FEMeshType::node_elem_ids_type::HostMirror & node_elem_ids )
  {
    typedef typename FEMeshType::size_type         size_type ;
    //typedef typename node_coords_type::value_type  coords_type ; // unused

    const size_type node_count_total = node_coords.dimension_0();
    const size_type elem_count_total = elem_node_ids.dimension_0();

    const ElementSpec element ;

    for ( size_type node_index = 0 ;
                    node_index < node_count_total ; ++node_index ) {

      for ( size_type
              j = node_elem_ids.row_map[ node_index ] ;
              j < node_elem_ids.row_map[ node_index + 1 ] ; ++j ) {

        const size_type elem_index = node_elem_ids.entries(j,0);
        const size_type node_local = node_elem_ids.entries(j,1);
        const size_type en_id      = elem_node_ids(elem_index,node_local);

        if ( node_index != en_id ) {
          std::ostringstream msg ;
          msg << "BoxMeshFixture node_elem_ids error"
              << " : node_index(" << node_index
              << ") entry(" << j
              << ") elem_index(" << elem_index
              << ") node_local(" << node_local
              << ") elem_node_id(" << en_id
              << ")" ;
          throw std::runtime_error( msg.str() );
        }
      }
    }

    for ( size_type elem_index = 0 ;
                    elem_index < elem_count_total; ++elem_index ) {

      coordinate_scalar_type elem_node_coord[ element_node_count ][3] ;

      for ( size_type nn = 0 ; nn < element_node_count ; ++nn ) {
        const size_type node_index = elem_node_ids( elem_index , nn );

        for ( size_type nc = 0 ; nc < 3 ; ++nc ) {
          elem_node_coord[nn][nc] = node_coords( node_index , nc );
        }
      }


      for ( size_type nn = 0 ; nn < element_node_count ; ++nn ) {

        const unsigned ix = element.elem_data.eval_map[nn][0] ;
        const unsigned iy = element.elem_data.eval_map[nn][1] ;
        const unsigned iz = element.elem_data.eval_map[nn][2] ;

        if ( elem_node_coord[nn][0] != elem_node_coord[0][0] + ix ||
             elem_node_coord[nn][1] != elem_node_coord[0][1] + iy ||
             elem_node_coord[nn][2] != elem_node_coord[0][2] + iz ) {

          std::ostringstream msg ;
          msg << "BoxMeshFixture elem_node_coord mapping failure { "
              << elem_node_coord[nn][0] << " "
              << elem_node_coord[nn][1] << " "
              << elem_node_coord[nn][2] << " } != { "
              << elem_node_coord[ 0][0] + ix << " "
              << elem_node_coord[ 0][1] + iy << " "
              << elem_node_coord[ 0][2] + iz
              << " }" ;
          throw std::runtime_error( msg.str() );
        }
      }
    }
  }

  //------------------------------------
  // Initialize element-node connectivity:
  // Order elements that only depend on owned nodes first.
  // These elements could be computed while waiting for
  // received node data.

  static void layout_elements_interior_exterior(
    const BoxType                vertex_box_local_used ,
    const BoxType                vertex_box_local_owned ,
    const BoxType                node_box_local_used ,
    const std::vector<size_t> &  node_used_id_map ,
    const ElementSpec            element_fixture ,
    const size_t                 elem_count_interior ,
    const typename elem_node_ids_type::HostMirror elem_node_ids )
  {
    size_t elem_index_interior = 0 ;
    size_t elem_index_boundary = elem_count_interior ;

    for ( size_t iz = vertex_box_local_used[2][0] ;
                 iz < vertex_box_local_used[2][1] - 1 ; ++iz ) {
    for ( size_t iy = vertex_box_local_used[1][0] ;
                 iy < vertex_box_local_used[1][1] - 1 ; ++iy ) {
    for ( size_t ix = vertex_box_local_used[0][0] ;
                 ix < vertex_box_local_used[0][1] - 1 ; ++ix ) {

      size_t elem_index ;

      // If lower and upper vertices are owned then element is interior
      if ( contain( vertex_box_local_owned, ix,   iy,   iz ) &&
           contain( vertex_box_local_owned, ix+1, iy+1, iz+1 ) ) {
        elem_index = elem_index_interior++ ;
      }
      else {
        elem_index = elem_index_boundary++ ;
      }

      for ( size_t nn = 0 ; nn < element_node_count ; ++nn ) {
        unsigned coord[3] = { static_cast<unsigned>(ix) , static_cast<unsigned>(iy) , static_cast<unsigned>(iz) };

        element_fixture.elem_to_node( nn , coord );

        const size_t node_local_id =
          box_map_id( node_box_local_used ,
                      node_used_id_map ,
                      coord[0] , coord[1] , coord[2] );

        elem_node_ids( elem_index , nn ) = node_local_id ;
      }
    }}}
  }

  //------------------------------------
  // Nested partitioning of elements by number of thread 'gangs'

  static void layout_elements_partitioned(
    const BoxType                vertex_box_local_used ,
    const BoxType                /*vertex_box_local_owned*/ ,
    const BoxType                node_box_local_used ,
    const std::vector<size_t> &  node_used_id_map ,
    const ElementSpec            element_fixture ,
    const size_t                 thread_gang_count ,
    const typename elem_node_ids_type::HostMirror elem_node_ids )
  {
    std::vector< BoxType > element_box_gangs( thread_gang_count );

    BoxType element_box_local_used = vertex_box_local_used ;

    element_box_local_used[0][1] -= 1 ;
    element_box_local_used[1][1] -= 1 ;
    element_box_local_used[2][1] -= 1 ;

    box_partition_rcb( element_box_local_used , element_box_gangs );

    size_t elem_index = 0 ;

    for ( size_t ig = 0 ; ig < thread_gang_count ; ++ig ) {

      const BoxType box = element_box_gangs[ig] ;

      for ( size_t iz = box[2][0] ; iz < box[2][1] ; ++iz ) {
      for ( size_t iy = box[1][0] ; iy < box[1][1] ; ++iy ) {
      for ( size_t ix = box[0][0] ; ix < box[0][1] ; ++ix , ++elem_index ) {

        for ( size_t nn = 0 ; nn < element_node_count ; ++nn ) {
          unsigned coord[3] = { static_cast<unsigned>(ix) , static_cast<unsigned>(iy) , static_cast<unsigned>(iz) };

          element_fixture.elem_to_node( nn , coord );

          const size_t node_local_id =
            box_map_id( node_box_local_used ,
                        node_used_id_map ,
                        coord[0] , coord[1] , coord[2] );

          elem_node_ids( elem_index , nn ) = node_local_id ;
        }
      }}}
    }
  }

  //------------------------------------

  static FEMeshType create( const size_t proc_count ,
                            const size_t proc_local ,
                            const size_t gang_count ,
                            const size_t elems_x ,
                            const size_t elems_y ,
                            const size_t elems_z ,
                            const double x_coord_curve = 1 ,
                            const double y_coord_curve = 1 ,
                            const double z_coord_curve = 1 )
  {
    const size_t vertices_x = elems_x + 1 ;
    const size_t vertices_y = elems_y + 1 ;
    const size_t vertices_z = elems_z + 1 ;

    const BoxBoundsLinear vertex_box_bounds ;
    const ElementSpec element ;

    // Partition based upon vertices:

    BoxType vertex_box_global ;
    std::vector< BoxType > vertex_box_parts( proc_count );

    vertex_box_global[0][0] = 0 ; vertex_box_global[0][1] = vertices_x ;
    vertex_box_global[1][0] = 0 ; vertex_box_global[1][1] = vertices_y ;
    vertex_box_global[2][0] = 0 ; vertex_box_global[2][1] = vertices_z ;

    box_partition_rcb( vertex_box_global , vertex_box_parts );

    const BoxType vertex_box_local_owned = vertex_box_parts[ proc_local ];

    // Determine interior and used vertices:

    BoxType vertex_box_local_interior ;
    BoxType vertex_box_local_used ;

    vertex_box_bounds.apply( vertex_box_global ,
                             vertex_box_local_owned ,
                             vertex_box_local_interior ,
                             vertex_box_local_used );

    // Element counts:

    const long local_elems_x =
      ( vertex_box_local_used[0][1] - vertex_box_local_used[0][0] ) - 1 ;
    const long local_elems_y =
      ( vertex_box_local_used[1][1] - vertex_box_local_used[1][0] ) - 1 ;
    const long local_elems_z =
      ( vertex_box_local_used[2][1] - vertex_box_local_used[2][0] ) - 1 ;

    const size_t elem_count_total = std::max( long(0) , local_elems_x ) *
                                    std::max( long(0) , local_elems_y ) *
                                    std::max( long(0) , local_elems_z );

    const long interior_elems_x =
      ( vertex_box_local_owned[0][1] - vertex_box_local_owned[0][0] ) - 1 ;
    const long interior_elems_y =
      ( vertex_box_local_owned[1][1] - vertex_box_local_owned[1][0] ) - 1 ;
    const long interior_elems_z =
      ( vertex_box_local_owned[2][1] - vertex_box_local_owned[2][0] ) - 1 ;

    const size_t elem_count_interior = std::max( long(0) , interior_elems_x ) *
                                       std::max( long(0) , interior_elems_y ) *
                                       std::max( long(0) , interior_elems_z );

    // Expand vertex boxes to node boxes:

    BoxType node_box_global ;
    BoxType node_box_local_used ;
    std::vector< BoxType > node_box_parts ;

    element.create_node_boxes_from_vertex_boxes(
      vertex_box_global , vertex_box_parts ,
      node_box_global , node_box_parts );

    // Node communication maps:

    size_t node_count_interior = 0 ;
    size_t node_count_owned    = 0 ;
    size_t node_count_total    = 0 ;
    std::vector<size_t>                 node_used_id_map ;
    std::vector<size_t>                 node_part_counts ;
    std::vector< std::vector<size_t> >  node_send_map ;

    box_partition_maps( node_box_global ,
                        node_box_parts ,
                        element.box_bounds ,
                        proc_local ,
                        node_box_local_used ,
                        node_used_id_map ,
                        node_count_interior ,
                        node_count_owned ,
                        node_count_total ,
                        node_part_counts ,
                        node_send_map );

    size_t node_count_send = 0 ;
    for ( size_t i = 0 ; i < node_send_map.size() ; ++i ) {
      node_count_send += node_send_map[i].size();
    }

    size_t recv_msg_count = 0 ;
    size_t send_msg_count = 0 ;
    size_t send_count = 0 ;

    for ( size_t i = 1 ; i < proc_count ; ++i ) {
      if ( node_part_counts[i] ) ++recv_msg_count ;
      if ( node_send_map[i].size() ) {
        ++send_msg_count ;
        send_count += node_send_map[i].size();
      }
    }

    // Finite element mesh:

    FEMeshType mesh ;

    if ( node_count_total ) {
      mesh.node_coords = node_coords_type( "node_coords", node_count_total );
    }

    if ( elem_count_total ) {
      mesh.elem_node_ids =
        elem_node_ids_type( "elem_node_ids", elem_count_total );
    }

    mesh.parallel_data_map.assign( node_count_interior ,
                                   node_count_owned ,
                                   node_count_total ,
                                   recv_msg_count ,
                                   send_msg_count ,
                                   send_count );

    typename node_coords_type::HostMirror node_coords =
      Kokkos::create_mirror( mesh.node_coords );

    typename elem_node_ids_type::HostMirror elem_node_ids =
      Kokkos::create_mirror( mesh.elem_node_ids );

    //------------------------------------
    // set node coordinates to grid location for subsequent verification

    for ( size_t iz = node_box_local_used[2][0] ;
                 iz < node_box_local_used[2][1] ; ++iz ) {

    for ( size_t iy = node_box_local_used[1][0] ;
                 iy < node_box_local_used[1][1] ; ++iy ) {

    for ( size_t ix = node_box_local_used[0][0] ;
                 ix < node_box_local_used[0][1] ; ++ix ) {

      const size_t node_local_id =
        box_map_id( node_box_local_used , node_used_id_map , ix , iy , iz );

      node_coords( node_local_id , 0 ) = ix ;
      node_coords( node_local_id , 1 ) = iy ;
      node_coords( node_local_id , 2 ) = iz ;
    }}}

    //------------------------------------
    // Initialize element-node connectivity:

    if ( 1 < gang_count ) {
      layout_elements_partitioned( vertex_box_local_used ,
                                   vertex_box_local_owned ,
                                   node_box_local_used ,
                                   node_used_id_map ,
                                   element ,
                                   gang_count ,
                                   elem_node_ids );
    }
    else {
      layout_elements_interior_exterior( vertex_box_local_used ,
                                         vertex_box_local_owned ,
                                         node_box_local_used ,
                                         node_used_id_map ,
                                         element ,
                                         elem_count_interior ,
                                         elem_node_ids );
    }

    //------------------------------------
    // Populate node->element connectivity:

    std::vector<size_t> node_elem_work( node_count_total , (size_t) 0 );

    for ( size_t i = 0 ; i < elem_count_total ; ++i ) {
      for ( size_t n = 0 ; n < element_node_count  ; ++n ) {
        ++node_elem_work[ elem_node_ids(i,n) ];
      }
    }

    mesh.node_elem_ids =
      Kokkos::create_staticcrsgraph< node_elem_ids_type >( "node_elem_ids" , node_elem_work );

    typename node_elem_ids_type::HostMirror
      node_elem_ids = Kokkos::create_mirror( mesh.node_elem_ids );

    for ( size_t i = 0 ; i < node_count_total ; ++i ) {
      node_elem_work[i] = node_elem_ids.row_map[i];
    }

    // Looping in element order insures the list of elements
    // is sorted by element index.

    for ( size_t i = 0 ; i < elem_count_total ; ++i ) {
      for ( size_t n = 0 ; n < element_node_count ; ++n ) {
        const unsigned nid = elem_node_ids(i, n);
        const unsigned j = node_elem_work[nid] ; ++node_elem_work[nid] ;

        node_elem_ids.entries( j , 0 ) = i ;
        node_elem_ids.entries( j , 1 ) = n ;
      }
    }
    //------------------------------------
    // Verify setup with node coordinates matching grid indices.
    verify( node_coords , elem_node_ids , node_elem_ids );

    //------------------------------------
    // Scale node coordinates to problem extent with
    // nonlinear mapping.
    {
      const double problem_extent[3] =
        { static_cast<double>( vertex_box_global[0][1] - 1 ) ,
          static_cast<double>( vertex_box_global[1][1] - 1 ) ,
          static_cast<double>( vertex_box_global[2][1] - 1 ) };

      const double grid_extent[3] =
        { static_cast<double>( node_box_global[0][1] - 1 ) ,
          static_cast<double>( node_box_global[1][1] - 1 ) ,
          static_cast<double>( node_box_global[2][1] - 1 ) };

      for ( size_t i = 0 ; i < node_count_total ; ++i ) {
        const double x_unit = node_coords(i,0) / grid_extent[0] ;
        const double y_unit = node_coords(i,1) / grid_extent[1] ;
        const double z_unit = node_coords(i,2) / grid_extent[2] ;

        node_coords(i,0) = coordinate_scalar_type( problem_extent[0] * std::pow( x_unit , x_coord_curve ) );
        node_coords(i,1) = coordinate_scalar_type( problem_extent[1] * std::pow( y_unit , y_coord_curve ) );
        node_coords(i,2) = coordinate_scalar_type( problem_extent[2] * std::pow( z_unit , z_coord_curve ) );
      }
    }

    Kokkos::deep_copy( mesh.node_coords ,   node_coords );
    Kokkos::deep_copy( mesh.elem_node_ids , elem_node_ids );
    Kokkos::deep_copy( mesh.node_elem_ids.entries , node_elem_ids.entries );

    //------------------------------------
    // Communication lists:
    {
      recv_msg_count = 0 ;
      send_msg_count = 0 ;
      send_count = 0 ;

      for ( size_t i = 1 ; i < proc_count ; ++i ) {

        // Order sending starting with the local processor rank
        // to try to smooth out the amount of messages simultaneously
        // send to a particular processor.

        const int proc = ( proc_local + i ) % proc_count ;
        if ( node_part_counts[i] ) {
          mesh.parallel_data_map.host_recv(recv_msg_count,0) = proc ;
          mesh.parallel_data_map.host_recv(recv_msg_count,1) = node_part_counts[i] ;
          ++recv_msg_count ;
        }
        if ( node_send_map[i].size() ) {
          mesh.parallel_data_map.host_send(send_msg_count,0) = proc ;
          mesh.parallel_data_map.host_send(send_msg_count,1) = node_send_map[i].size() ;
          for ( size_t j = 0 ; j < node_send_map[i].size() ; ++j , ++send_count ) {
            mesh.parallel_data_map.host_send_item(send_count) = node_send_map[i][j] - node_count_interior ;
          }
          ++send_msg_count ;
        }
      }
    }

    return mesh ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_BOXMESHFIXTURE_HPP */


