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

#ifndef SPARSELINEARSYSTEMFILL_HPP
#define SPARSELINEARSYSTEMFILL_HPP

#include <vector>
#include <algorithm>
#include <limits>

#include <FEMesh.hpp>
#include <SparseLinearSystem.hpp>

//----------------------------------------------------------------------------

namespace HybridFEM {

template< class MatrixType , class MeshType ,
          class elem_matrices_type ,
          class elem_vectors_type > struct GatherFill ;


template< typename ScalarType ,
          class    DeviceType ,
          unsigned ElemNode ,
          typename CoordScalarType ,
          class elem_matrices_type ,
          class elem_vectors_type >
struct GatherFill< 
  Kokkos::CrsMatrix< ScalarType , DeviceType > ,
  FEMesh< CoordScalarType , ElemNode , DeviceType > ,
  elem_matrices_type , elem_vectors_type >
{
  typedef DeviceType     execution_space ;
  typedef typename execution_space::size_type  size_type ;

  static const size_type ElemNodeCount = ElemNode ;

  typedef Kokkos::CrsMatrix< ScalarType , execution_space >    matrix_type ;
  typedef typename matrix_type::coefficients_type   coefficients_type ;
  typedef Kokkos::View< ScalarType[] , execution_space >  vector_type ;
  typedef Kokkos::View< size_type[][ElemNodeCount][ElemNodeCount] , execution_space >       elem_graph_type ;

  typedef FEMesh< CoordScalarType , ElemNodeCount , execution_space > mesh_type ;
  typedef typename mesh_type::node_elem_ids_type node_elem_ids_type ;

private:

  node_elem_ids_type  node_elem_ids ;
  elem_graph_type     elem_graph ;
  elem_matrices_type  elem_matrices ;
  elem_vectors_type   elem_vectors ;
  coefficients_type   system_coeff ;
  vector_type         system_rhs ;

public:

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type irow ) const
  {
    const size_type node_elem_begin = node_elem_ids.row_map[irow];
    const size_type node_elem_end   = node_elem_ids.row_map[irow+1];

    //  for each element that a node belongs to 

    for ( size_type i = node_elem_begin ; i < node_elem_end ; i++ ) {

      const size_type elem_id   = node_elem_ids.entries( i, 0);
      const size_type row_index = node_elem_ids.entries( i, 1);

      system_rhs(irow) += elem_vectors(elem_id, row_index);

      //  for each node in a particular related element  
      //  gather the contents of the element stiffness
      //  matrix that belong in irow

      for ( size_type j = 0 ; j < ElemNodeCount ; ++j ){
        const size_type A_index = elem_graph( elem_id , row_index , j );

        system_coeff( A_index ) += elem_matrices( elem_id, row_index, j );
      }
    }
  }


  static void apply( const matrix_type & matrix ,
                     const vector_type & rhs ,
                     const mesh_type   & mesh ,
                     const elem_graph_type    & elem_graph ,
                     const elem_matrices_type & elem_matrices ,
                     const elem_vectors_type  & elem_vectors )
  {
    const size_t row_count = matrix.graph.row_map.dimension_0() - 1 ;
    GatherFill op ;
    op.node_elem_ids = mesh.node_elem_ids ;
    op.elem_graph    = elem_graph ;
    op.elem_matrices = elem_matrices ;
    op.elem_vectors  = elem_vectors ;
    op.system_coeff  = matrix.coefficients ;
    op.system_rhs    = rhs ;

    parallel_for( row_count , op );
  }
};

} /* namespace HybridFEM */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace HybridFEM {

template< class GraphType , class MeshType >
struct GraphFactory {
  typedef GraphType                         graph_type ;
  typedef MeshType                          mesh_type ;
  typedef typename graph_type::execution_space  execution_space ;
  typedef typename execution_space::size_type   size_type  ;

  static const unsigned ElemNodeCount = mesh_type::element_node_count ;

  typedef Kokkos::View< size_type[][ElemNodeCount][ElemNodeCount] , execution_space >         element_map_type ;

  static
  void
  create( const mesh_type & mesh ,
          graph_type & graph ,
          element_map_type & elem_map )
  {
    typename mesh_type::node_elem_ids_type::HostMirror
      node_elem_ids = create_mirror( mesh.node_elem_ids );

    typename mesh_type::elem_node_ids_type::HostMirror
      elem_node_ids = create_mirror( mesh.elem_node_ids );

    typedef typename element_map_type::HostMirror element_map_host_type ;

    deep_copy( elem_node_ids , mesh.elem_node_ids );
    deep_copy( node_elem_ids.entries , mesh.node_elem_ids.entries );

    const size_t owned_node = mesh.parallel_data_map.count_owned ;
    const size_t total_elem = mesh.elem_node_ids.dimension_0();

    if ( total_elem ) {
      elem_map = element_map_type( std::string("element_map"), total_elem );
    }

    element_map_host_type elem_map_host = create_mirror( elem_map );

    //------------------------------------
    //  Node->node mapping for the CrsMatrix graph

    std::vector< std::vector< unsigned > > node_node_ids( owned_node );
    std::vector< unsigned > node_node_begin( owned_node );

    size_t offset = 0 ;
    for ( size_t i = 0 ; i < owned_node ; ++i ) {
      const size_t j_end = node_elem_ids.row_map[i+1];
            size_t j     = node_elem_ids.row_map[i];

      node_node_begin[i] = offset ;

      std::vector< unsigned > & work = node_node_ids[i] ;

      for ( ; j < j_end ; ++j ) {
        const size_t elem_id = node_elem_ids.entries(j,0);
        for ( size_t k = 0 ; k < ElemNodeCount ; ++k ) {
          work.push_back( elem_node_ids( elem_id , k ) );
        }
      }

      std::sort( work.begin() , work.end() );

      work.erase( std::unique( work.begin() , work.end() ) , work.end() );

      offset += work.size();
    }

    graph = Kokkos::create_staticcrsgraph< graph_type >( "node_node_ids" , node_node_ids );

    //------------------------------------
    // ( element , node_row , node_column ) -> matrix_crs_column

    for ( size_t elem_id = 0 ; elem_id < total_elem ; ++elem_id ) {
      for ( size_t i = 0 ; i < ElemNodeCount ; ++i ) {

        const size_t node_row = elem_node_ids( elem_id , i );
        const size_t node_row_begin = node_node_begin[ node_row ];
        const std::vector< unsigned > & column = node_node_ids[ node_row ] ;

        if ( owned_node <= node_row ) {
          for ( unsigned j = 0 ; j < ElemNodeCount ; ++j ) {
            elem_map_host( elem_id , i , j ) = std::numeric_limits<size_type>::max();
          }
        }
        else {

          for ( unsigned j = 0 ; j < ElemNodeCount ; ++j ) {
            const size_type node_col = elem_node_ids( elem_id , j );

            int col_search = 0 ;

            for ( int len = column.size() ; 0 < len ; ) {

              const int half = len >> 1;
              const int middle = col_search + half ;

              if ( column[middle] < node_col ){
                col_search = middle + 1 ;
                len -= half + 1 ;
              }
              else {
                len = half ;
              }
            }
if ( node_col != column[col_search] ) {
  throw std::runtime_error(std::string("Failed"));
}
            elem_map_host( elem_id , i , j ) = col_search + node_row_begin ;
          }
        }
      }
    }

    deep_copy( elem_map , elem_map_host );
  }
};

} // namespace HybridFEM


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef SPARSELINEARSYSTEMFILL_HPP */

