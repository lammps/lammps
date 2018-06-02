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

#ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP
#define KOKKOS_EXAMPLE_FENLFUNCTORS_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <Kokkos_Pair.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <impl/Kokkos_Timer.hpp>

#include <BoxElemFixture.hpp>
#include <HexElement.hpp>
#include <CGSolve.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template< class ElemNodeIdView , class CrsGraphType , unsigned ElemNode >
class NodeNodeGraph {
public:

  typedef typename ElemNodeIdView::execution_space  execution_space ;
  typedef pair<unsigned,unsigned> key_type ;

  typedef Kokkos::UnorderedMap< key_type, void , execution_space >  SetType ;
  typedef typename CrsGraphType::row_map_type::non_const_type       RowMapType ;
  typedef Kokkos::View< unsigned ,  execution_space >               UnsignedValue ;

  // Static dimensions of 0 generate compiler warnings or errors.
  typedef Kokkos::View< unsigned*[ElemNode][ElemNode] , execution_space >
    ElemGraphType ;

  struct TagFillNodeSet {};
  struct TagScanNodeCount {};
  struct TagFillGraphEntries {};
  struct TagSortGraphEntries {};
  struct TagFillElementGraph {};

private:

  enum PhaseType { FILL_NODE_SET ,
                   SCAN_NODE_COUNT ,
                   FILL_GRAPH_ENTRIES ,
                   SORT_GRAPH_ENTRIES ,
                   FILL_ELEMENT_GRAPH };

  const unsigned        node_count ;
  const ElemNodeIdView  elem_node_id ;
  UnsignedValue         row_total ;
  RowMapType            row_count ;
  RowMapType            row_map ;
  SetType               node_node_set ;
  PhaseType             phase ;

public:

  CrsGraphType          graph ;
  ElemGraphType         elem_graph ;

  struct Times
  {
    double ratio;
    double fill_node_set;
    double scan_node_count;
    double fill_graph_entries;
    double sort_graph_entries;
    double fill_element_graph;
  };

  NodeNodeGraph( const ElemNodeIdView & arg_elem_node_id ,
                 const unsigned         arg_node_count,
                 Times & results
               )
    : node_count(arg_node_count)
    , elem_node_id( arg_elem_node_id )
    , row_total( "row_total" )
    , row_count(Kokkos::ViewAllocateWithoutInitializing("row_count") , node_count ) // will deep_copy to 0 inside loop
    , row_map( "graph_row_map" , node_count + 1 )
    , node_node_set()
    , phase( FILL_NODE_SET )
    , graph()
    , elem_graph()
   {
      //--------------------------------
      // Guess at capacity required for the map:

      Kokkos::Timer wall_clock ;

      wall_clock.reset();
      phase = FILL_NODE_SET ;

      // upper bound on the capacity
      size_t set_capacity = (28ull * node_count) / 2;
      unsigned failed_insert_count = 0 ;

      do {
        // Zero the row count to restart the fill
        Kokkos::deep_copy( row_count , 0u );

        node_node_set = SetType( ( set_capacity += failed_insert_count ) );

        // May be larger that requested:
        set_capacity = node_node_set.capacity();

        Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space,TagFillNodeSet>(0,elem_node_id.extent(0))
                               , *this
                               , failed_insert_count );

      } while ( failed_insert_count );

      execution_space::fence();
      results.ratio = (double)node_node_set.size() / (double)node_node_set.capacity();
      results.fill_node_set = wall_clock.seconds();
      //--------------------------------

      wall_clock.reset();
      phase = SCAN_NODE_COUNT ;

      // Exclusive scan of row_count into row_map
      // including the final total in the 'node_count + 1' position.
      // Zero the 'row_count' values.
      Kokkos::parallel_scan( node_count , *this );

      // Zero the row count for the fill:
      Kokkos::deep_copy( row_count , 0u );

      unsigned graph_entry_count = 0 ;

      Kokkos::deep_copy( graph_entry_count , row_total );

      // Assign graph's row_map and allocate graph's entries
      graph.row_map = row_map ;
      graph.entries = typename CrsGraphType::entries_type( "graph_entries" , graph_entry_count );

      //--------------------------------
      // Fill graph's entries from the (node,node) set.

      execution_space::fence();
      results.scan_node_count = wall_clock.seconds();

      wall_clock.reset();
      phase = FILL_GRAPH_ENTRIES ;
      Kokkos::parallel_for( node_node_set.capacity() , *this );

      execution_space::fence();
      results.fill_graph_entries = wall_clock.seconds();

      //--------------------------------
      // Done with the temporary sets and arrays
      wall_clock.reset();
      phase = SORT_GRAPH_ENTRIES ;

      row_total = UnsignedValue();
      row_count = RowMapType();
      row_map   = RowMapType();
      node_node_set.clear();

      //--------------------------------

      Kokkos::parallel_for( node_count , *this );

      execution_space::fence();
      results.sort_graph_entries = wall_clock.seconds();

      //--------------------------------
      // Element-to-graph mapping:
      wall_clock.reset();
      phase = FILL_ELEMENT_GRAPH ;
      elem_graph = ElemGraphType("elem_graph", elem_node_id.extent(0) );
      Kokkos::parallel_for( elem_node_id.extent(0) , *this );

      execution_space::fence();
      results.fill_element_graph = wall_clock.seconds();
    }

  //------------------------------------
  // parallel_for: create map and count row length

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagFillNodeSet & , unsigned ielem , unsigned & count ) const
  {
    // Loop over element's (row_local_node,col_local_node) pairs:
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = row_local_node ; col_local_node < elem_node_id.extent(1) ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        // If either node is locally owned then insert the pair into the unordered map:

        if ( row_node < row_count.extent(0) || col_node < row_count.extent(0) ) {

          const key_type key = (row_node < col_node) ? make_pair( row_node, col_node ) : make_pair( col_node, row_node ) ;

          const typename SetType::insert_result result = node_node_set.insert( key );

          // A successfull insert: the first time this pair was added
          if ( result.success() ) {

            // If row node is owned then increment count
            if ( row_node < row_count.extent(0) ) { atomic_fetch_add( & row_count( row_node ) , 1 ); }

            // If column node is owned and not equal to row node then increment count
            if ( col_node < row_count.extent(0) && col_node != row_node ) { atomic_fetch_add( & row_count( col_node ) , 1 ); }
          }
          else if ( result.failed() ) {
            ++count ;
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_graph_entries( const unsigned iset ) const
  {
    if ( node_node_set.valid_at(iset) ) {
      // Add each entry to the graph entries.

      const key_type key = node_node_set.key_at(iset) ;
      const unsigned row_node = key.first ;
      const unsigned col_node = key.second ;

      if ( row_node < row_count.extent(0) ) {
        const unsigned offset = graph.row_map( row_node ) + atomic_fetch_add( & row_count( row_node ) , 1 );
        graph.entries( offset ) = col_node ;
      }

      if ( col_node < row_count.extent(0) && col_node != row_node ) {
        const unsigned offset = graph.row_map( col_node ) + atomic_fetch_add( & row_count( col_node ) , 1 );
        graph.entries( offset ) = row_node ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void sort_graph_entries( const unsigned irow ) const
  {
    const unsigned row_beg = graph.row_map( irow );
    const unsigned row_end = graph.row_map( irow + 1 );
    for ( unsigned i = row_beg + 1 ; i < row_end ; ++i ) {
      const unsigned col = graph.entries(i);
      unsigned j = i ;
      for ( ; row_beg < j && col < graph.entries(j-1) ; --j ) {
        graph.entries(j) = graph.entries(j-1);
      }
      graph.entries(j) = col ;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_elem_graph_map( const unsigned ielem ) const
  {
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = 0 ; col_local_node < elem_node_id.extent(1) ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        unsigned entry = ~0u ;

        if ( row_node + 1 < graph.row_map.extent(0) ) {

          const unsigned entry_end = graph.row_map( row_node + 1 );

          entry = graph.row_map( row_node );

          for ( ; entry < entry_end && graph.entries(entry) != col_node ; ++entry );

          if ( entry == entry_end ) entry = ~0u ;
        }

        elem_graph( ielem , row_local_node , col_local_node ) = entry ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned iwork ) const
  {
/*
    if ( phase == FILL_NODE_SET ) {
      operator()( TagFillNodeSet() , iwork );
    }
    else */
    if ( phase == FILL_GRAPH_ENTRIES ) {
      fill_graph_entries( iwork );
    }
    else if ( phase == SORT_GRAPH_ENTRIES ) {
      sort_graph_entries( iwork );
    }
    else if ( phase == FILL_ELEMENT_GRAPH ) {
      fill_elem_graph_map( iwork );
    }
  }

  //------------------------------------
  // parallel_scan: row offsets

  typedef unsigned value_type ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned irow , unsigned & update , const bool final ) const
  {
    // exclusive scan
    if ( final ) { row_map( irow ) = update ; }

    update += row_count( irow );

    if ( final ) {
      if ( irow + 1 == row_count.extent(0) ) {
        row_map( irow + 1 ) = update ;
        row_total()         = update ;
      }
    }
  }

  // For the reduce phase:
  KOKKOS_INLINE_FUNCTION
  void init( const TagFillNodeSet & , unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( const TagFillNodeSet &
           , volatile       unsigned & update
           , volatile const unsigned & input ) const { update += input ; }

  // For the scan phase::
  KOKKOS_INLINE_FUNCTION
  void init( unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile       unsigned & update
           , volatile const unsigned & input ) const { update += input ; }

  //------------------------------------
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template< class ElemCompType >
class NodeElemGatherFill {
public:

  typedef typename ElemCompType::execution_space         execution_space ;
  typedef typename ElemCompType::vector_type         vector_type ;
  typedef typename ElemCompType::sparse_matrix_type  sparse_matrix_type ;
  typedef typename ElemCompType::elem_node_type      elem_node_type ;
  typedef typename ElemCompType::elem_vectors_type   elem_vectors_type ;
  typedef typename ElemCompType::elem_matrices_type  elem_matrices_type ;
  typedef typename ElemCompType::elem_graph_type     elem_graph_type ;

  static const unsigned ElemNodeCount = ElemCompType::ElemNodeCount ;

  //------------------------------------

private:

  typedef Kokkos::StaticCrsGraph< unsigned[2] , execution_space >  CrsGraphType ;
  typedef typename CrsGraphType::row_map_type::non_const_type  RowMapType ;
  typedef Kokkos::View< unsigned ,  execution_space >              UnsignedValue ;

  enum PhaseType { FILL_NODE_COUNT ,
                   SCAN_NODE_COUNT ,
                   FILL_GRAPH_ENTRIES ,
                   SORT_GRAPH_ENTRIES ,
                   GATHER_FILL };

  const elem_node_type  elem_node_id ;
  const elem_graph_type elem_graph ;
  UnsignedValue         row_total ;
  RowMapType            row_count ;
  RowMapType            row_map ;
  CrsGraphType          graph ;
  vector_type           residual ;
  sparse_matrix_type    jacobian ;
  elem_vectors_type     elem_residual ;
  elem_matrices_type    elem_jacobian ;
  PhaseType             phase ;

public:

  NodeElemGatherFill()
    : elem_node_id()
    , elem_graph()
    , row_total()
    , row_count()
    , row_map()
    , graph()
    , residual()
    , jacobian()
    , elem_residual()
    , elem_jacobian()
    , phase( FILL_NODE_COUNT )
    {}

  NodeElemGatherFill( const NodeElemGatherFill & rhs )
    : elem_node_id(  rhs.elem_node_id )
    , elem_graph(    rhs.elem_graph )
    , row_total(     rhs.row_total )
    , row_count(     rhs.row_count )
    , row_map(       rhs.row_map )
    , graph(         rhs.graph )
    , residual(      rhs.residual )
    , jacobian(      rhs.jacobian )
    , elem_residual( rhs.elem_residual )
    , elem_jacobian( rhs.elem_jacobian )
    , phase(         rhs.phase )
    {}

  NodeElemGatherFill( const elem_node_type     & arg_elem_node_id ,
                      const elem_graph_type    & arg_elem_graph ,
                      const vector_type        & arg_residual ,
                      const sparse_matrix_type & arg_jacobian ,
                      const elem_vectors_type  & arg_elem_residual ,
                      const elem_matrices_type & arg_elem_jacobian )
    : elem_node_id( arg_elem_node_id )
    , elem_graph( arg_elem_graph )
    , row_total( "row_total" )
    , row_count( "row_count" , arg_residual.extent(0) )
    , row_map( "graph_row_map" , arg_residual.extent(0) + 1 )
    , graph()
    , residual( arg_residual )
    , jacobian( arg_jacobian )
    , elem_residual( arg_elem_residual )
    , elem_jacobian( arg_elem_jacobian )
    , phase( FILL_NODE_COUNT )
    {
      //--------------------------------
      // Count node->element relations

      phase = FILL_NODE_COUNT ;

      Kokkos::parallel_for( elem_node_id.extent(0) , *this );

      //--------------------------------

      phase = SCAN_NODE_COUNT ;

      // Exclusive scan of row_count into row_map
      // including the final total in the 'node_count + 1' position.
      // Zero the 'row_count' values.
      Kokkos::parallel_scan( residual.extent(0) , *this );

      // Zero the row count for the fill:
      Kokkos::deep_copy( row_count , typename RowMapType::value_type(0) );

      unsigned graph_entry_count = 0 ;

      Kokkos::deep_copy( graph_entry_count , row_total );

      // Assign graph's row_map and allocate graph's entries
      graph.row_map = row_map ;

      typedef typename CrsGraphType::entries_type graph_entries_type ;

      graph.entries = graph_entries_type( "graph_entries" , graph_entry_count );

      //--------------------------------
      // Fill graph's entries from the (node,node) set.

      phase = FILL_GRAPH_ENTRIES ;

      Kokkos::deep_copy( row_count , 0u );
      Kokkos::parallel_for( elem_node_id.extent(0) , *this );

      execution_space::fence();

      //--------------------------------
      // Done with the temporary sets and arrays

      row_total = UnsignedValue();
      row_count = RowMapType();
      row_map   = RowMapType();

      //--------------------------------

      phase = SORT_GRAPH_ENTRIES ;
      Kokkos::parallel_for( residual.extent(0) , *this );

      execution_space::fence();

      phase = GATHER_FILL ;
    }

  void apply() const
  {
    Kokkos::parallel_for( residual.extent(0) , *this );
  }

  //------------------------------------
  //------------------------------------
  // parallel_for: Count node->element pairs

  KOKKOS_INLINE_FUNCTION
  void fill_node_count( const unsigned ielem ) const
  {
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      if ( row_node < row_count.extent(0) ) {
        atomic_fetch_add( & row_count( row_node ) , 1 );
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_graph_entries( const unsigned ielem ) const
  {
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      if ( row_node < row_count.extent(0) ) {

        const unsigned offset = graph.row_map( row_node ) + atomic_fetch_add( & row_count( row_node ) , 1 );

        graph.entries( offset , 0 ) = ielem ;
        graph.entries( offset , 1 ) = row_local_node ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void sort_graph_entries( const unsigned irow ) const
  {
    const unsigned row_beg = graph.row_map( irow );
    const unsigned row_end = graph.row_map( irow + 1 );
    for ( unsigned i = row_beg + 1 ; i < row_end ; ++i ) {
      const unsigned elem  = graph.entries(i,0);
      const unsigned local = graph.entries(i,1);
      unsigned j = i ;
      for ( ; row_beg < j && elem < graph.entries(j-1,0) ; --j ) {
        graph.entries(j,0) = graph.entries(j-1,0);
        graph.entries(j,1) = graph.entries(j-1,1);
      }
      graph.entries(j,0) = elem ;
      graph.entries(j,1) = local ;
    }
  }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  void gather_fill( const unsigned irow ) const
  {
    const unsigned node_elem_begin = graph.row_map(irow);
    const unsigned node_elem_end   = graph.row_map(irow+1);

    //  for each element that a node belongs to

    for ( unsigned i = node_elem_begin ; i < node_elem_end ; i++ ) {

      const unsigned elem_id   = graph.entries( i, 0);
      const unsigned row_index = graph.entries( i, 1);

      residual(irow) += elem_residual(elem_id, row_index);

      //  for each node in a particular related element
      //  gather the contents of the element stiffness
      //  matrix that belong in irow

      for ( unsigned j = 0 ; j < ElemNodeCount ; ++j ) {
        const unsigned A_index = elem_graph( elem_id , row_index , j );

        jacobian.coeff( A_index ) += elem_jacobian( elem_id, row_index, j );
      }
    }
  }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned iwork ) const
  {
    if ( phase == FILL_NODE_COUNT ) {
      fill_node_count( iwork );
    }
    else if ( phase == FILL_GRAPH_ENTRIES ) {
      fill_graph_entries( iwork );
    }
    else if ( phase == SORT_GRAPH_ENTRIES ) {
      sort_graph_entries( iwork );
    }
    else if ( phase == GATHER_FILL ) {
      gather_fill( iwork );
    }
  }

  //------------------------------------
  // parallel_scan: row offsets

  typedef unsigned value_type ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned irow , unsigned & update , const bool final ) const
  {
    // exclusive scan
    if ( final ) { row_map( irow ) = update ; }

    update += row_count( irow );

    if ( final ) {
      if ( irow + 1 == row_count.extent(0) ) {
        row_map( irow + 1 ) = update ;
        row_total()         = update ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile unsigned & update , const volatile unsigned & input ) const { update += input ; }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template< class FiniteElementMeshType , class SparseMatrixType >
class ElementComputation ;


template< class ExecSpace , BoxElemPart::ElemOrder Order , class CoordinateMap , typename ScalarType >
class ElementComputation<
  Kokkos::Example::BoxElemFixture< ExecSpace , Order , CoordinateMap > ,
  Kokkos::Example::CrsMatrix< ScalarType , ExecSpace > >
{
public:

  typedef Kokkos::Example::BoxElemFixture< ExecSpace, Order, CoordinateMap >  mesh_type ;
  typedef Kokkos::Example::HexElement_Data< mesh_type::ElemNode >             element_data_type ;

  typedef Kokkos::Example::CrsMatrix< ScalarType , ExecSpace >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType       sparse_graph_type ;

  typedef ExecSpace   execution_space ;
  typedef ScalarType  scalar_type ;

  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;
  static const unsigned ElemNodeCount    = element_data_type::element_node_count ;
  static const unsigned FunctionCount    = element_data_type::function_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;

  //------------------------------------

  typedef typename mesh_type::node_coord_type                                      node_coord_type ;
  typedef typename mesh_type::elem_node_type                                       elem_node_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount][FunctionCount] , execution_space > elem_matrices_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount] ,                execution_space > elem_vectors_type ;
  typedef Kokkos::View< scalar_type* ,                               execution_space > vector_type ;

  typedef typename NodeNodeGraph< elem_node_type , sparse_graph_type , ElemNodeCount >::ElemGraphType elem_graph_type ;

  //------------------------------------


  //------------------------------------
  // Computational data:

  const element_data_type   elem_data ;
  const elem_node_type      elem_node_ids ;
  const node_coord_type     node_coords ;
  const elem_graph_type     elem_graph ;
  const elem_matrices_type  elem_jacobians ;
  const elem_vectors_type   elem_residuals ;
  const vector_type         solution ;
  const vector_type         residual ;
  const sparse_matrix_type  jacobian ;
  const scalar_type         coeff_K ;

  ElementComputation( const ElementComputation & rhs )
    : elem_data()
    , elem_node_ids( rhs.elem_node_ids )
    , node_coords(   rhs.node_coords )
    , elem_graph(    rhs.elem_graph )
    , elem_jacobians( rhs.elem_jacobians )
    , elem_residuals( rhs.elem_residuals )
    , solution( rhs.solution )
    , residual( rhs.residual )
    , jacobian( rhs.jacobian )
    , coeff_K( rhs.coeff_K )
    {}

  // If the element->sparse_matrix graph is provided then perform atomic updates
  // Otherwise fill per-element contributions for subequent gather-add into a residual and jacobian.
  ElementComputation( const mesh_type          & arg_mesh ,
	              const scalar_type          arg_coeff_K ,
                      const vector_type        & arg_solution ,
                      const elem_graph_type    & arg_elem_graph ,
                      const sparse_matrix_type & arg_jacobian ,
                      const vector_type        & arg_residual )
    : elem_data()
    , elem_node_ids( arg_mesh.elem_node() )
    , node_coords(   arg_mesh.node_coord() )
    , elem_graph(    arg_elem_graph )
    , elem_jacobians()
    , elem_residuals()
    , solution( arg_solution )
    , residual( arg_residual )
    , jacobian( arg_jacobian )
    , coeff_K( arg_coeff_K )
    {}

  ElementComputation( const mesh_type    & arg_mesh ,
	              const scalar_type    arg_coeff_K ,
                      const vector_type  & arg_solution )
    : elem_data()
    , elem_node_ids( arg_mesh.elem_node() )
    , node_coords(   arg_mesh.node_coord() )
    , elem_graph()
    , elem_jacobians( "elem_jacobians" , arg_mesh.elem_count() )
    , elem_residuals( "elem_residuals" , arg_mesh.elem_count() )
    , solution( arg_solution )
    , residual()
    , jacobian()
    , coeff_K( arg_coeff_K )
    {}

  //------------------------------------

  void apply() const
  {
    parallel_for( elem_node_ids.extent(0) , *this );
  }

  //------------------------------------

  static const unsigned FLOPS_transform_gradients =
     /* Jacobian */           FunctionCount * TensorDim * 2 +
     /* Inverse jacobian */   TensorDim * 6 + 6 +
     /* Gradient transform */ FunctionCount * 15 ;

  KOKKOS_INLINE_FUNCTION
  float transform_gradients(
    const float grad[][ FunctionCount ] , // Gradient of bases master element
    const double x[] ,
    const double y[] ,
    const double z[] ,
    float dpsidx[] ,
    float dpsidy[] ,
    float dpsidz[] ) const
  {
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Jacobian accumulation:

    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const double x1 = x[i] ;
      const double x2 = y[i] ;
      const double x3 = z[i] ;

      const float g1 = grad[0][i] ;
      const float g2 = grad[1][i] ;
      const float g3 = grad[2][i] ;

      J[j11] += g1 * x1 ;
      J[j12] += g1 * x2 ;
      J[j13] += g1 * x3 ;

      J[j21] += g2 * x1 ;
      J[j22] += g2 * x2 ;
      J[j23] += g2 * x3 ;

      J[j31] += g3 * x1 ;
      J[j32] += g3 * x2 ;
      J[j33] += g3 * x3 ;
    }

    // Inverse jacobian:

    float invJ[ TensorDim ] = {
      static_cast<float>( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      static_cast<float>( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      static_cast<float>( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      static_cast<float>( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      static_cast<float>( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      static_cast<float>( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      static_cast<float>( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      static_cast<float>( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      static_cast<float>( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const float detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    const float detJinv = 1.0 / detJ ;

    for ( unsigned i = 0 ; i < TensorDim ; ++i ) { invJ[i] *= detJinv ; }

    // Transform gradients:

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const float g0 = grad[0][i];
      const float g1 = grad[1][i];
      const float g2 = grad[2][i];

      dpsidx[i] = g0 * invJ[j11] + g1 * invJ[j12] + g2 * invJ[j13];
      dpsidy[i] = g0 * invJ[j21] + g1 * invJ[j22] + g2 * invJ[j23];
      dpsidz[i] = g0 * invJ[j31] + g1 * invJ[j32] + g2 * invJ[j33];
    }

    return detJ ;
  }

  KOKKOS_INLINE_FUNCTION
  void contributeResidualJacobian(
    const float coeff_k ,
    const double dof_values[] ,
    const float dpsidx[] ,
    const float dpsidy[] ,
    const float dpsidz[] ,
    const float detJ ,
    const float integ_weight ,
    const float bases_vals[] ,
    double elem_res[] ,
    double elem_mat[][ FunctionCount ] ) const
  {
    double value_at_pt = 0 ;
    double gradx_at_pt = 0 ;
    double grady_at_pt = 0 ;
    double gradz_at_pt = 0 ;

    for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
      value_at_pt += dof_values[m] * bases_vals[m] ;
      gradx_at_pt += dof_values[m] * dpsidx[m] ;
      grady_at_pt += dof_values[m] * dpsidy[m] ;
      gradz_at_pt += dof_values[m] * dpsidz[m] ;
    }

    const scalar_type k_detJ_weight = coeff_k        * detJ * integ_weight ;
    const double res_val = value_at_pt * value_at_pt * detJ * integ_weight ;
    const double mat_val = 2.0 * value_at_pt         * detJ * integ_weight ;

    // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$
    // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$

    for ( unsigned m = 0; m < FunctionCount; ++m) {
      double * const mat = elem_mat[m] ;
      const float bases_val_m = bases_vals[m];
      const float dpsidx_m    = dpsidx[m] ;
      const float dpsidy_m    = dpsidy[m] ;
      const float dpsidz_m    = dpsidz[m] ;

      elem_res[m] += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                       dpsidy_m * grady_at_pt +
                                       dpsidz_m * gradz_at_pt ) +
                     res_val * bases_val_m ;

      for( unsigned n = 0; n < FunctionCount; n++) {

        mat[n] += k_detJ_weight * ( dpsidx_m * dpsidx[n] +
                                    dpsidy_m * dpsidy[n] +
                                    dpsidz_m * dpsidz[n] ) +
                  mat_val * bases_val_m * bases_vals[n];
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    // Gather nodal coordinates and solution vector:

    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    double val[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = node_coords( ni , 0 );
      y[i] = node_coords( ni , 1 );
      z[i] = node_coords( ni , 2 );

      val[i] = solution( ni );
    }


    double elem_vec[ FunctionCount ] ;
    double elem_mat[ FunctionCount ][ FunctionCount ] ;

    for( unsigned i = 0; i < FunctionCount ; i++ ) {
      elem_vec[i] = 0 ;
      for( unsigned j = 0; j < FunctionCount ; j++){
        elem_mat[i][j] = 0 ;
      }
    }


    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {
      float dpsidx[ FunctionCount ] ;
      float dpsidy[ FunctionCount ] ;
      float dpsidz[ FunctionCount ] ;

      const float detJ =
        transform_gradients( elem_data.gradients[i] , x , y , z ,
                             dpsidx , dpsidy , dpsidz );

      contributeResidualJacobian( coeff_K ,
                                  val , dpsidx , dpsidy , dpsidz ,
                                  detJ ,
                                  elem_data.weights[i] ,
                                  elem_data.values[i] ,
                                  elem_vec , elem_mat );
    }

#if 0

if ( 1 == ielem ) {
  printf("ElemResidual { %f %f %f %f %f %f %f %f }\n",
         elem_vec[0], elem_vec[1], elem_vec[2], elem_vec[3],
         elem_vec[4], elem_vec[5], elem_vec[6], elem_vec[7]);

  printf("ElemJacobian {\n");

  for ( unsigned j = 0 ; j < FunctionCount ; ++j ) {
  printf("  { %f %f %f %f %f %f %f %f }\n",
         elem_mat[j][0], elem_mat[j][1], elem_mat[j][2], elem_mat[j][3],
         elem_mat[j][4], elem_mat[j][5], elem_mat[j][6], elem_mat[j][7]);
  }
  printf("}\n");
}

#endif

    if ( ! residual.extent(0) ) {
      for( unsigned i = 0; i < FunctionCount ; i++){
        elem_residuals(ielem, i) = elem_vec[i] ;
        for( unsigned j = 0; j < FunctionCount ; j++){
          elem_jacobians(ielem, i, j) = elem_mat[i][j] ;
        }
      }
    }
    else {
      for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
        const unsigned row = node_index[i] ;
        if ( row < residual.extent(0) ) {
          atomic_fetch_add( & residual( row ) , elem_vec[i] );

          for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
            const unsigned entry = elem_graph( ielem , i , j );
            if ( entry != ~0u ) {
              atomic_fetch_add( & jacobian.coeff( entry ) , elem_mat[i][j] );
            }
          }
        }
      }
    }
  }
}; /* ElementComputation */

//----------------------------------------------------------------------------

template< class FixtureType , class SparseMatrixType >
class DirichletComputation ;

template< class ExecSpace , BoxElemPart::ElemOrder Order , class CoordinateMap , typename ScalarType >
class DirichletComputation<
  Kokkos::Example::BoxElemFixture< ExecSpace , Order , CoordinateMap > ,
  Kokkos::Example::CrsMatrix< ScalarType , ExecSpace > >
{
public:

  typedef Kokkos::Example::BoxElemFixture< ExecSpace, Order, CoordinateMap >  mesh_type ;
  typedef typename mesh_type::node_coord_type                                 node_coord_type ;
  typedef typename node_coord_type::value_type                                scalar_coord_type ;

  typedef Kokkos::Example::CrsMatrix< ScalarType , ExecSpace >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType       sparse_graph_type ;

  typedef ExecSpace   execution_space ;
  typedef ScalarType  scalar_type ;

  //------------------------------------

  typedef Kokkos::View< scalar_type* , execution_space > vector_type ;

  //------------------------------------
  // Computational data:

  const node_coord_type     node_coords ;
  const vector_type         solution ;
  const sparse_matrix_type  jacobian ;
  const vector_type         residual ;
  const scalar_type         bc_lower_value ;
  const scalar_type         bc_upper_value ;
  const scalar_coord_type   bc_lower_limit ;
  const scalar_coord_type   bc_upper_limit ;
  const unsigned            bc_plane ;
  const unsigned            node_count ;
        bool                init ;


  DirichletComputation( const mesh_type          & arg_mesh ,
                        const vector_type        & arg_solution ,
                        const sparse_matrix_type & arg_jacobian ,
                        const vector_type        & arg_residual ,
                        const unsigned             arg_bc_plane ,
                        const scalar_type          arg_bc_lower_value ,
                        const scalar_type          arg_bc_upper_value )
    : node_coords( arg_mesh.node_coord() )
    , solution(    arg_solution )
    , jacobian(    arg_jacobian )
    , residual(    arg_residual )
    , bc_lower_value( arg_bc_lower_value )
    , bc_upper_value( arg_bc_upper_value )
    , bc_lower_limit( std::numeric_limits<scalar_coord_type>::epsilon() )
    , bc_upper_limit( scalar_coord_type(1) - std::numeric_limits<scalar_coord_type>::epsilon() )
    , bc_plane(       arg_bc_plane )
    , node_count( arg_mesh.node_count_owned() )
    , init( false )
    {
      parallel_for( node_count , *this );
      init = true ;
    }

  void apply() const
  {
    parallel_for( node_count , *this );
  }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned inode ) const
  {
    //  Apply dirichlet boundary condition on the Solution and Residual vectors.
    //  To maintain the symmetry of the original global stiffness matrix,
    //  zero out the columns that correspond to boundary conditions, and
    //  update the residual vector accordingly

    const unsigned iBeg = jacobian.graph.row_map[inode];
    const unsigned iEnd = jacobian.graph.row_map[inode+1];

    const scalar_coord_type c = node_coords(inode,bc_plane);
    const bool bc_lower = c <= bc_lower_limit ;
    const bool bc_upper = bc_upper_limit <= c ;

    if ( ! init ) {
      solution(inode) = bc_lower ? bc_lower_value : (
                        bc_upper ? bc_upper_value : 0 );
    }
    else {
      if ( bc_lower || bc_upper ) {

        residual(inode) = 0 ;

        //  zero each value on the row, and leave a one
        //  on the diagonal

        for( unsigned i = iBeg ; i < iEnd ; ++i ) {
          jacobian.coeff(i) = int(inode) == int(jacobian.graph.entries(i)) ? 1 : 0 ;
        }
      }
      else {

        //  Find any columns that are boundary conditions.
        //  Clear them and adjust the residual vector

        for( unsigned i = iBeg ; i < iEnd ; ++i ) {
          const unsigned       cnode = jacobian.graph.entries(i) ;
          const scalar_coord_type cc = node_coords(cnode,bc_plane);

          if ( ( cc <= bc_lower_limit ) || ( bc_upper_limit <= cc ) ) {
            jacobian.coeff(i) = 0 ;
          }
        }
      }
    }
  }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

/* A Cuda-specific specialization for the element computation functor. */
#if defined( __CUDACC__ )
// #include <NonlinearElement_Cuda.hpp>
#endif

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP */

