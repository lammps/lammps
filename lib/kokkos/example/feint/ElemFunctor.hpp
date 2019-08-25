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

#ifndef KOKKOS_EXAMPLE_FEINT_FUNCTORS_HPP
#define KOKKOS_EXAMPLE_FEINT_FUNCTORS_HPP

#include <cstdio>
#include <Kokkos_Core.hpp>
#include <BoxElemFixture.hpp>

namespace Kokkos {
namespace Example {

/** \brief  Numerically integrate a function on a finite element mesh and
 *          project the integrated values to nodes.
 */
template< class FixtureType ,
          class FunctionType ,
          bool PerformScatterAddWithAtomic >
struct FiniteElementIntegration ;

// Specialized for an 'Example::BoxElemFixture' finite element mesh
template< class Device , BoxElemPart::ElemOrder ElemOrder , class GridMap ,
          class FunctionType ,
          bool PerformScatterAddWithAtomic >
struct FiniteElementIntegration<
  Kokkos::Example::BoxElemFixture< Device , ElemOrder , GridMap > ,
  FunctionType ,
  PerformScatterAddWithAtomic >
{
  // Element mesh types:
  typedef Kokkos::Example::BoxElemFixture< Device , ElemOrder >
    BoxFixtureType ;

  typedef Kokkos::Example::HexElement_Data< BoxFixtureType::ElemNode >
    HexElemDataType ;

  enum { ElemNodeCount    = HexElemDataType::element_node_count  };
  enum { IntegrationCount = HexElemDataType::integration_count };
  enum { ValueCount       = FunctionType::value_count };

  // Dictionary of view types:
  typedef View<int*,                              Device> ElemErrorType ;
  typedef View<double*[ElemNodeCount][ValueCount],Device> ElemValueType ;
  typedef View<double*[ValueCount],               Device> NodeValueType ;

  // Data members for this Functor:
  const HexElemDataType  m_hex_elem_data ; ///< Master element
  const BoxFixtureType   m_box_fixture ;   ///< Unstructured mesh data
  const FunctionType     m_function ;      ///< Function to integrate
  const ElemErrorType    m_elem_error ;    ///< Flags for element errors
  const ElemValueType    m_elem_integral ; ///< Per-element quantities
  const NodeValueType    m_node_lumped ;   ///< Quantities lumped to nodes

  //----------------------------------------

  FiniteElementIntegration(
    const BoxFixtureType & box_fixture ,
    const FunctionType   & function )
    : m_hex_elem_data()
    , m_box_fixture( box_fixture ) // Shallow copy of the mesh fixture
    , m_function( function )
    , m_elem_error(    "elem_error"    , box_fixture.elem_count() )
    , m_elem_integral( "elem_integral" , box_fixture.elem_count() )
    , m_node_lumped(   "node_lumped"   , box_fixture.node_count() )
    {}

  //----------------------------------------
  // Device for parallel dispatch.
  typedef typename Device::execution_space execution_space;

  // Value type for global parallel reduction.
  struct value_type {
    double value[ ValueCount ]; ///< Integrated quantitie
    int    error ;              ///< Element inversion flag
  };

  //----------------------------------------
  // Transform element interpolation function gradients and
  // compute determinant of spatial jacobian.
  KOKKOS_INLINE_FUNCTION
  float transform_gradients(
    const float  grad[][  ElemNodeCount ] , // Gradient of bases master element
    const double coord[][ ElemNodeCount ] ,
          float  dpsi[][  ElemNodeCount ] ) const
  {
    enum { TensorDim = 9 };
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Temporary for jacobian accumulation is double for summation accuracy.
    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( int i = 0; i < ElemNodeCount ; ++i ) {
      J[j11] += grad[0][i] * coord[0][i] ;
      J[j12] += grad[0][i] * coord[1][i] ;
      J[j13] += grad[0][i] * coord[2][i] ;

      J[j21] += grad[1][i] * coord[0][i] ;
      J[j22] += grad[1][i] * coord[1][i] ;
      J[j23] += grad[1][i] * coord[2][i] ;

      J[j31] += grad[2][i] * coord[0][i] ;
      J[j32] += grad[2][i] * coord[1][i] ;
      J[j33] += grad[2][i] * coord[2][i] ;
    }

    // Inverse jacobian, compute as double and store as float.
    float invJ[ TensorDim ] = {
      float( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      float( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      float( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      float( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      float( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      float( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      float( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      float( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      float( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const float detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    {
      const float detJinv = 1.0 / detJ ;
      for ( int i = 0 ; i < TensorDim ; ++i ) { invJ[i] *= detJinv ; }
    }

    // Transform gradients:
    for ( int i = 0; i < ElemNodeCount ; ++i ) {
      dpsi[0][i] = grad[0][i] * invJ[j11] +
                   grad[1][i] * invJ[j12] +
                   grad[2][i] * invJ[j13];
      dpsi[1][i] = grad[0][i] * invJ[j21] +
                   grad[1][i] * invJ[j22] +
                   grad[2][i] * invJ[j23];
      dpsi[2][i] = grad[0][i] * invJ[j31] +
                   grad[1][i] * invJ[j32] +
                   grad[2][i] * invJ[j33];
    }

    return detJ ;
  }

  // Functor's function called for each element in the mesh
  // to numerically integrate the function and add element quantities
  // to the global integral.
  KOKKOS_INLINE_FUNCTION
  void operator()( const int ielem , value_type & update ) const
  {
    // Local temporaries for gathering nodal data.
    double node_coord[3][ ElemNodeCount ];

    int inode[ ElemNodeCount ] ;

    // Gather indices of element's node from global memory to local memory.
    for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
      inode[i] = m_box_fixture.elem_node( ielem , i );
    }

    // Gather coordinates of element's nodes from global memory to local memory.
    for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
      node_coord[0][i] = m_box_fixture.node_coord( inode[i] , 0 );
      node_coord[1][i] = m_box_fixture.node_coord( inode[i] , 1 );
      node_coord[2][i] = m_box_fixture.node_coord( inode[i] , 2 );
    }

    // Local temporary to accumulate numerical integration
    // of vector valued function.
    double accum[ ValueCount ];

    for ( int j = 0 ; j < ValueCount ; ++j ) { accum[j] = 0 ; }

    int error = 0 ;

    // Numerical integration loop for this element:
    for ( int k = 0 ; k < IntegrationCount ; ++k ) {

      // Integration point in space as interpolated from nodal coordinates:
      double point[3] = { 0 , 0 , 0 };
      for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
        point[0] += node_coord[0][i] * m_hex_elem_data.values[k][i] ;
        point[1] += node_coord[1][i] * m_hex_elem_data.values[k][i] ;
        point[2] += node_coord[2][i] * m_hex_elem_data.values[k][i] ;
      }

      // Example function vector value at cubature point:
      double val_at_pt[ ValueCount ];
      m_function( point , val_at_pt );

      // Temporary array for transformed element basis functions' gradient.
      // Not used in this example, but computed anyway by the more general
      // deformation function.
      float dpsi[3][ ElemNodeCount ];

      // Compute deformation jacobian, transform basis function gradient,
      // and return determinant of deformation jacobian.
      float detJ = transform_gradients( m_hex_elem_data.gradients[k] ,
                                        node_coord , dpsi );

      // Check for inverted spatial jacobian
      if ( detJ <= 0 ) { error = 1 ; detJ = 0 ; }

      // Integration weight.
      const float w = m_hex_elem_data.weights[k] * detJ ;

      // Cubature of function.
      for ( int j = 0 ; j < ValueCount ; ++j ) {
        accum[j] += val_at_pt[j] * w ;
      }
    }

    m_elem_error(ielem) = error ;


    // Element contribution to global integral:

    if ( error ) { update.error = 1 ; }

    for ( int j = 0 ; j < ValueCount ; ++j ) { update.value[j] += accum[j] ; }

    // Element-node quantity for lumping to nodes:
    for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
      for ( int j = 0 ; j < ValueCount ; ++j ) {
        // Save element's integral apportionment to nodes to global memory
        m_elem_integral( ielem , i , j ) = accum[j] / ElemNodeCount ;
      }
    }

    if ( PerformScatterAddWithAtomic ) {
      // Option to immediately scatter-add the integrated quantities to nodes.
      // This is a race condition as two or more threads could attempt
      // concurrent update of nodal values.  The atomic_fetch_add (+=)
      // function guarantees that the summation will occur correctly;
      // however, there can be no guarantee for the order of summation.
      // Due to non-associativity of floating point arithmetic the result
      // is non-deterministic within bounds of floating point round-off.

      for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
        for ( int j = 0 ; j < ValueCount ; ++j ) {
          Kokkos::atomic_fetch_add( & m_node_lumped( inode[i] , j ) ,
                                    m_elem_integral( ielem , i , j ) );
        }
      }
    }
  }
  //--------------------------------------------------------------------------

  // Initialization of the global reduction value.
  KOKKOS_INLINE_FUNCTION
  void init( value_type & update ) const
  {
    for ( int j = 0 ; j < ValueCount ; ++j ) update.value[j] = 0 ;
    update.error = 0 ;
  }

  // Join two contributions to global reduction value.
  KOKKOS_INLINE_FUNCTION
  void join( volatile       value_type & update ,
             volatile const value_type & input ) const
  {
    for ( int j = 0 ; j < ValueCount ; ++j ) update.value[j] += input.value[j] ;
    if ( input.error ) update.error = 1 ;
  }
};

} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template< class ViewElemNode ,
          class ViewNodeScan ,
          class ViewNodeElem >
void map_node_to_elem( const ViewElemNode & elem_node ,
                       const ViewNodeScan & node_scan ,
                       const ViewNodeElem & node_elem );

/** \brief  Functor to gather-sum elements' per-node quantities
 *          to element nodes.  Gather-sum is thread safe and
 *          does not require atomic updates.
 */
template< class ViewNodeValue ,
          class ViewElemValue ,
          bool  AlreadyUsedAtomic >
struct LumpElemToNode {

  typedef typename ViewElemValue::execution_space execution_space ;

  // In this example we know that the ViewElemValue
  // array specification is < double*[nNode][nValue] >

  enum { value_count = ViewElemValue::dimension::N2 };

  ViewNodeValue             m_node_value ; ///< Integrated values at nodes
  ViewElemValue             m_elem_value ; ///< Values apportioned to nodes
  View<int*,   execution_space> m_node_scan ;  ///< Offsets for nodes->element
  View<int*[2],execution_space> m_node_elem ;  ///< Node->element connectivity

  // Only allocate node->element connectivity if have
  // not already used atomic updates for the nodes.
  template< class ViewElemNode >
  LumpElemToNode( const ViewNodeValue & node_value ,
                  const ViewElemValue & elem_value ,
                  const ViewElemNode  & elem_node )
    : m_node_value( node_value )
    , m_elem_value( elem_value )
    , m_node_scan( "node_scan" ,
                   AlreadyUsedAtomic ? 0 : node_value.extent(0) + 1 )
    , m_node_elem( "node_elem" ,
                   AlreadyUsedAtomic ? 0 : elem_node.extent(0) *
                                           elem_node.extent(1) )
    {
      if ( ! AlreadyUsedAtomic ) {
        map_node_to_elem( elem_node , m_node_scan , m_node_elem );
      }
    }

  //----------------------------------------

  struct value_type { double value[ value_count ]; };

  KOKKOS_INLINE_FUNCTION
  void operator()( const int inode , value_type & update ) const
  {
    if ( ! AlreadyUsedAtomic ) {
      // Sum element quantities to a local variable.
      value_type local ;
      for ( int j = 0 ; j < value_count ; ++j ) { local.value[j] = 0 ; }

      {
        // nodes' element ids span [i,end)
        int i = m_node_scan(inode);
        const int end = m_node_scan(inode+1);

        for ( ; i < end ; ++i ) {
          // element #ielem , local node #ielem_node is this node:
          const int ielem      = m_node_elem(i,0);
          const int ielem_node = m_node_elem(i,1);
          // Sum the vector-values quantity
          for ( int j = 0 ; j < value_count ; ++j ) {
            local.value[j] += m_elem_value( ielem , ielem_node , j );
          }
        }
      }

      // Assign nodal quantity (no race condition).
      // Sum global value.
      for ( int j = 0 ; j < value_count ; ++j ) {
        m_node_value( inode , j ) = local.value[j] ;
        update.value[j] += local.value[j] ;
      }
    }
    else {
      // Already used atomic update of the nodal quantity,
      // query and sum the value.
      for ( int j = 0 ; j < value_count ; ++j ) {
        update.value[j] += m_node_value( inode , j );
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & update ) const
    { for ( int j = 0 ; j < value_count ; ++j ) { update.value[j] = 0 ; } }

  KOKKOS_INLINE_FUNCTION
  void join( volatile       value_type & update ,
             volatile const value_type & input ) const
    {
      for ( int j = 0 ; j < value_count ; ++j ) {
        update.value[j] += input.value[j] ;
      }
    }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class ViewElemNode ,
          class ViewNodeScan ,
          class ViewNodeElem >
void map_node_to_elem( const ViewElemNode & elem_node ,
                       const ViewNodeScan & node_scan ,
                       const ViewNodeElem & node_elem )
{
  typedef typename ViewElemNode::host_mirror_space host_mirror_space ;

  const typename ViewElemNode::HostMirror host_elem_node =
    Kokkos::create_mirror_view(elem_node);

  const typename ViewNodeScan::HostMirror host_node_scan =
    Kokkos::create_mirror_view(node_scan);

  const typename ViewNodeElem::HostMirror host_node_elem =
    Kokkos::create_mirror_view(node_elem);

  const int elem_count      = host_elem_node.extent(0);
  const int elem_node_count = host_elem_node.extent(1);
  const int node_count      = host_node_scan.extent(0) - 1 ;

  const View<int*, host_mirror_space >
    node_elem_count( "node_elem_count" , node_count );

  Kokkos::deep_copy( host_elem_node , elem_node );

  for ( int i = 0 ; i < elem_count ; ++i ) {
    for ( int j = 0 ; j < elem_node_count ; ++j ) {
      ++node_elem_count( host_elem_node(i,j) );
    }
  }

  for ( int i = 0 ; i < node_count ; ++i ) {
    host_node_scan(i+1) += host_node_scan(i) + node_elem_count(i);
    node_elem_count(i) = 0 ;
  }

  for ( int i = 0 ; i < elem_count ; ++i ) {
    for ( int j = 0 ; j < elem_node_count ; ++j ) {
      const int inode  = host_elem_node(i,j);
      const int offset = host_node_scan(inode) + node_elem_count(inode);

      host_node_elem( offset , 0 ) = i ;
      host_node_elem( offset , 1 ) = j ;

      ++node_elem_count(inode);
    }
  }

  Kokkos::deep_copy( node_scan , host_node_scan );
  Kokkos::deep_copy( node_elem , host_node_elem );
}

} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_EXAMPLE_FEINT_FUNCTORS_HPP */

