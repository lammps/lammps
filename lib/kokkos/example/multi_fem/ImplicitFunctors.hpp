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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

namespace HybridFEM {
namespace Implicit {

//----------------------------------------------------------------------------

template< typename Scalar , unsigned Dim , unsigned N >
struct TensorIntegration ;

template<typename Scalar >
struct TensorIntegration<Scalar,1,1> {
  Scalar pts[1] ;
  Scalar wts[1] ;

  TensorIntegration() { pts[0] = 0 ; wts[0] = 2 ; }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,2>
{
  Scalar pts[2] ;
  Scalar wts[2] ;

  TensorIntegration()
  {
    const Scalar x2 = 0.577350269 ;
    pts[0] = -x2; wts[0] = 1.0;
    pts[1] =  x2; wts[1] = 1.0;
  }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,3>
{
  Scalar pts[3] ;
  Scalar wts[3] ;

  TensorIntegration()
  {
    const Scalar x3 = 0.774596669 ;
    const Scalar w1 = 0.555555556 ;
    const Scalar w2 = 0.888888889 ;
    pts[0] =  -x3 ;  wts[0] = w1 ;
    pts[1] =    0 ;  wts[1] = w2 ;
    pts[2] =   x3 ;  wts[2] = w1 ;
  }
};

template< typename Scalar , unsigned Order >
struct TensorIntegration<Scalar,3,Order>
{
  static const unsigned N = Order * Order * Order ;

  Scalar pts[N][3] ;
  Scalar wts[N];

  TensorIntegration()
  {
    TensorIntegration<Scalar,1,Order> oneD ;

    unsigned n = 0 ;
    for ( unsigned k = 0 ; k < Order ; ++k ) {
    for ( unsigned j = 0 ; j < Order ; ++j ) {
    for ( unsigned i = 0 ; i < Order ; ++i , ++n ) {
      pts[n][0] = oneD.pts[i] ;
      pts[n][1] = oneD.pts[j] ;
      pts[n][2] = oneD.pts[k] ;
      wts[n] = oneD.wts[i] * oneD.wts[j] * oneD.wts[k] ;
    }}}
  }
};

//----------------------------------------------------------------------------

template< typename Scalar >
struct ShapeFunctionEvaluation {

  static const unsigned FunctionCount = 8 ;
  static const unsigned SpatialDimension = 3 ;
  static const unsigned IntegrationOrder = 2 ;

  typedef TensorIntegration< Scalar , SpatialDimension , IntegrationOrder > 
    TensorIntegrationType ;

  static const unsigned PointCount = TensorIntegrationType::N ;

  Scalar value   [ PointCount ][ FunctionCount ] ;
  Scalar gradient[ PointCount ][ FunctionCount * SpatialDimension ];
  Scalar weight  [ PointCount ];

  ShapeFunctionEvaluation()
  {
    const TensorIntegration< Scalar , SpatialDimension , IntegrationOrder > 
      integration ;

    const Scalar ONE8TH = 0.125 ;

    for ( unsigned i = 0 ; i < PointCount ; ++i ) {

      const Scalar u = 1.0 - integration.pts[i][0];
      const Scalar v = 1.0 - integration.pts[i][1];
      const Scalar w = 1.0 - integration.pts[i][2];

      const Scalar up1 = 1.0 + integration.pts[i][0];
      const Scalar vp1 = 1.0 + integration.pts[i][1];
      const Scalar wp1 = 1.0 + integration.pts[i][2];

      weight[i] = integration.wts[i] ;

      // Vaues:
      value[i][0] = ONE8TH *   u *   v *  w ;
      value[i][1] = ONE8TH * up1 *   v *  w ;
      value[i][2] = ONE8TH * up1 * vp1 *  w ;
      value[i][3] = ONE8TH *   u * vp1 *  w ;

      value[i][4] = ONE8TH *   u *   v *  wp1 ;
      value[i][5] = ONE8TH * up1 *   v *  wp1 ;
      value[i][6] = ONE8TH * up1 * vp1 *  wp1 ;
      value[i][7] = ONE8TH *   u * vp1 *  wp1 ;

      //fn 0 = u * v * w
      gradient[i][ 0] = ONE8TH * -1  *  v  *  w  ;
      gradient[i][ 1] = ONE8TH *  u  * -1  *  w  ;
      gradient[i][ 2] = ONE8TH *  u  *  v  * -1  ;

      //fn 1 = up1 * v * w
      gradient[i][ 3] = ONE8TH *  1  *  v  *  w  ;
      gradient[i][ 4] = ONE8TH * up1 * -1  *  w  ;
      gradient[i][ 5] = ONE8TH * up1 *  v  * -1  ;

      //fn 2 = up1 * vp1 * w
      gradient[i][ 6] = ONE8TH *  1  * vp1 *  w ;
      gradient[i][ 7] = ONE8TH * up1 *  1  *  w ;
      gradient[i][ 8] = ONE8TH * up1 * vp1 * -1 ;

      //fn 3 = u * vp1 * w
      gradient[i][ 9] = ONE8TH * -1 * vp1 *  w ;
      gradient[i][10] = ONE8TH *  u *  1  *  w ;
      gradient[i][11] = ONE8TH *  u * vp1 * -1 ;

      //fn 4 = u * v * wp1
      gradient[i][12] = ONE8TH * -1  *  v  * wp1 ;
      gradient[i][13] = ONE8TH *  u  * -1  * wp1 ;
      gradient[i][14] = ONE8TH *  u  *  v  *  1  ;

      //fn 5 = up1 * v * wp1
      gradient[i][15] = ONE8TH *  1  *  v  * wp1 ;
      gradient[i][16] = ONE8TH * up1 * -1  * wp1 ;
      gradient[i][17] = ONE8TH * up1 *  v  *  1  ;

      //fn 6 = up1 * vp1 * wp1
      gradient[i][18] = ONE8TH *  1  * vp1 * wp1 ;
      gradient[i][19] = ONE8TH * up1 *  1  * wp1 ;
      gradient[i][20] = ONE8TH * up1 * vp1 *  1 ;

      //fn 7 = u * vp1 * wp1
      gradient[i][21] = ONE8TH * -1 * vp1 * wp1 ;
      gradient[i][22] = ONE8TH *  u *  1  * wp1 ;
      gradient[i][23] = ONE8TH *  u * vp1 *  1 ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename ScalarType , typename ScalarCoordType , class DeviceType >
struct ElementComputation
{
  typedef DeviceType     execution_space;
  typedef ScalarType              scalar_type ;
  typedef typename execution_space::size_type  size_type ;

  static const size_type ElementNodeCount = 8 ;

  typedef FEMesh< ScalarCoordType , ElementNodeCount , execution_space > mesh_type ;
  typedef Kokkos::View< scalar_type[][ElementNodeCount][ElementNodeCount] , execution_space > elem_matrices_type ;
  typedef Kokkos::View< scalar_type[][ElementNodeCount] , execution_space > elem_vectors_type ;

  typedef ShapeFunctionEvaluation< scalar_type > shape_function_data ;

  static const unsigned SpatialDim    = shape_function_data::SpatialDimension ;
  static const unsigned FunctionCount = shape_function_data::FunctionCount ;

private:

  const shape_function_data               shape_eval ;
  typename mesh_type::elem_node_ids_type  elem_node_ids ;
  typename mesh_type::node_coords_type    node_coords ;
  elem_matrices_type                      element_matrices ;
  elem_vectors_type                       element_vectors ;
  scalar_type                             coeff_K ;
  scalar_type                             coeff_Q ;

  ElementComputation( const mesh_type   & arg_mesh ,
                      const elem_matrices_type  & arg_element_matrices , 
                      const elem_vectors_type   & arg_element_vectors ,
                      const scalar_type   arg_coeff_K ,
                      const scalar_type   arg_coeff_Q )
  : shape_eval()
  , elem_node_ids( arg_mesh.elem_node_ids )
  , node_coords(   arg_mesh.node_coords )
  , element_matrices( arg_element_matrices )
  , element_vectors( arg_element_vectors )
  , coeff_K( arg_coeff_K )
  , coeff_Q( arg_coeff_Q )
  {}

public:

  static void apply( const mesh_type  & mesh ,
                     const elem_matrices_type & elem_matrices ,
                     const elem_vectors_type  & elem_vectors ,
                     const scalar_type  elem_coeff_K ,
                     const scalar_type  elem_coeff_Q )
  {
    ElementComputation comp( mesh , elem_matrices , elem_vectors , elem_coeff_K , elem_coeff_Q );
    const size_t elem_count = mesh.elem_node_ids.dimension_0();

    parallel_for( elem_count , comp );
  }

  //------------------------------------

  static const unsigned FLOPS_jacobian =
    FunctionCount * SpatialDim * SpatialDim * 2 ;

  KOKKOS_INLINE_FUNCTION
  void jacobian( const ScalarCoordType * x, 
                 const ScalarCoordType * y, 
                 const ScalarCoordType * z, 
                 const scalar_type * grad_vals, 
                 scalar_type * J) const
  {
    int i_grad = 0 ;

    for( unsigned i = 0; i < ElementNodeCount ; ++i , i_grad += SpatialDim ) {
      const scalar_type g0 = grad_vals[ i_grad ];
      const scalar_type g1 = grad_vals[ i_grad + 1 ];
      const scalar_type g2 = grad_vals[ i_grad + 2 ];
      const scalar_type x0 = x[i] ;
      const scalar_type x1 = y[i] ;
      const scalar_type x2 = z[i] ;

      J[0] += g0 * x0 ;
      J[1] += g0 * x1 ;
      J[2] += g0 * x2 ;

      J[3] += g1 * x0 ;
      J[4] += g1 * x1 ;
      J[5] += g1 * x2 ;

      J[6] += g2 * x0 ;
      J[7] += g2 * x1 ;
      J[8] += g2 * x2 ;
    }
  }

  //------------------------------------

  static const unsigned FLOPS_inverse_and_det = 46 ;

  KOKKOS_INLINE_FUNCTION
  scalar_type inverse_and_determinant3x3( scalar_type * const J ) const
  {
    const scalar_type J00 = J[0];
    const scalar_type J01 = J[1];
    const scalar_type J02 = J[2];

    const scalar_type J10 = J[3];
    const scalar_type J11 = J[4];
    const scalar_type J12 = J[5];

    const scalar_type J20 = J[6];
    const scalar_type J21 = J[7];
    const scalar_type J22 = J[8];

    const scalar_type term0 = J22*J11 - J21*J12;
    const scalar_type term1 = J22*J01 - J21*J02;
    const scalar_type term2 = J12*J01 - J11*J02;

    const scalar_type detJ = J00*term0 - J10*term1 + J20*term2;
    const scalar_type inv_detJ = 1.0/detJ;

    J[0] =  term0*inv_detJ;
    J[1] = -term1*inv_detJ;
    J[2] =  term2*inv_detJ;

    J[3] = -(J22*J10 - J20*J12)*inv_detJ;
    J[4] =  (J22*J00 - J20*J02)*inv_detJ;
    J[5] = -(J12*J00 - J10*J02)*inv_detJ;

    J[6] =  (J21*J10 - J20*J11)*inv_detJ;
    J[7] = -(J21*J00 - J20*J01)*inv_detJ;
    J[8] =  (J11*J00 - J10*J01)*inv_detJ;

    return detJ ;
  }

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  void matTransMat3x3_X_3xn( const scalar_type * A, int n,
                             const scalar_type * B,
                             scalar_type * C ) const
  {
    //A is 3x3, B is 3xn. So C is also 3xn.
    //A,B,C are all assumed to be ordered such that columns are contiguous.

    scalar_type * Cj = C;
    const scalar_type * Bj = B;

    for(int j=0; j<n; ++j) {
      Cj[0] = A[0]*Bj[0] + A[1]*Bj[1] + A[2]*Bj[2];
      Cj[1] = A[3]*Bj[0] + A[4]*Bj[1] + A[5]*Bj[2];
      Cj[2] = A[6]*Bj[0] + A[7]*Bj[1] + A[8]*Bj[2];
      Bj += 3;
      Cj += 3;
    }

  }
  //------------------------------------

  static const unsigned FLOPS_contributeDiffusionMatrix = FunctionCount * ( 3 * 5 + FunctionCount * 7 ) ;

  KOKKOS_INLINE_FUNCTION
  void contributeDiffusionMatrix(
    const scalar_type weight ,
    const scalar_type grad_vals[] ,
    const scalar_type invJ[] ,
    scalar_type elem_mat[][8] ) const
  {
    scalar_type dpsidx[8], dpsidy[8], dpsidz[8];

    int i_grad = 0 ;
    for( unsigned i = 0; i < FunctionCount ; ++i , i_grad += 3 ) {
      const scalar_type g0 = grad_vals[i_grad+0];
      const scalar_type g1 = grad_vals[i_grad+1];
      const scalar_type g2 = grad_vals[i_grad+2];

      dpsidx[i] = g0 * invJ[0] + g1 * invJ[1] + g2 * invJ[2];
      dpsidy[i] = g0 * invJ[3] + g1 * invJ[4] + g2 * invJ[5];
      dpsidz[i] = g0 * invJ[6] + g1 * invJ[7] + g2 * invJ[8];
    }

    for( unsigned m = 0; m < FunctionCount; m++) {
      for( unsigned n = 0; n < FunctionCount; n++) {

        elem_mat[m][n] += weight * 
          ((dpsidx[m] * dpsidx[n]) + 
           (dpsidy[m] * dpsidy[n]) +
           (dpsidz[m] * dpsidz[n]));            
      }
    }
  }

  //------------------------------------

  static const unsigned FLOPS_contributeSourceVector = FunctionCount * 2 ;

  KOKKOS_INLINE_FUNCTION
  void contributeSourceVector( const scalar_type term ,
                               const scalar_type psi[] ,
                               scalar_type elem_vec[] ) const
  {
     for( unsigned i=0; i< FunctionCount ; ++i) {
       elem_vec[i] += psi[i] * term ;
     }
  }


  static const unsigned FLOPS_operator =
           shape_function_data::PointCount * ( 3
             + FLOPS_jacobian
             + FLOPS_inverse_and_det
             + FLOPS_contributeDiffusionMatrix
             + FLOPS_contributeSourceVector ) ;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ielem )const {

    scalar_type elem_vec[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
    scalar_type elem_mat[8][8] =
      { { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } };

    ScalarCoordType x[8], y[8], z[8];

    for ( int i = 0 ; i < 8 ; ++i ) {
      const int node_index = elem_node_ids( ielem , i );
      x[i] = node_coords( node_index , 0 );
      y[i] = node_coords( node_index , 1 );
      z[i] = node_coords( node_index , 2 );
    }

    // This loop could be parallelized; however,
    // it would require additional per-thread temporaries
    // of 'elem_vec' and 'elem_mat' which would
    // consume more local memory and have to be reduced.

    for ( unsigned i = 0 ; i < shape_function_data::PointCount ; ++i ) {

      scalar_type J[SpatialDim*SpatialDim] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

      jacobian( x, y, z, shape_eval.gradient[i] , J );

      // Overwrite J with its inverse to save scratch memory space.
      const scalar_type detJ_w   = shape_eval.weight[i] * inverse_and_determinant3x3(J);
      const scalar_type k_detJ_w = coeff_K * detJ_w ;
      const scalar_type Q_detJ_w = coeff_Q * detJ_w ;

      contributeDiffusionMatrix( k_detJ_w , shape_eval.gradient[i] , J , elem_mat );

      contributeSourceVector( Q_detJ_w , shape_eval.value[i] , elem_vec );
    }

    for( size_type i=0; i< ElementNodeCount ; ++i) {
      element_vectors(ielem, i) = elem_vec[i] ;
    }

    for( size_type i = 0; i < ElementNodeCount ; i++){
      for( size_type j = 0; j < ElementNodeCount ; j++){
        element_matrices(ielem, i, j) = elem_mat[i][j] ;
      }
    }
  }
}; /* ElementComputation */

//----------------------------------------------------------------------------

template< typename ScalarType , typename ScalarCoordType , class DeviceType >
struct DirichletBoundary
{
  typedef DeviceType     execution_space;
  typedef typename execution_space::size_type  size_type ;

  static const size_type ElementNodeCount = 8 ;

  typedef Kokkos::CrsMatrix< ScalarType , execution_space >    matrix_type ;
  typedef Kokkos::View< ScalarType[] , execution_space >  vector_type ;

  typedef FEMesh< ScalarCoordType , ElementNodeCount , execution_space > mesh_type ;

  typename mesh_type::node_coords_type node_coords ;
  matrix_type     matrix ;
  vector_type     rhs ;
  ScalarCoordType bc_lower_z ;
  ScalarCoordType bc_upper_z ;
  ScalarType      bc_lower_value ;
  ScalarType      bc_upper_value ;

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type inode ) const
  {
    //  Apply a dirichlet boundary condition to 'irow'
    //  to maintain the symmetry of the original 
    //  global stiffness matrix, zero out the columns
    //  that correspond to boundary conditions, and
    //  adjust the load vector accordingly

    const size_type iBeg = matrix.graph.row_map[inode];
    const size_type iEnd = matrix.graph.row_map[inode+1];

    const ScalarCoordType z = node_coords(inode,2);
    const bool bc_lower = z <= bc_lower_z ;
    const bool bc_upper = bc_upper_z <= z ;

    if ( bc_lower || bc_upper ) {
      const ScalarType bc_value = bc_lower ? bc_lower_value
                                           : bc_upper_value ;

      rhs(inode) = bc_value ; //  set the rhs vector

      //  zero each value on the row, and leave a one
      //  on the diagonal

      for( size_type i = iBeg ; i < iEnd ; i++) {
        matrix.coefficients(i) =
          (int) inode == matrix.graph.entries(i) ? 1 : 0 ;
      }
    }
    else {
      //  Find any columns that are boundary conditions.
      //  Clear them and adjust the load vector

      for( size_type i = iBeg ; i < iEnd ; i++ ) {
        const size_type cnode = matrix.graph.entries(i) ;

        const ScalarCoordType zc = node_coords(cnode,2);
        const bool c_bc_lower = zc <= bc_lower_z ;
        const bool c_bc_upper = bc_upper_z <= zc ;

        if ( c_bc_lower || c_bc_upper ) {

          const ScalarType c_bc_value = c_bc_lower ? bc_lower_value
                                                   : bc_upper_value ;

          rhs( inode ) -= c_bc_value * matrix.coefficients(i);

          matrix.coefficients(i) = 0 ;
        }
      }
    }
  }


  static void apply( const matrix_type & linsys_matrix ,
                     const vector_type & linsys_rhs ,
                     const mesh_type   & mesh ,
                     const ScalarCoordType  bc_lower_z ,
                     const ScalarCoordType  bc_upper_z ,
                     const ScalarType       bc_lower_value ,
                     const ScalarType       bc_upper_value )
  {
    const size_t row_count = linsys_matrix.graph.row_map.dimension_0() - 1 ;
    DirichletBoundary op ;
    op.node_coords    = mesh.node_coords ;
    op.matrix         = linsys_matrix ;
    op.rhs            = linsys_rhs ;
    op.bc_lower_z     = bc_lower_z ;
    op.bc_upper_z     = bc_upper_z ;
    op.bc_lower_value = bc_lower_value ;
    op.bc_upper_value = bc_upper_value ;
    parallel_for( row_count , op );
  }
};

//----------------------------------------------------------------------------

} /* namespace Implicit */
} /* namespace HybridFEM */

