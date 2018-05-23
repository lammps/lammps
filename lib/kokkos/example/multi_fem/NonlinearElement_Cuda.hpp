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

#include <cstdio>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <HexElement.hpp>
#include <FEMesh.hpp>

namespace HybridFEM {
namespace Nonlinear {

template< class MeshType , typename ScalarType > struct ElementComputation ;

//----------------------------------------------------------------------------

template<>
struct ElementComputation< FEMesh< double , 27 , Kokkos::Cuda > , double >
{
  typedef Kokkos::Cuda    execution_space ;

  static const unsigned ElementNodeCount = 27 ;

  typedef HexElement_Data< ElementNodeCount >                element_data_type ;
  typedef FEMesh< double , ElementNodeCount , execution_space >  mesh_type ;

  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned FunctionCount    = element_data_type::function_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;

  typedef Kokkos::View< double[][FunctionCount][FunctionCount] , execution_space > elem_matrices_type ;
  typedef Kokkos::View< double[][FunctionCount] , execution_space > elem_vectors_type ;
  typedef Kokkos::View< double[] , execution_space > value_vector_type ;

private:

  const element_data_type                       elem_data ;
  const typename mesh_type::elem_node_ids_type  elem_node_ids ;
  const typename mesh_type::node_coords_type    node_coords ;
  const value_vector_type                       nodal_values ;
  const elem_matrices_type                      element_matrices ;
  const elem_vectors_type                       element_vectors ;
  const float                                   coeff_K ;
  const unsigned                                elem_count ;
        unsigned                                invJacIndex[9][4] ;

  static const unsigned j11 = 0 , j12 = 1 , j13 = 2 ,
                        j21 = 3 , j22 = 4 , j23 = 5 ,
                        j31 = 6 , j32 = 7 , j33 = 8 ;

  // Can only handle up to 16 warps:
  static const unsigned BlockDimX = 32 ;
  static const unsigned BlockDimY = 7 ;

  struct WorkSpace {
    double sum[ BlockDimY ][ BlockDimX ];

    double  value_at_integ[ IntegrationCount ];
    double  gradx_at_integ[ IntegrationCount ];
    double  grady_at_integ[ IntegrationCount ];
    double  gradz_at_integ[ IntegrationCount ];

    float  spaceJac[    BlockDimY ][ 9 ];
    float  spaceInvJac[ BlockDimY ][ 9 ];

    float  detJweight[ IntegrationCount ];

    float  dpsidx[ FunctionCount ][ IntegrationCount ];
    float  dpsidy[ FunctionCount ][ IntegrationCount ];
    float  dpsidz[ FunctionCount ][ IntegrationCount ];
  };

public:

  ElementComputation ( const mesh_type          & arg_mesh ,
                       const elem_matrices_type & arg_element_matrices ,
                       const elem_vectors_type  & arg_element_vectors ,
                       const value_vector_type  & arg_nodal_values ,
                       const float                arg_coeff_K )
  : elem_data()
  , elem_node_ids(    arg_mesh.elem_node_ids )
  , node_coords(      arg_mesh.node_coords )
  , nodal_values(     arg_nodal_values )
  , element_matrices( arg_element_matrices )
  , element_vectors(  arg_element_vectors )
  , coeff_K(          arg_coeff_K )
  , elem_count(       arg_mesh.elem_node_ids.dimension_0() )
  {
    const unsigned jInvJ[9][4] =
     { { j22 , j33 , j23 , j32 } ,
       { j13 , j32 , j12 , j33 } ,
       { j12 , j23 , j13 , j22 } ,

       { j23 , j31 , j21 , j33 } ,
       { j11 , j33 , j13 , j31 } ,
       { j13 , j21 , j11 , j23 } ,

       { j21 , j32 , j22 , j31 } ,
       { j12 , j31 , j11 , j32 } ,
       { j11 , j22 , j12 , j21 } };

    for ( unsigned i = 0 ; i < 9 ; ++i ) {
    for ( unsigned j = 0 ; j < 4 ; ++j ) {
      invJacIndex[i][j] = jInvJ[i][j] ;
    }
    }

    const unsigned shmem = sizeof(WorkSpace);
    const unsigned grid_max = 65535 ;
    const unsigned grid_count = std::min( grid_max , elem_count );

    // For compute capability 2.x up to 1024 threads per block
    const dim3 block( BlockDimX , BlockDimY , 1 );
    const dim3 grid( grid_count , 1 , 1 );

    Kokkos::Impl::CudaParallelLaunch< ElementComputation >( *this , grid , block , shmem );
  }

public:

  //------------------------------------
  // Sum among the threadIdx.x

  template< typename Type >
  __device__ inline static
  void sum_x( Type & result , const double value )
  {
    extern __shared__ WorkSpace work_data[] ;

    volatile double * const base_sum =
      & work_data->sum[ threadIdx.y ][ threadIdx.x ] ;

    base_sum[ 0] = value ;

    if ( threadIdx.x < 16 ) {
      base_sum[0] += base_sum[16];
      base_sum[0] += base_sum[ 8];
      base_sum[0] += base_sum[ 4];
      base_sum[0] += base_sum[ 2];
      base_sum[0] += base_sum[ 1];
    }

    if ( 0 == threadIdx.x ) {
      result = base_sum[0] ;
    }
  }

  __device__ inline static
  void sum_x_clear()
  {
    extern __shared__ WorkSpace work_data[] ;

    work_data->sum[ threadIdx.y ][ threadIdx.x ] = 0 ;
  }

  //------------------------------------
  //------------------------------------

  __device__ inline
  void evaluateFunctions( const unsigned ielem ) const
  {
    extern __shared__ WorkSpace work_data[] ;

    // Each warp (threadIdx.y) computes an integration point
    // Each thread is responsible for a node / function.

    const unsigned iFunc = threadIdx.x ;
    const bool     hasFunc = iFunc < FunctionCount ;

    //------------------------------------
    // Each warp gathers a different variable into 'elem_mat' shared memory.

    if ( hasFunc ) {

      const unsigned node = elem_node_ids( ielem , iFunc );

      for ( unsigned iy = threadIdx.y ; iy < 4 ; iy += blockDim.y ) {
      switch( iy ) {
      case 0 : work_data->sum[0][iFunc] = node_coords(node,0); break ;
      case 1 : work_data->sum[1][iFunc] = node_coords(node,1); break ;
      case 2 : work_data->sum[2][iFunc] = node_coords(node,2); break ;
      case 3 : work_data->sum[3][iFunc] = nodal_values(node); break ;
      default: break ;
      }
      }
    }

    __syncthreads(); // Wait for all warps to finish gathering

    // now get local 'const' copies in register space:

    const double x       = work_data->sum[0][ iFunc ];
    const double y       = work_data->sum[1][ iFunc ];
    const double z       = work_data->sum[2][ iFunc ];
    const double dof_val = work_data->sum[3][ iFunc ];

    __syncthreads(); // Wait for all warps to finish extracting

    sum_x_clear(); // Make sure summation scratch is zero

    //------------------------------------
    // Each warp is now on its own computing an integration point
    // so no further explicit synchronizations are required.

    if ( hasFunc ) {

      float * const J    = work_data->spaceJac[    threadIdx.y ];
      float * const invJ = work_data->spaceInvJac[ threadIdx.y ];

      for ( unsigned iInt = threadIdx.y ;
                     iInt < IntegrationCount ; iInt += blockDim.y ) {

        const float val = elem_data.values[iInt][iFunc] ;
        const float gx  = elem_data.gradients[iInt][0][iFunc] ;
        const float gy  = elem_data.gradients[iInt][1][iFunc] ;
        const float gz  = elem_data.gradients[iInt][2][iFunc] ;

        sum_x( J[j11], gx * x );
        sum_x( J[j12], gx * y );
        sum_x( J[j13], gx * z );

        sum_x( J[j21], gy * x );
        sum_x( J[j22], gy * y );
        sum_x( J[j23], gy * z );

        sum_x( J[j31], gz * x );
        sum_x( J[j32], gz * y );
        sum_x( J[j33], gz * z );

        // Inverse jacobian, only enough parallel work for 9 threads in the warp

        if ( iFunc < TensorDim ) {

          invJ[ iFunc ] =
            J[ invJacIndex[iFunc][0] ] * J[ invJacIndex[iFunc][1] ] -
            J[ invJacIndex[iFunc][2] ] * J[ invJacIndex[iFunc][3] ] ;

          // Let all threads in the warp compute determinant into a register

          const float detJ = J[j11] * invJ[j11] +
                             J[j21] * invJ[j12] +
                             J[j31] * invJ[j13] ;

          invJ[ iFunc ] /= detJ ;

          if ( 0 == iFunc ) {
            work_data->detJweight[ iInt ] = detJ * elem_data.weights[ iInt ] ;
          }
        }

        // Transform bases gradients and compute value and gradient

        const float dx = gx * invJ[j11] + gy * invJ[j12] + gz * invJ[j13];
        const float dy = gx * invJ[j21] + gy * invJ[j22] + gz * invJ[j23];
        const float dz = gx * invJ[j31] + gy * invJ[j32] + gz * invJ[j33];

        work_data->dpsidx[iFunc][iInt] = dx ;
        work_data->dpsidy[iFunc][iInt] = dy ;
        work_data->dpsidz[iFunc][iInt] = dz ;

        sum_x( work_data->gradx_at_integ[iInt] , dof_val * dx );
        sum_x( work_data->grady_at_integ[iInt] , dof_val * dy );
        sum_x( work_data->gradz_at_integ[iInt] , dof_val * dz );
        sum_x( work_data->value_at_integ[iInt] , dof_val * val );
      }
    }

    __syncthreads(); // All shared data must be populated at return.
  }

  __device__ inline
  void contributeResidualJacobian( const unsigned ielem ) const
  {
    extern __shared__ WorkSpace work_data[] ;

    sum_x_clear(); // Make sure summation scratch is zero

    // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$
    // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$

    const unsigned iInt = threadIdx.x ;

    if ( iInt < IntegrationCount ) {

      const double value_at_integ = work_data->value_at_integ[ iInt ] ;
      const double gradx_at_integ = work_data->gradx_at_integ[ iInt ] ;
      const double grady_at_integ = work_data->grady_at_integ[ iInt ] ;
      const double gradz_at_integ = work_data->gradz_at_integ[ iInt ] ;

      const float detJweight     = work_data->detJweight[ iInt ] ;
      const float coeff_K_detJweight = coeff_K * detJweight ;

      for ( unsigned iRow = threadIdx.y ;
                     iRow < FunctionCount ; iRow += blockDim.y ) {

        const float value_row  = elem_data.values[ iInt ][ iRow ] * detJweight ;
        const float dpsidx_row = work_data->dpsidx[ iRow ][ iInt ] * coeff_K_detJweight ;
        const float dpsidy_row = work_data->dpsidy[ iRow ][ iInt ] * coeff_K_detJweight ;
        const float dpsidz_row = work_data->dpsidz[ iRow ][ iInt ] * coeff_K_detJweight ;

        const double res_del = dpsidx_row * gradx_at_integ +
                               dpsidy_row * grady_at_integ +
                               dpsidz_row * gradz_at_integ ;

        const double res_val = value_at_integ * value_at_integ * value_row ;
        const double jac_val_row = 2 * value_at_integ * value_row ;

        sum_x( element_vectors( ielem , iRow ) , res_del + res_val );

        for ( unsigned iCol = 0 ; iCol < FunctionCount ; ++iCol ) {

          const float jac_del =
            dpsidx_row * work_data->dpsidx[iCol][iInt] +
            dpsidy_row * work_data->dpsidy[iCol][iInt] +
            dpsidz_row * work_data->dpsidz[iCol][iInt] ;

          const double jac_val =
            jac_val_row * elem_data.values[ iInt ][ iCol ] ;

          sum_x( element_matrices( ielem , iRow , iCol ) , jac_del + jac_val );
        }
      }
    }

    __syncthreads(); // All warps finish before refilling shared data
  }

  __device__ inline
  void operator()(void) const
  {
    extern __shared__ WorkSpace work_data[] ;

    for ( unsigned ielem = blockIdx.x ; ielem < elem_count ; ielem += gridDim.x ) {

      evaluateFunctions( ielem );

      contributeResidualJacobian( ielem );
    }
  }

}; /* ElementComputation */

} /* namespace Nonlinear */
} /* namespace HybridFEM */

