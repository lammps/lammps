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

#ifndef KOKKOS_EXPLICITFUNCTORS_HPP
#define KOKKOS_EXPLICITFUNCTORS_HPP

#include <cmath>
#include <Kokkos_Core.hpp>
#include <FEMesh.hpp>

namespace Explicit {

template<typename Scalar , class Device >
struct Fields {

  static const int NumStates     = 2 ;
  static const int SpatialDim    = 3 ;
  static const int ElemNodeCount = 8 ;

  // Indices for full 3x3 tensor:

  static const int K_F_XX = 0 ;
  static const int K_F_YY = 1 ;
  static const int K_F_ZZ = 2 ;
  static const int K_F_XY = 3 ;
  static const int K_F_YZ = 4 ;
  static const int K_F_ZX = 5 ;
  static const int K_F_YX = 6 ;
  static const int K_F_ZY = 7 ;
  static const int K_F_XZ = 8 ;

  //  Indexes into a 3 by 3 symmetric tensor stored as a length 6 vector

  static const int K_S_XX = 0 ;
  static const int K_S_YY = 1 ;
  static const int K_S_ZZ = 2 ;
  static const int K_S_XY = 3 ;
  static const int K_S_YZ = 4 ;
  static const int K_S_ZX = 5 ;
  static const int K_S_YX = 3 ;
  static const int K_S_ZY = 4 ;
  static const int K_S_XZ = 5 ;

  //  Indexes into a 3 by 3 skew symmetric tensor stored as a length 3 vector

  static const int K_V_XY = 0 ;
  static const int K_V_YZ = 1 ;
  static const int K_V_ZX = 2 ;


  typedef Device                           execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef HybridFEM::FEMesh<double,ElemNodeCount,execution_space>  FEMesh ;

  typedef typename FEMesh::node_coords_type    node_coords_type ;
  typedef typename FEMesh::elem_node_ids_type  elem_node_ids_type ;
  typedef typename FEMesh::node_elem_ids_type  node_elem_ids_type ;
  typedef typename Kokkos::ParallelDataMap   parallel_data_map ;

  typedef Kokkos::View< double[][ SpatialDim ][ NumStates ] , execution_space > geom_state_array_type ;
  typedef Kokkos::View< Scalar[][ SpatialDim ] , execution_space > geom_array_type ;
  typedef Kokkos::View< Scalar[] ,               execution_space > array_type ;
  typedef Kokkos::View< Scalar ,                 execution_space >  scalar_type ;

  typedef Kokkos::View< Scalar[][  6 ] ,    execution_space >  elem_sym_tensor_type ;
  typedef Kokkos::View< Scalar[][  9 ] ,    execution_space >  elem_tensor_type ;
  typedef Kokkos::View< Scalar[][  9 ][ NumStates ] , execution_space >  elem_tensor_state_type ;
  typedef Kokkos::View< Scalar[][ SpatialDim ][ ElemNodeCount ] , execution_space > elem_node_geom_type ;

  // Parameters:
  const int num_nodes ;
  const int num_nodes_owned ;
  const int num_elements ;

  const Scalar  lin_bulk_visc;
  const Scalar  quad_bulk_visc;
  const Scalar  two_mu;
  const Scalar  bulk_modulus;
  const Scalar  density;

  // Mesh:
  const elem_node_ids_type  elem_node_connectivity ;
  const node_elem_ids_type  node_elem_connectivity ;
  const node_coords_type    model_coords ;

  // Compute:
  const scalar_type                dt ;
  const scalar_type                prev_dt ;
  const geom_state_array_type      displacement ;
  const geom_state_array_type      velocity ;
  const geom_array_type            acceleration ;
  const geom_array_type            internal_force ;
  const array_type                 nodal_mass ;
  const array_type                 elem_mass ;
  const array_type                 internal_energy ;
  const elem_sym_tensor_type       stress_new ;
  const elem_tensor_state_type     rotation ;
  const elem_node_geom_type        element_force ;
  const elem_tensor_type           vel_grad ;
  const elem_sym_tensor_type       stretch ;
  const elem_sym_tensor_type       rot_stretch ;

  Fields(
      const FEMesh & mesh,
      Scalar arg_lin_bulk_visc,
      Scalar arg_quad_bulk_visc,
      Scalar youngs_modulus,
      Scalar poissons_ratio,
      Scalar arg_density )
    : num_nodes(       mesh.parallel_data_map.count_owned +
                       mesh.parallel_data_map.count_receive )
    , num_nodes_owned( mesh.parallel_data_map.count_owned )
    , num_elements(    mesh.elem_node_ids.dimension_0() )
    , lin_bulk_visc(  arg_lin_bulk_visc )
    , quad_bulk_visc( arg_quad_bulk_visc )
    , two_mu(youngs_modulus/(1.0+poissons_ratio))
    , bulk_modulus(youngs_modulus/(3*(1.0-2.0*poissons_ratio)))
    , density(arg_density)

    // mesh

    , elem_node_connectivity( mesh.elem_node_ids ) // ( num_elements , ElemNodeCount )
    , node_elem_connectivity( mesh.node_elem_ids ) // ( num_nodes , ... )
    , model_coords(  mesh.node_coords )            // ( num_nodes , 3 )

    // compute with input/output

    , dt(              "dt" )
    , prev_dt(         "prev_dt" )
    , displacement(    "displacement" ,   num_nodes )
    , velocity(        "velocity" ,       num_nodes )
    , acceleration(    "acceleration" ,   num_nodes_owned )
    , internal_force(  "internal_force" , num_nodes_owned )
    , nodal_mass(      "nodal_mass" ,     num_nodes_owned )
    , elem_mass(       "elem_mass" ,       num_elements )
    , internal_energy( "internal_energy" , num_elements )
    , stress_new(      "stress_new" ,      num_elements )

    // temporary arrays

    , rotation(      "rotation" ,  num_elements )
    , element_force( "element_force" ,  num_elements )
    , vel_grad(      "vel_grad" , num_elements )
    , stretch(       "stretch" , num_elements )
    , rot_stretch(   "rot_stretch" , num_elements )
  { }
};


//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
KOKKOS_INLINE_FUNCTION
Scalar dot8( const Scalar * a , const Scalar * b )
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] +
         a[4] * b[4] + a[5] * b[5] + a[6] * b[6] + a[7] * b[7] ; }

template< typename Scalar , class DeviceType >
KOKKOS_INLINE_FUNCTION
void comp_grad( const Scalar * const x ,
                const Scalar * const y ,
                const Scalar * const z,
                Scalar * const grad_x ,
                Scalar * const grad_y ,
                Scalar * const grad_z )
{
  //  calc X difference vectors

  Scalar R42=(x[3] - x[1]);
  Scalar R52=(x[4] - x[1]);
  Scalar R54=(x[4] - x[3]);

  Scalar R63=(x[5] - x[2]);
  Scalar R83=(x[7] - x[2]);
  Scalar R86=(x[7] - x[5]);

  Scalar R31=(x[2] - x[0]);
  Scalar R61=(x[5] - x[0]);
  Scalar R74=(x[6] - x[3]);

  Scalar R72=(x[6] - x[1]);
  Scalar R75=(x[6] - x[4]);
  Scalar R81=(x[7] - x[0]);

  Scalar t1=(R63 + R54);
  Scalar t2=(R61 + R74);
  Scalar t3=(R72 + R81);

  Scalar t4 =(R86 + R42);
  Scalar t5 =(R83 + R52);
  Scalar t6 =(R75 + R31);

  //  Calculate Y gradient from X and Z data

  grad_y[0] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
  grad_y[1] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
  grad_y[2] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
  grad_y[3] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
  grad_y[4] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
  grad_y[5] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
  grad_y[6] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
  grad_y[7] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

  //   calc Z difference vectors

  R42=(z[3] - z[1]);
  R52=(z[4] - z[1]);
  R54=(z[4] - z[3]);

  R63=(z[5] - z[2]);
  R83=(z[7] - z[2]);
  R86=(z[7] - z[5]);

  R31=(z[2] - z[0]);
  R61=(z[5] - z[0]);
  R74=(z[6] - z[3]);

  R72=(z[6] - z[1]);
  R75=(z[6] - z[4]);
  R81=(z[7] - z[0]);

  t1=(R63 + R54);
  t2=(R61 + R74);
  t3=(R72 + R81);

  t4 =(R86 + R42);
  t5 =(R83 + R52);
  t6 =(R75 + R31);

  //  Calculate X gradient from Y and Z data

  grad_x[0] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
  grad_x[1] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
  grad_x[2] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
  grad_x[3] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
  grad_x[4] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
  grad_x[5] = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
  grad_x[6] = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
  grad_x[7] = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

  //  calc Y difference vectors

  R42=(y[3] - y[1]);
  R52=(y[4] - y[1]);
  R54=(y[4] - y[3]);

  R63=(y[5] - y[2]);
  R83=(y[7] - y[2]);
  R86=(y[7] - y[5]);

  R31=(y[2] - y[0]);
  R61=(y[5] - y[0]);
  R74=(y[6] - y[3]);

  R72=(y[6] - y[1]);
  R75=(y[6] - y[4]);
  R81=(y[7] - y[0]);

  t1=(R63 + R54);
  t2=(R61 + R74);
  t3=(R72 + R81);

  t4 =(R86 + R42);
  t5 =(R83 + R52);
  t6 =(R75 + R31);

  //  Calculate Z gradient from X and Y data

  grad_z[0] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
  grad_z[1] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
  grad_z[2] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
  grad_z[3] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
  grad_z[4] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
  grad_z[5] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
  grad_z[6] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
  grad_z[7] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);
}

//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct initialize_element
{
  typedef DeviceType     execution_space ;

  typedef Explicit::Fields< Scalar , execution_space > Fields ;

  typename Fields::elem_node_ids_type      elem_node_connectivity ;
  typename Fields::node_coords_type        model_coords ;
  typename Fields::elem_sym_tensor_type    stretch ;
  typename Fields::elem_tensor_state_type  rotation ;
  typename Fields::array_type              elem_mass ;

  const Scalar density ;

  initialize_element( const Fields & mesh_fields )
    : elem_node_connectivity( mesh_fields.elem_node_connectivity )
    , model_coords(           mesh_fields.model_coords )
    , stretch(                mesh_fields.stretch )
    , rotation(               mesh_fields.rotation )
    , elem_mass(              mesh_fields.elem_mass )
    , density(                mesh_fields.density )
    {}

  KOKKOS_INLINE_FUNCTION
  void operator()( int ielem )const
  {
    const int K_XX = 0 ;
    const int K_YY = 1 ;
    const int K_ZZ = 2 ;
    const Scalar ONE12TH = 1.0 / 12.0 ;

    Scalar x[ Fields::ElemNodeCount ];
    Scalar y[ Fields::ElemNodeCount ];
    Scalar z[ Fields::ElemNodeCount ];
    Scalar grad_x[ Fields::ElemNodeCount ];
    Scalar grad_y[ Fields::ElemNodeCount ];
    Scalar grad_z[ Fields::ElemNodeCount ];

    for ( int i = 0 ; i < Fields::ElemNodeCount ; ++i ) {
      const int n = elem_node_connectivity( ielem , i );

      x[i]  = model_coords( n , 0 );
      y[i]  = model_coords( n , 1 );
      z[i]  = model_coords( n , 2 );
    }

    comp_grad<Scalar,execution_space>( x, y, z, grad_x, grad_y, grad_z);

    stretch(ielem,K_XX) = 1 ;
    stretch(ielem,K_YY) = 1 ;
    stretch(ielem,K_ZZ) = 1 ;

    rotation(ielem,K_XX,0) = 1 ;
    rotation(ielem,K_YY,0) = 1 ;
    rotation(ielem,K_ZZ,0) = 1 ;

    rotation(ielem,K_XX,1) = 1 ;
    rotation(ielem,K_YY,1) = 1 ;
    rotation(ielem,K_ZZ,1) = 1 ;

    elem_mass(ielem) = ONE12TH * density *
                                 dot8<Scalar,execution_space>( x , grad_x );
  }

  static void apply( const Fields & mesh_fields )
  {
    initialize_element op( mesh_fields );
    Kokkos::parallel_for( mesh_fields.num_elements , op );
  }
};


template<typename Scalar , class DeviceType >
struct initialize_node
{
  typedef DeviceType     execution_space ;

  typedef Explicit::Fields< Scalar , execution_space > Fields ;

  typename Fields::node_elem_ids_type      node_elem_connectivity ;
  typename Fields::array_type              nodal_mass ;
  typename Fields::array_type              elem_mass ;

  static const int ElemNodeCount = Fields::ElemNodeCount ;

  initialize_node( const Fields & mesh_fields )
    : node_elem_connectivity( mesh_fields.node_elem_connectivity )
    , nodal_mass(             mesh_fields.nodal_mass )
    , elem_mass(              mesh_fields.elem_mass )
    {}


  KOKKOS_INLINE_FUNCTION
  void operator()( int inode )const
  {
    const int begin = node_elem_connectivity.row_map[inode];
    const int end   = node_elem_connectivity.row_map[inode+1];

    Scalar node_mass = 0;

    for(int i = begin; i != end; ++i) {
      const int elem_id = node_elem_connectivity.entries( i , 0 );
      node_mass += elem_mass(elem_id);
    }

    nodal_mass(inode) = node_mass / ElemNodeCount ;
  }

  static void apply( const Fields & mesh_fields )
  {
    initialize_node op( mesh_fields );
    Kokkos::parallel_for( mesh_fields.num_nodes_owned , op );
  }
};

//----------------------------------------------------------------------------


template<typename Scalar, class DeviceType >
struct grad
{
  typedef DeviceType execution_space ;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  static const int ElemNodeCount = Fields::ElemNodeCount ;

  static const int K_F_XX = Fields::K_F_XX ;
  static const int K_F_YY = Fields::K_F_YY ;
  static const int K_F_ZZ = Fields::K_F_ZZ ;
  static const int K_F_XY = Fields::K_F_XY ;
  static const int K_F_YZ = Fields::K_F_YZ ;
  static const int K_F_ZX = Fields::K_F_ZX ;
  static const int K_F_YX = Fields::K_F_YX ;
  static const int K_F_ZY = Fields::K_F_ZY ;
  static const int K_F_XZ = Fields::K_F_XZ ;

  // Global arrays used by this functor.

  const typename Fields::elem_node_ids_type     elem_node_connectivity ;
  const typename Fields::node_coords_type       model_coords ;
  const typename Fields::geom_state_array_type  displacement ;
  const typename Fields::geom_state_array_type  velocity ;
  const typename Fields::elem_tensor_type       vel_grad ;
  const typename Fields::scalar_type            dt ;

  const int  current_state;
  const int  previous_state;

  // Constructor on the Host to populate this device functor.
  // All array view copies are shallow.
  grad( const Fields &  fields,
        const int arg_current_state,
        const int arg_previous_state)
    : elem_node_connectivity( fields.elem_node_connectivity)
    , model_coords( fields.model_coords)
    , displacement( fields.displacement)
    , velocity( fields.velocity)
    , vel_grad( fields.vel_grad)
    , dt(  fields.dt)
    , current_state(arg_current_state)
    , previous_state(arg_previous_state)
    { }

  //--------------------------------------------------------------------------

    //   Calculate Velocity Gradients
    KOKKOS_INLINE_FUNCTION
    void v_grad(  int ielem,
      Scalar * vx,       Scalar * vy,       Scalar * vz,
      Scalar * grad_x,     Scalar * grad_y,     Scalar * grad_z,
      Scalar inv_vol) const
    {
      const int K_F_XX = Fields::K_F_XX ;
      const int K_F_YY = Fields::K_F_YY ;
      const int K_F_ZZ = Fields::K_F_ZZ ;
      const int K_F_XY = Fields::K_F_XY ;
      const int K_F_YZ = Fields::K_F_YZ ;
      const int K_F_ZX = Fields::K_F_ZX ;
      const int K_F_YX = Fields::K_F_YX ;
      const int K_F_ZY = Fields::K_F_ZY ;
      const int K_F_XZ = Fields::K_F_XZ ;

      vel_grad(ielem, K_F_XX) = inv_vol * dot8<Scalar,execution_space>( vx , grad_x );
      vel_grad(ielem, K_F_YX) = inv_vol * dot8<Scalar,execution_space>( vy , grad_x );
      vel_grad(ielem, K_F_ZX) = inv_vol * dot8<Scalar,execution_space>( vz , grad_x );

      vel_grad(ielem, K_F_XY) = inv_vol * dot8<Scalar,execution_space>( vx , grad_y );
      vel_grad(ielem, K_F_YY) = inv_vol * dot8<Scalar,execution_space>( vy , grad_y );
      vel_grad(ielem, K_F_ZY) = inv_vol * dot8<Scalar,execution_space>( vz , grad_y );

      vel_grad(ielem, K_F_XZ) = inv_vol * dot8<Scalar,execution_space>( vx , grad_z );
      vel_grad(ielem, K_F_YZ) = inv_vol * dot8<Scalar,execution_space>( vy , grad_z );
      vel_grad(ielem, K_F_ZZ) = inv_vol * dot8<Scalar,execution_space>( vz , grad_z );
    }

  //--------------------------------------------------------------------------
  // Functor operator() which calls the three member functions.


  KOKKOS_INLINE_FUNCTION
  void operator()( int ielem )const
  {
    const int X = 0 ;
    const int Y = 1 ;
    const int Z = 2 ;
    const Scalar dt_scale = -0.5 * *dt;

    //  declare and reuse local data for frequently accessed data to
    //  reduce global memory reads and writes.

    Scalar      x[8],      y[8],      z[8];
    Scalar     vx[8],     vy[8],     vz[8];
    Scalar grad_x[8], grad_y[8], grad_z[8];

    // Read global velocity once and use many times
    // via local registers / L1 cache.
    //  store the velocity information in local memory before using,
    //  so it can be returned for other functions to use

    // Read global coordinates and velocity once and use many times
    // via local registers / L1 cache.
    // load X coordinate information and move by half time step

    for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
      const int n = elem_node_connectivity( ielem , i );

      vx[i] = velocity( n , X , current_state );
      vy[i] = velocity( n , Y , current_state );
      vz[i] = velocity( n , Z , current_state );

      x[i]  = model_coords( n , X ) +
              displacement( n , X , current_state ) +
              dt_scale * vx[i];

      y[i]  = model_coords( n , Y ) +
              displacement( n , Y , current_state ) +
              dt_scale * vy[i];

      z[i]  = model_coords( n , Z ) +
              displacement( n , Z , current_state ) +
              dt_scale * vz[i];
    }

    comp_grad<Scalar,execution_space>( x, y, z, grad_x, grad_y, grad_z);

    //  Calculate hexahedral volume from x model_coords and gradient information

    const Scalar inv_vol = 1.0 / dot8<Scalar,execution_space>( x , grad_x );

    v_grad(ielem, vx, vy, vz, grad_x, grad_y, grad_z, inv_vol);
  }

  static void apply( const Fields & fields ,
                     const int arg_current_state ,
                     const int arg_previous_state )
  {
    grad op( fields, arg_current_state , arg_previous_state );
    Kokkos::parallel_for( fields.num_elements , op );
  }
};

//----------------------------------------------------------------------------

template<typename Scalar, class DeviceType >
struct decomp_rotate
{
  typedef DeviceType execution_space ;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  static const int ElemNodeCount = Fields::ElemNodeCount ;

  static const int K_F_XX = Fields::K_F_XX ;
  static const int K_F_YY = Fields::K_F_YY ;
  static const int K_F_ZZ = Fields::K_F_ZZ ;
  static const int K_F_XY = Fields::K_F_XY ;
  static const int K_F_YZ = Fields::K_F_YZ ;
  static const int K_F_ZX = Fields::K_F_ZX ;
  static const int K_F_YX = Fields::K_F_YX ;
  static const int K_F_ZY = Fields::K_F_ZY ;
  static const int K_F_XZ = Fields::K_F_XZ ;

  static const int K_S_XX = Fields::K_S_XX ;
  static const int K_S_YY = Fields::K_S_YY ;
  static const int K_S_ZZ = Fields::K_S_ZZ ;
  static const int K_S_XY = Fields::K_S_XY ;
  static const int K_S_YZ = Fields::K_S_YZ ;
  static const int K_S_ZX = Fields::K_S_ZX ;
  static const int K_S_YX = Fields::K_S_YX ;
  static const int K_S_ZY = Fields::K_S_ZY ;
  static const int K_S_XZ = Fields::K_S_XZ ;

  static const int K_V_XY = Fields::K_V_XY ;
  static const int K_V_YZ = Fields::K_V_YZ ;
  static const int K_V_ZX = Fields::K_V_ZX ;

  // Global arrays used by this functor.

  const typename Fields::elem_tensor_state_type     rotation ;
  const typename Fields::elem_tensor_type           vel_grad ;
  const typename Fields::elem_sym_tensor_type       stretch ;
  const typename Fields::elem_sym_tensor_type       rot_stretch ;
  const typename Fields::scalar_type                dt_value ;

  const int  current_state;
  const int  previous_state;

  decomp_rotate( const Fields & mesh_fields ,
                 const int arg_current_state,
                 const int arg_previous_state)
    : rotation(    mesh_fields.rotation )
    , vel_grad(    mesh_fields.vel_grad )
    , stretch(     mesh_fields.stretch )
    , rot_stretch( mesh_fields.rot_stretch )
    , dt_value(    mesh_fields.dt)
    , current_state( arg_current_state)
    , previous_state(arg_previous_state)
    {}

  static void apply( const Fields & mesh_fields ,
                     const int arg_current_state ,
                     const int arg_previous_state )
  {
    decomp_rotate op( mesh_fields , arg_current_state , arg_previous_state );
    Kokkos::parallel_for( mesh_fields.num_elements , op );
  }


  KOKKOS_INLINE_FUNCTION
  void additive_decomp(int ielem, Scalar * v_gr, Scalar * str_ten) const
  {
    //  In addition to calculating stretching_tensor,
    //  use this as an opportunity to load global
    //  variables into a local space

    for ( int i = 0 ; i < 9 ; ++i ) {
      v_gr[i] = vel_grad( ielem , i );
    }

    //
    //  Symmetric part
    //
    str_ten[K_S_XX] = v_gr[K_F_XX];
    str_ten[K_S_YY] = v_gr[K_F_YY];
    str_ten[K_S_ZZ] = v_gr[K_F_ZZ];
    str_ten[K_S_XY] = 0.5*(v_gr[K_F_XY] + v_gr[K_F_YX]);
    str_ten[K_S_YZ] = 0.5*(v_gr[K_F_YZ] + v_gr[K_F_ZY]);
    str_ten[K_S_ZX] = 0.5*(v_gr[K_F_ZX] + v_gr[K_F_XZ]);
  }

  KOKKOS_INLINE_FUNCTION
  void polar_decomp(int ielem, Scalar * v_gr, Scalar * str_ten, Scalar * str, Scalar * vort, Scalar * rot_old, Scalar * rot_new)const
  {
    const Scalar dt = *dt_value;
    const Scalar dt_half = 0.5 * dt;

    //  Skew Symmetric part
    vort[K_V_XY] = 0.5*(v_gr[K_F_XY] - v_gr[K_F_YX]);
    vort[K_V_YZ] = 0.5*(v_gr[K_F_YZ] - v_gr[K_F_ZY]);
    vort[K_V_ZX] = 0.5*(v_gr[K_F_ZX] - v_gr[K_F_XZ]);

    //   calculate the rates of rotation via gauss elimination.
    for ( int i = 0 ; i < 6 ; ++i ) {
      str[i] = stretch(ielem, i);
    }

    Scalar z1 = str_ten[K_S_XY] * str[K_S_ZX] -
                str_ten[K_S_ZX] * str[K_S_XY] +
                str_ten[K_S_YY] * str[K_S_YZ] -
                str_ten[K_S_YZ] * str[K_S_YY] +
                str_ten[K_S_YZ] * str[K_S_ZZ] -
                str_ten[K_S_ZZ] * str[K_S_YZ];

    Scalar z2 = str_ten[K_S_ZX] * str[K_S_XX] -
                str_ten[K_S_XX] * str[K_S_ZX] +
                str_ten[K_S_YZ] * str[K_S_XY] -
                str_ten[K_S_XY] * str[K_S_YZ] +
                str_ten[K_S_ZZ] * str[K_S_ZX] -
                str_ten[K_S_ZX] * str[K_S_ZZ];

    Scalar z3 = str_ten[K_S_XX] * str[K_S_XY] -
                str_ten[K_S_XY] * str[K_S_XX] +
                str_ten[K_S_XY] * str[K_S_YY] -
                str_ten[K_S_YY] * str[K_S_XY] +
                str_ten[K_S_ZX] * str[K_S_YZ] -
                str_ten[K_S_YZ] * str[K_S_ZX];

  //   forward elimination
    const Scalar a1inv = 1.0 / (str[K_S_YY] + str[K_S_ZZ]);

    const Scalar a4BYa1 = -1 * str[K_S_XY] * a1inv;

    const Scalar a2inv = 1.0 / (str[K_S_ZZ] + str[K_S_XX] + str[K_S_XY] * a4BYa1);

    const Scalar a5 =  -str[K_S_YZ] + str[K_S_ZX] * a4BYa1;

    z2 -= z1 * a4BYa1;
    Scalar a6BYa1 = -1 * str[K_S_ZX] * a1inv;
    const Scalar a5BYa2 = a5 * a2inv;
    z3 -= z1 * a6BYa1 - z2 * a5BYa2;

  //   backward substitution -
    z3 /= (str[K_S_XX] + str[K_S_YY] + str[K_S_ZX] * a6BYa1 + a5 * a5BYa2);
    z2 = (z2 - a5 * z3) * a2inv;
    z1 = (z1*a1inv - a6BYa1 * z3 -a4BYa1 * z2);

  //   calculate rotation rates - recall that spin_rate is an asymmetric tensor,
  //   so compute spin rate vector as dual of spin rate tensor,
  //   i.e   w_i = e_ijk * spin_rate_jk
    z1 += vort[K_V_YZ];
    z2 += vort[K_V_ZX];
    z3 += vort[K_V_XY];

  //   update rotation tensor:
  //  1) premultiply old rotation tensor to get right-hand side.

    for ( int i = 0 ; i < 9 ; ++i ) {
      rot_old[i] = rotation(ielem, i, previous_state);
    }

    Scalar r_XX = rot_old[K_F_XX] + dt_half*( z3 * rot_old[K_F_YX] - z2 * rot_old[K_F_ZX] );
    Scalar r_YX = rot_old[K_F_YX] + dt_half*( z1 * rot_old[K_F_ZX] - z3 * rot_old[K_F_XX] );
    Scalar r_ZX = rot_old[K_F_ZX] + dt_half*( z2 * rot_old[K_F_XX] - z1 * rot_old[K_F_YX] );
    Scalar r_XY = rot_old[K_F_XY] + dt_half*( z3 * rot_old[K_F_YY] - z2 * rot_old[K_F_ZY] );
    Scalar r_YY = rot_old[K_F_YY] + dt_half*( z1 * rot_old[K_F_ZY] - z3 * rot_old[K_F_XY] );
    Scalar r_ZY = rot_old[K_F_ZY] + dt_half*( z2 * rot_old[K_F_XY] - z1 * rot_old[K_F_YY] );
    Scalar r_XZ = rot_old[K_F_XZ] + dt_half*( z3 * rot_old[K_F_YZ] - z2 * rot_old[K_F_ZZ] );
    Scalar r_YZ = rot_old[K_F_YZ] + dt_half*( z1 * rot_old[K_F_ZZ] - z3 * rot_old[K_F_XZ] );
    Scalar r_ZZ = rot_old[K_F_ZZ] + dt_half*( z2 * rot_old[K_F_XZ] - z1 * rot_old[K_F_YZ] );


  //  2) solve for new rotation tensor via gauss elimination.
  //   forward elimination -
    Scalar a12 = - dt_half * z3;
    Scalar a13 =   dt_half * z2;
    Scalar b32 = - dt_half * z1;
    Scalar a22inv = 1.0 / (1.0 + a12 * a12);

    Scalar a13a12 = a13*a12;
    Scalar a23 = b32 + a13a12;
    r_YX += r_XX * a12;
    r_YY += r_XY * a12;
    r_YZ += r_XZ * a12;


    b32 = (b32 - a13a12) * a22inv;
    r_ZX += r_XX * a13 + r_YX * b32;
    r_ZY += r_XY * a13 + r_YY * b32;
    r_ZZ += r_XZ * a13 + r_YZ * b32;


  //   backward substitution -
    const Scalar a33inv = 1.0 / (1.0 + a13 * a13 + a23 * b32);

    rot_new[K_F_ZX] = r_ZX * a33inv;
    rot_new[K_F_ZY] = r_ZY * a33inv;
    rot_new[K_F_ZZ] = r_ZZ * a33inv;
    rot_new[K_F_YX] = ( r_YX - rot_new[K_F_ZX] * a23 ) * a22inv;
    rot_new[K_F_YY] = ( r_YY - rot_new[K_F_ZY] * a23 ) * a22inv;
    rot_new[K_F_YZ] = ( r_YZ - rot_new[K_F_ZZ] * a23 ) * a22inv;
    rot_new[K_F_XX] = r_XX - rot_new[K_F_ZX] * a13 - rot_new[K_F_YX] * a12;
    rot_new[K_F_XY] = r_XY - rot_new[K_F_ZY] * a13 - rot_new[K_F_YY] * a12;
    rot_new[K_F_XZ] = r_XZ - rot_new[K_F_ZZ] * a13 - rot_new[K_F_YZ] * a12;

    for ( int i = 0 ; i < 9 ; ++i ) {
      rotation(ielem, i, current_state) = rot_new[i] ;
    }

  //   update stretch tensor in the new configuration -
    const Scalar a1 = str_ten[K_S_XY] + vort[K_V_XY];
    const Scalar a2 = str_ten[K_S_YZ] + vort[K_V_YZ];
    const Scalar a3 = str_ten[K_S_ZX] + vort[K_V_ZX];
    const Scalar b1 = str_ten[K_S_ZX] - vort[K_V_ZX];
    const Scalar b2 = str_ten[K_S_XY] - vort[K_V_XY];
    const Scalar b3 = str_ten[K_S_YZ] - vort[K_V_YZ];

    const Scalar s_XX = str[K_S_XX];
    const Scalar s_YY = str[K_S_YY];
    const Scalar s_ZZ = str[K_S_ZZ];
    const Scalar s_XY = str[K_S_XY];
    const Scalar s_YZ = str[K_S_YZ];
    const Scalar s_ZX = str[K_S_ZX];

    str[K_S_XX] += dt * (str_ten[K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
    str[K_S_YY] += dt * (str_ten[K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
    str[K_S_ZZ] += dt * (str_ten[K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
    str[K_S_XY] += dt * (str_ten[K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
    str[K_S_YZ] += dt * (str_ten[K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
    str[K_S_ZX] += dt * (str_ten[K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);

  }


  KOKKOS_INLINE_FUNCTION
  void rotate_tensor(int ielem, Scalar * str_ten, Scalar * str, Scalar * rot_new)const {

    Scalar t[9];
    Scalar rot_str[6]; // Rotated stretch

    t[0] = str_ten[K_S_XX]*rot_new[K_F_XX] +
           str_ten[K_S_XY]*rot_new[K_F_YX] +
           str_ten[K_S_XZ]*rot_new[K_F_ZX];

    t[1] = str_ten[K_S_YX]*rot_new[K_F_XX] +
           str_ten[K_S_YY]*rot_new[K_F_YX] +
           str_ten[K_S_YZ]*rot_new[K_F_ZX];

    t[2] = str_ten[K_S_ZX]*rot_new[K_F_XX] +
           str_ten[K_S_ZY]*rot_new[K_F_YX] +
           str_ten[K_S_ZZ]*rot_new[K_F_ZX];

    t[3] = str_ten[K_S_XX]*rot_new[K_F_XY] +
           str_ten[K_S_XY]*rot_new[K_F_YY] +
           str_ten[K_S_XZ]*rot_new[K_F_ZY];

    t[4] = str_ten[K_S_YX]*rot_new[K_F_XY] +
           str_ten[K_S_YY]*rot_new[K_F_YY] +
           str_ten[K_S_YZ]*rot_new[K_F_ZY];

    t[5] = str_ten[K_S_ZX]*rot_new[K_F_XY] +
           str_ten[K_S_ZY]*rot_new[K_F_YY] +
           str_ten[K_S_ZZ]*rot_new[K_F_ZY];

    t[6] = str_ten[K_S_XX]*rot_new[K_F_XZ] +
           str_ten[K_S_XY]*rot_new[K_F_YZ] +
           str_ten[K_S_XZ]*rot_new[K_F_ZZ];

    t[7] = str_ten[K_S_YX]*rot_new[K_F_XZ] +
           str_ten[K_S_YY]*rot_new[K_F_YZ] +
           str_ten[K_S_YZ]*rot_new[K_F_ZZ];

    t[8] = str_ten[K_S_ZX]*rot_new[K_F_XZ] +
           str_ten[K_S_ZY]*rot_new[K_F_YZ] +
           str_ten[K_S_ZZ]*rot_new[K_F_ZZ];


    rot_str[ K_S_XX ] = rot_new[K_F_XX] * t[0] +
                        rot_new[K_F_YX] * t[1] +
                        rot_new[K_F_ZX] * t[2];
    rot_str[ K_S_YY ] = rot_new[K_F_XY] * t[3] +
                        rot_new[K_F_YY] * t[4] +
                        rot_new[K_F_ZY] * t[5];
    rot_str[ K_S_ZZ ] = rot_new[K_F_XZ] * t[6] +
                        rot_new[K_F_YZ] * t[7] +
                        rot_new[K_F_ZZ] * t[8];

    rot_str[ K_S_XY ] = rot_new[K_F_XX] * t[3] +
                        rot_new[K_F_YX] * t[4] +
                        rot_new[K_F_ZX] * t[5];
    rot_str[ K_S_YZ ] = rot_new[K_F_XY] * t[6] +
                        rot_new[K_F_YY] * t[7] +
                        rot_new[K_F_ZY] * t[8];
    rot_str[ K_S_ZX ] = rot_new[K_F_XZ] * t[0] +
                        rot_new[K_F_YZ] * t[1] +
                        rot_new[K_F_ZZ] * t[2];

    for ( int i = 0 ; i < 6 ; ++i ) {
      rot_stretch(ielem, i) = rot_str[i] ;
    }

    for ( int i = 0 ; i < 6 ; ++i ) {
      stretch(ielem, i) = str[i] ;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( int ielem )const {

    //   Local scratch space to avoid multiple
    //   accesses to global memory.
    Scalar str_ten[6]; // Stretching tensor
    Scalar str[6];     // Stretch
    Scalar rot_old[9]; // Rotation old
    Scalar rot_new[9]; // Rotation new
    Scalar vort[3];    // Vorticity
    Scalar v_gr[9];    // Velocity gradient

    additive_decomp(ielem, v_gr, str_ten);

    polar_decomp(ielem, v_gr, str_ten, str, vort, rot_old, rot_new);

    rotate_tensor(ielem, str_ten, str, rot_new);
  }
};

//----------------------------------------------------------------------------

template<typename Scalar, class DeviceType >
struct internal_force
{
  typedef DeviceType execution_space ;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  static const int ElemNodeCount = Fields::ElemNodeCount ;

  static const int K_F_XX = Fields::K_F_XX ;
  static const int K_F_YY = Fields::K_F_YY ;
  static const int K_F_ZZ = Fields::K_F_ZZ ;
  static const int K_F_XY = Fields::K_F_XY ;
  static const int K_F_YZ = Fields::K_F_YZ ;
  static const int K_F_ZX = Fields::K_F_ZX ;
  static const int K_F_YX = Fields::K_F_YX ;
  static const int K_F_ZY = Fields::K_F_ZY ;
  static const int K_F_XZ = Fields::K_F_XZ ;

  static const int K_S_XX = Fields::K_S_XX ;
  static const int K_S_YY = Fields::K_S_YY ;
  static const int K_S_ZZ = Fields::K_S_ZZ ;
  static const int K_S_XY = Fields::K_S_XY ;
  static const int K_S_YZ = Fields::K_S_YZ ;
  static const int K_S_ZX = Fields::K_S_ZX ;
  static const int K_S_YX = Fields::K_S_YX ;
  static const int K_S_ZY = Fields::K_S_ZY ;
  static const int K_S_XZ = Fields::K_S_XZ ;

  //--------------------------------------------------------------------------
  // Reduction:

  typedef Scalar value_type;

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) {
    update = 1.0e32;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update,
                    const volatile value_type & source )
  {
    update = update < source ? update : source;
  }

  // Final serial processing of reduction value:
  KOKKOS_INLINE_FUNCTION
  void final( value_type & result ) const
  {
    *prev_dt = *dt ;
    *dt = result ;
  };

  //--------------------------------------------------------------------------

  // Global arrays used by this functor.

  const typename Fields::elem_node_ids_type      elem_node_connectivity ;
  const typename Fields::node_coords_type        model_coords ;
  const typename Fields::scalar_type             dt ;
  const typename Fields::scalar_type             prev_dt ;
  const typename Fields::geom_state_array_type   displacement ;
  const typename Fields::geom_state_array_type   velocity ;
  const typename Fields::array_type              elem_mass ;
  const typename Fields::array_type              internal_energy ;
  const typename Fields::elem_sym_tensor_type    stress_new ;
  const typename Fields::elem_node_geom_type     element_force ;
  const typename Fields::elem_tensor_state_type  rotation ;
  const typename Fields::elem_sym_tensor_type    rot_stretch ;

  const Scalar     two_mu;
  const Scalar     bulk_modulus;
  const Scalar     lin_bulk_visc;
  const Scalar     quad_bulk_visc;
  const Scalar     user_dt;
  const int        current_state;

  internal_force( const Fields & mesh_fields,
                  const Scalar arg_user_dt,
                  const int arg_current_state )
    : elem_node_connectivity( mesh_fields.elem_node_connectivity )
    , model_coords(           mesh_fields.model_coords )
    , dt(                     mesh_fields.dt )
    , prev_dt(                mesh_fields.prev_dt )
    , displacement(           mesh_fields.displacement )
    , velocity(               mesh_fields.velocity )
    , elem_mass(              mesh_fields.elem_mass )
    , internal_energy(        mesh_fields.internal_energy )
    , stress_new(             mesh_fields.stress_new )
    , element_force(          mesh_fields.element_force )
    , rotation(               mesh_fields.rotation )
    , rot_stretch(            mesh_fields.rot_stretch )
    , two_mu(                 mesh_fields.two_mu )
    , bulk_modulus(           mesh_fields.bulk_modulus )
    , lin_bulk_visc(          mesh_fields.lin_bulk_visc )
    , quad_bulk_visc(         mesh_fields.quad_bulk_visc )
    , user_dt(       arg_user_dt )
    , current_state( arg_current_state )
  {}

  static void apply( const Fields & mesh_fields ,
                     const Scalar arg_user_dt,
                     const int arg_current_state )
  {
    internal_force  op_force( mesh_fields , arg_user_dt , arg_current_state );

    Kokkos::parallel_reduce( mesh_fields.num_elements, op_force );
  }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void rotate_tensor_backward(int ielem ,
    const Scalar * const s_n ,
    Scalar * const rot_stress )const
  {
    const int rot_state = current_state ; // 1 ;

    //   t : temporary variables
    //   s_n : stress_new in local memory space
    //   r_n : rotation_new in local memory space
    Scalar t[9], r_n[9];

    r_n[0] = rotation(ielem, 0, rot_state );
    r_n[1] = rotation(ielem, 1, rot_state );
    r_n[2] = rotation(ielem, 2, rot_state );
    r_n[3] = rotation(ielem, 3, rot_state );
    r_n[4] = rotation(ielem, 4, rot_state );
    r_n[5] = rotation(ielem, 5, rot_state );
    r_n[6] = rotation(ielem, 6, rot_state );
    r_n[7] = rotation(ielem, 7, rot_state );
    r_n[8] = rotation(ielem, 8, rot_state );

    t[0] = s_n[K_S_XX]*r_n[K_F_XX]+ s_n[K_S_XY]*r_n[K_F_XY]+ s_n[K_S_XZ]*r_n[K_F_XZ];
    t[1] = s_n[K_S_YX]*r_n[K_F_XX]+ s_n[K_S_YY]*r_n[K_F_XY]+ s_n[K_S_YZ]*r_n[K_F_XZ];
    t[2] = s_n[K_S_ZX]*r_n[K_F_XX]+ s_n[K_S_ZY]*r_n[K_F_XY]+ s_n[K_S_ZZ]*r_n[K_F_XZ];
    t[3] = s_n[K_S_XX]*r_n[K_F_YX]+ s_n[K_S_XY]*r_n[K_F_YY]+ s_n[K_S_XZ]*r_n[K_F_YZ];
    t[4] = s_n[K_S_YX]*r_n[K_F_YX]+ s_n[K_S_YY]*r_n[K_F_YY]+ s_n[K_S_YZ]*r_n[K_F_YZ];
    t[5] = s_n[K_S_ZX]*r_n[K_F_YX]+ s_n[K_S_ZY]*r_n[K_F_YY]+ s_n[K_S_ZZ]*r_n[K_F_YZ];
    t[6] = s_n[K_S_XX]*r_n[K_F_ZX]+ s_n[K_S_XY]*r_n[K_F_ZY]+ s_n[K_S_XZ]*r_n[K_F_ZZ];
    t[7] = s_n[K_S_YX]*r_n[K_F_ZX]+ s_n[K_S_YY]*r_n[K_F_ZY]+ s_n[K_S_YZ]*r_n[K_F_ZZ];
    t[8] = s_n[K_S_ZX]*r_n[K_F_ZX]+ s_n[K_S_ZY]*r_n[K_F_ZY]+ s_n[K_S_ZZ]*r_n[K_F_ZZ];

    rot_stress[ K_S_XX ] = r_n[K_F_XX]*t[0] + r_n[K_F_XY]*t[1] + r_n[K_F_XZ]*t[2];
    rot_stress[ K_S_YY ] = r_n[K_F_YX]*t[3] + r_n[K_F_YY]*t[4] + r_n[K_F_YZ]*t[5];
    rot_stress[ K_S_ZZ ] = r_n[K_F_ZX]*t[6] + r_n[K_F_ZY]*t[7] + r_n[K_F_ZZ]*t[8];

    rot_stress[ K_S_XY ] = r_n[K_F_XX]*t[3] + r_n[K_F_XY]*t[4] + r_n[K_F_XZ]*t[5];
    rot_stress[ K_S_YZ ] = r_n[K_F_YX]*t[6] + r_n[K_F_YY]*t[7] + r_n[K_F_YZ]*t[8];
    rot_stress[ K_S_ZX ] = r_n[K_F_ZX]*t[0] + r_n[K_F_ZY]*t[1] + r_n[K_F_ZZ]*t[2];
  }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void comp_force(int ielem,
     const Scalar * const vx ,
     const Scalar * const vy ,
     const Scalar * const vz ,
     const Scalar * const grad_x ,
     const Scalar * const grad_y ,
     const Scalar * const grad_z ,
     Scalar * total_stress12th ) const
  {
    Scalar internal_energy_inc = 0 ;

    for(int inode = 0; inode < 8; ++inode) {

      const Scalar fx =
        total_stress12th[K_S_XX] * grad_x[inode] +
        total_stress12th[K_S_XY] * grad_y[inode] +
        total_stress12th[K_S_XZ] * grad_z[inode] ;

      element_force(ielem, 0, inode) = fx ;

      const Scalar fy =
        total_stress12th[K_S_YX] * grad_x[inode] +
        total_stress12th[K_S_YY] * grad_y[inode] +
        total_stress12th[K_S_YZ] * grad_z[inode] ;

      element_force(ielem, 1, inode) = fy ;

      const Scalar fz =
        total_stress12th[K_S_ZX] * grad_x[inode] +
        total_stress12th[K_S_ZY] * grad_y[inode] +
        total_stress12th[K_S_ZZ] * grad_z[inode] ;

      element_force(ielem, 2, inode) = fz ;

      internal_energy_inc +=
        fx * vx[inode] +
        fy * vy[inode] +
        fz * vz[inode] ;
    }

    internal_energy(ielem) = internal_energy_inc ;
  }

  //----------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void get_stress(int ielem , Scalar * const s_n ) const
    {
      const int kxx = 0;
      const int kyy = 1;
      const int kzz = 2;
      const int kxy = 3;
      const int kyz = 4;
      const int kzx = 5;

      const Scalar e = (rot_stretch(ielem,kxx)+rot_stretch(ielem,kyy)+rot_stretch(ielem,kzz))/3.0;

      s_n[kxx] = stress_new(ielem,kxx) += *dt * (two_mu * (rot_stretch(ielem,kxx)-e)+3*bulk_modulus*e);
      s_n[kyy] = stress_new(ielem,kyy) += *dt * (two_mu * (rot_stretch(ielem,kyy)-e)+3*bulk_modulus*e);
      s_n[kzz] = stress_new(ielem,kzz) += *dt * (two_mu * (rot_stretch(ielem,kzz)-e)+3*bulk_modulus*e);

      s_n[kxy] = stress_new(ielem,kxy) += *dt * two_mu * rot_stretch(ielem,kxy);
      s_n[kyz] = stress_new(ielem,kyz) += *dt * two_mu * rot_stretch(ielem,kyz);
      s_n[kzx] = stress_new(ielem,kzx) += *dt * two_mu * rot_stretch(ielem,kzx);
    }

  //----------------------------------------------------------------------------


  KOKKOS_INLINE_FUNCTION
  void operator()( int ielem, value_type & update )const
  {
    const Scalar ONE12TH = 1.0 / 12.0 ;

    Scalar x[8], y[8], z[8] ;
    Scalar vx[8], vy[8], vz[8];
    Scalar grad_x[8], grad_y[8], grad_z[8];

    // Position and velocity:

    for ( int i = 0 ; i < ElemNodeCount ; ++i ) {
      const int n = elem_node_connectivity(ielem,i);

      x[i] = model_coords(n, 0) + displacement(n, 0, current_state) ;
      y[i] = model_coords(n, 1) + displacement(n, 1, current_state) ;
      z[i] = model_coords(n, 2) + displacement(n, 2, current_state) ;

      vx[i] = velocity(n, 0, current_state);
      vy[i] = velocity(n, 1, current_state);
      vz[i] = velocity(n, 2, current_state);
    }

    // Gradient:

    comp_grad<Scalar,execution_space>( x , y , z , grad_x , grad_y , grad_z );


    const Scalar mid_vol = dot8<Scalar,execution_space>( x , grad_x );

    const Scalar shr = two_mu ;
    const Scalar dil = bulk_modulus + ((2.0*shr)/3.0);

    const Scalar aspect = 6.0 * mid_vol /
                          ( dot8<Scalar,execution_space>( grad_x , grad_x ) +
                            dot8<Scalar,execution_space>( grad_y , grad_y ) +
                            dot8<Scalar,execution_space>( grad_z , grad_z ) );

    const Scalar dtrial = std::sqrt(elem_mass(ielem) * aspect / dil);
    const Scalar traced = (rot_stretch(ielem, 0) + rot_stretch(ielem, 1) + rot_stretch(ielem, 2));

    const Scalar eps = traced < 0 ? (lin_bulk_visc - quad_bulk_visc * traced * dtrial) : lin_bulk_visc ;

    const Scalar bulkq = eps * dil * dtrial * traced;

    Scalar cur_time_step = dtrial * ( std::sqrt( 1.0 + eps * eps) - eps);

    // force fixed time step if input

    cur_time_step = user_dt > 0 ? user_dt : cur_time_step;

    update = update < cur_time_step ? update : cur_time_step;


    Scalar s_n[ 6 ];

    get_stress( ielem, s_n );

    Scalar total_stress12th[6];

    // Get rotated stress:

    rotate_tensor_backward(ielem, s_n , total_stress12th );

    total_stress12th[0] = ONE12TH*( total_stress12th[ 0 ] + bulkq );
    total_stress12th[1] = ONE12TH*( total_stress12th[ 1 ] + bulkq );
    total_stress12th[2] = ONE12TH*( total_stress12th[ 2 ] + bulkq );
    total_stress12th[3] = ONE12TH*( total_stress12th[ 3 ] );
    total_stress12th[4] = ONE12TH*( total_stress12th[ 4 ] );
    total_stress12th[5] = ONE12TH*( total_stress12th[ 5 ] );

    comp_force(ielem, vx, vy, vz,
                      grad_x, grad_y, grad_z, total_stress12th);
  }
};

//----------------------------------------------------------------------------

template<typename Scalar, class DeviceType >
struct nodal_step
{
  typedef DeviceType     execution_space ;
  typedef typename execution_space::size_type  size_type;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  const typename Fields::scalar_type            dt ;
  const typename Fields::scalar_type            prev_dt ;
  const typename Fields::node_elem_ids_type     node_elem_connectivity ;
  const typename Fields::node_coords_type       model_coords ;
  const typename Fields::array_type             nodal_mass ;
  const typename Fields::geom_state_array_type  displacement ;
  const typename Fields::geom_state_array_type  velocity ;
  const typename Fields::geom_array_type        acceleration ;
  const typename Fields::geom_array_type        internal_force ;
  const typename Fields::elem_node_geom_type    element_force ;

  const Scalar   x_bc;
  const int      current_state;
  const int      next_state;


  nodal_step( const Fields  & mesh_fields ,
              const Scalar    arg_x_bc,
              const int       arg_current_state,
              const int       arg_next_state)
   : dt(       mesh_fields.dt )
   , prev_dt(  mesh_fields.prev_dt )
   , node_elem_connectivity( mesh_fields.node_elem_connectivity )
   , model_coords(   mesh_fields.model_coords )
   , nodal_mass(     mesh_fields.nodal_mass )
   , displacement(   mesh_fields.displacement )
   , velocity(       mesh_fields.velocity )
   , acceleration(   mesh_fields.acceleration )
   , internal_force( mesh_fields.internal_force )
   , element_force(  mesh_fields.element_force )
   , x_bc(          arg_x_bc )
   , current_state( arg_current_state )
   , next_state(    arg_next_state )
   {
        //std::cout << "finish_step dt: " << dt << std::endl;
        //std::cout << "finish_step prev_dt: " << prev_dt << std::endl;
   }

  static void apply( const Fields  & mesh_fields ,
                     const Scalar    arg_x_bc ,
                     const int       arg_current_state ,
                     const int       arg_next_state )
  {
    nodal_step op( mesh_fields, arg_x_bc, arg_current_state, arg_next_state );

    // Only update the owned nodes:

    Kokkos::parallel_for( mesh_fields.num_nodes_owned , op );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int inode) const
    {
      // Getting count as per 'CSR-like' data structure
      const int begin = node_elem_connectivity.row_map[inode];
      const int end   = node_elem_connectivity.row_map[inode+1];

      double local_force[] = {0.0, 0.0, 0.0};

      // Gather-sum internal force from
      // each element that a node is attached to.

      for ( int i = begin; i < end ; ++i ){

        //  node_elem_offset is a cumulative structure, so
        //  node_elem_offset(inode) should be the index where
        //  a particular row's elem_IDs begin
        const int nelem = node_elem_connectivity.entries( i, 0);

        //  find the row in an element's stiffness matrix
        //  that corresponds to inode
        const int elem_node_index = node_elem_connectivity.entries( i, 1);

        local_force[0] += element_force(nelem, 0, elem_node_index);
        local_force[1] += element_force(nelem, 1, elem_node_index);
        local_force[2] += element_force(nelem, 2, elem_node_index);
      }

      internal_force(inode, 0) = local_force[0];
      internal_force(inode, 1) = local_force[1];
      internal_force(inode, 2) = local_force[2];

      // Acceleration:

      Scalar v_new[3];
      Scalar a_current[3];

      const Scalar tol = 1.0e-7;

      // If not on the boundary then: a = F / m
      if ( tol < fabs(model_coords(inode,0)-x_bc) ) {

        const Scalar m = nodal_mass( inode );

        acceleration(inode,0) = a_current[0] = -local_force[0] / m ;
        acceleration(inode,1) = a_current[1] = -local_force[1] / m ;
        acceleration(inode,2) = a_current[2] = -local_force[2] / m ;
      }
      else { //enforce fixed BC
        acceleration(inode,0) = a_current[0] = 0;
        acceleration(inode,1) = a_current[1] = 0;
        acceleration(inode,2) = a_current[2] = 0;
      }

      // Central difference time integration:

      const Scalar dt_disp = *dt ;
      const Scalar dt_vel = ( *dt + *prev_dt ) / 2.0 ;

      velocity(inode,0,next_state) = v_new[0] =
        velocity(inode,0,current_state) + dt_vel * a_current[0];

      velocity(inode,1,next_state) = v_new[1] =
        velocity(inode,1,current_state) + dt_vel * a_current[1];

      velocity(inode,2,next_state) = v_new[2] =
        velocity(inode,2,current_state) + dt_vel * a_current[2];

      displacement(inode,0,next_state) =
        displacement(inode,0,current_state) + dt_disp * v_new[0];

      displacement(inode,1,next_state) =
        displacement(inode,1,current_state) + dt_disp * v_new[1];

      displacement(inode,2,next_state) =
        displacement(inode,2,current_state) + dt_disp * v_new[2];
    }
};

//----------------------------------------------------------------------------

template< typename Scalar , class DeviceType >
struct pack_state
{
  typedef DeviceType     execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  typedef typename Fields::geom_state_array_type::value_type  value_type ;
  typedef Kokkos::View< value_type* , execution_space >     buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::geom_state_array_type  displacement ;
  const typename Fields::geom_state_array_type  velocity ;
  const buffer_type  output ;
  const size_type    inode_base ;
  const size_type    state_next ;

  pack_state( const buffer_type & arg_output ,
              const Fields      & mesh_fields ,
              const size_type     arg_begin ,
              const size_type     arg_state )
   : displacement( mesh_fields.displacement )
   , velocity(     mesh_fields.velocity )
   , output(       arg_output )
   , inode_base(   arg_begin )
   , state_next(   arg_state )
   {}

  static void apply( const buffer_type & arg_output ,
                     const size_type     arg_begin ,
                     const size_type     arg_count ,
                     const Fields      & mesh_fields ,
                     const size_type     arg_state )
  {
    pack_state op( arg_output , mesh_fields , arg_begin , arg_state );

    Kokkos::parallel_for( arg_count , op );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
    const size_type inode = inode_base + i ;

    size_type j = i * value_count ;

    output[j++] = displacement( inode , 0 , state_next );
    output[j++] = displacement( inode , 1 , state_next );
    output[j++] = displacement( inode , 2 , state_next );
    output[j++] = velocity( inode , 0 , state_next );
    output[j++] = velocity( inode , 1 , state_next );
    output[j++] = velocity( inode , 2 , state_next );
  }
};

template< typename Scalar , class DeviceType >
struct unpack_state
{
  typedef DeviceType     execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef Explicit::Fields< Scalar , execution_space >  Fields ;

  typedef typename Fields::geom_state_array_type::value_type  value_type ;
  typedef Kokkos::View< value_type* , execution_space >     buffer_type ;

  static const unsigned value_count = 6 ;

  const typename Fields::geom_state_array_type  displacement ;
  const typename Fields::geom_state_array_type  velocity ;
  const buffer_type  input ;
  const size_type    inode_base ;
  const size_type    state_next ;

  unpack_state( const buffer_type & arg_input ,
                const Fields      & mesh_fields ,
                const size_type     arg_begin ,
                const size_type     arg_state )
   : displacement( mesh_fields.displacement )
   , velocity(     mesh_fields.velocity )
   , input(        arg_input )
   , inode_base(   arg_begin )
   , state_next(   arg_state )
   {}

  static void apply( const Fields      & mesh_fields ,
                     const size_type     arg_state ,
                     const buffer_type & arg_input ,
                     const size_type     arg_begin ,
                     const size_type     arg_count )
  {
    unpack_state op( arg_input , mesh_fields , arg_begin , arg_state );

    Kokkos::parallel_for( arg_count , op );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i ) const
  {
    const size_type inode = inode_base + i ;

    size_type j = i * value_count ;

    displacement( inode , 0 , state_next ) = input[j++] ;
    displacement( inode , 1 , state_next ) = input[j++] ;
    displacement( inode , 2 , state_next ) = input[j++] ;
    velocity( inode , 0 , state_next ) = input[j++] ;
    velocity( inode , 1 , state_next ) = input[j++] ;
    velocity( inode , 2 , state_next ) = input[j++] ;
  }
};

} /* namespace Explicit */

#endif /* #ifndef KOKKOS_EXPLICITFUNCTORS_HPP */


