/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reax_c.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

void Temperature_Control( control_params *control, simulation_data *data )
{
  real tmp;

  if( control->T_mode == 1 ) {// step-wise temperature control
    if((data->step-data->prev_steps) % ((int)(control->T_freq/control->dt))==0){
      if( fabs( control->T - control->T_final ) >= fabs( control->T_rate ) )
        control->T += control->T_rate;
      else control->T = control->T_final;
    }
  }
  else if( control->T_mode == 2 ) { // constant slope control
    tmp = control->T_rate * control->dt / control->T_freq;

    if( fabs( control->T - control->T_final ) >= fabs( tmp ) )
      control->T += tmp;
  }
}



void Compute_Kinetic_Energy( reax_system* system, simulation_data* data,
                             MPI_Comm comm )
{
  int i;
  rvec p;
  real m;

  data->my_en.e_kin = 0.0;
  data->sys_en.e_kin = 0.0;
  data->therm.T = 0;

  for( i = 0; i < system->n; i++ ) {
    m = system->reax_param.sbp[system->my_atoms[i].type].mass;

    rvec_Scale( p, m, system->my_atoms[i].v );
    data->my_en.e_kin += 0.5 * rvec_Dot( p, system->my_atoms[i].v );
  }

  MPI_Allreduce( &data->my_en.e_kin,  &data->sys_en.e_kin,
                 1, MPI_DOUBLE, MPI_SUM, comm );

  data->therm.T = (2. * data->sys_en.e_kin) / (data->N_f * K_B);

  // avoid T being an absolute zero, might cause F.P.E!
  if( fabs(data->therm.T) < ALMOST_ZERO )
    data->therm.T = ALMOST_ZERO;
}


void Compute_System_Energy( reax_system *system, simulation_data *data,
                            MPI_Comm comm )
{
  real my_en[15], sys_en[15];

  my_en[0] = data->my_en.e_bond;
  my_en[1] = data->my_en.e_ov;
  my_en[2] = data->my_en.e_un;
  my_en[3] = data->my_en.e_lp;
  my_en[4] = data->my_en.e_ang;
  my_en[5] = data->my_en.e_pen;
  my_en[6] = data->my_en.e_coa;
  my_en[7] = data->my_en.e_hb;
  my_en[8] = data->my_en.e_tor;
  my_en[9] = data->my_en.e_con;
  my_en[10] = data->my_en.e_vdW;
  my_en[11] = data->my_en.e_ele;
  my_en[12] = data->my_en.e_pol;
  my_en[13] = data->my_en.e_kin;
  MPI_Reduce( my_en, sys_en, 14, MPI_DOUBLE, MPI_SUM, MASTER_NODE, comm );

  data->my_en.e_pot = data->my_en.e_bond +
    data->my_en.e_ov + data->my_en.e_un  + data->my_en.e_lp +
    data->my_en.e_ang + data->my_en.e_pen + data->my_en.e_coa +
    data->my_en.e_hb +
    data->my_en.e_tor + data->my_en.e_con +
    data->my_en.e_vdW + data->my_en.e_ele + data->my_en.e_pol;

  data->my_en.e_tot = data->my_en.e_pot + E_CONV * data->my_en.e_kin;

  if( system->my_rank == MASTER_NODE ) {
    data->sys_en.e_bond = sys_en[0];
    data->sys_en.e_ov = sys_en[1];
    data->sys_en.e_un = sys_en[2];
    data->sys_en.e_lp = sys_en[3];
    data->sys_en.e_ang = sys_en[4];
    data->sys_en.e_pen = sys_en[5];
    data->sys_en.e_coa = sys_en[6];
    data->sys_en.e_hb = sys_en[7];
    data->sys_en.e_tor = sys_en[8];
    data->sys_en.e_con = sys_en[9];
    data->sys_en.e_vdW = sys_en[10];
    data->sys_en.e_ele = sys_en[11];
    data->sys_en.e_pol = sys_en[12];
    data->sys_en.e_kin = sys_en[13];

    data->sys_en.e_pot = data->sys_en.e_bond +
      data->sys_en.e_ov + data->sys_en.e_un  + data->sys_en.e_lp +
      data->sys_en.e_ang + data->sys_en.e_pen + data->sys_en.e_coa +
      data->sys_en.e_hb +
      data->sys_en.e_tor + data->sys_en.e_con +
      data->sys_en.e_vdW + data->sys_en.e_ele + data->sys_en.e_pol;

    data->sys_en.e_tot = data->sys_en.e_pot + E_CONV * data->sys_en.e_kin;
  }
}


void Compute_Total_Mass( reax_system *system, simulation_data *data,
                         MPI_Comm comm  )
{
  int  i;
  real tmp;

  tmp = 0;
  for( i = 0; i < system->n; i++ )
    tmp += system->reax_param.sbp[ system->my_atoms[i].type ].mass;

  MPI_Allreduce( &tmp, &data->M, 1, MPI_DOUBLE, MPI_SUM, comm );

  data->inv_M = 1. / data->M;
}



void Compute_Center_of_Mass( reax_system *system, simulation_data *data,
                             mpi_datatypes *mpi_data, MPI_Comm comm )
{
  int i;
  real m, det; //xx, xy, xz, yy, yz, zz;
  real tmp_mat[6], tot_mat[6];
  rvec my_xcm, my_vcm, my_amcm, my_avcm;
  rvec tvec, diff;
  rtensor mat, inv;

  rvec_MakeZero( my_xcm );  // position of CoM
  rvec_MakeZero( my_vcm );  // velocity of CoM
  rvec_MakeZero( my_amcm ); // angular momentum of CoM
  rvec_MakeZero( my_avcm ); // angular velocity of CoM

  /* Compute the position, vel. and ang. momentum about the centre of mass */
  for( i = 0; i < system->n; ++i ) {
    m = system->reax_param.sbp[ system->my_atoms[i].type ].mass;

    rvec_ScaledAdd( my_xcm, m, system->my_atoms[i].x );
    rvec_ScaledAdd( my_vcm, m, system->my_atoms[i].v );

    rvec_Cross( tvec, system->my_atoms[i].x, system->my_atoms[i].v );
    rvec_ScaledAdd( my_amcm, m, tvec );
  }

  MPI_Allreduce( my_xcm, data->xcm, 3, MPI_DOUBLE, MPI_SUM, comm );
  MPI_Allreduce( my_vcm, data->vcm, 3, MPI_DOUBLE, MPI_SUM, comm );
  MPI_Allreduce( my_amcm, data->amcm, 3, MPI_DOUBLE, MPI_SUM, comm );

  rvec_Scale( data->xcm, data->inv_M, data->xcm );
  rvec_Scale( data->vcm, data->inv_M, data->vcm );
  rvec_Cross( tvec, data->xcm, data->vcm );
  rvec_ScaledAdd( data->amcm, -data->M, tvec );
  data->etran_cm = 0.5 * data->M * rvec_Norm_Sqr( data->vcm );

  /* Calculate and then invert the inertial tensor */
  for( i = 0; i < 6; ++i )
    tmp_mat[i] = 0;
  //my_xx = my_xy = my_xz = my_yy = my_yz = my_zz = 0;

  for( i = 0; i < system->n; ++i ){
    m = system->reax_param.sbp[ system->my_atoms[i].type ].mass;
    rvec_ScaledSum( diff, 1., system->my_atoms[i].x, -1., data->xcm );

    tmp_mat[0]/*my_xx*/ += diff[0] * diff[0] * m;
    tmp_mat[1]/*my_xy*/ += diff[0] * diff[1] * m;
    tmp_mat[2]/*my_xz*/ += diff[0] * diff[2] * m;
    tmp_mat[3]/*my_yy*/ += diff[1] * diff[1] * m;
    tmp_mat[4]/*my_yz*/ += diff[1] * diff[2] * m;
    tmp_mat[5]/*my_zz*/ += diff[2] * diff[2] * m;
  }

  MPI_Reduce( tmp_mat, tot_mat, 6, MPI_DOUBLE, MPI_SUM, MASTER_NODE, comm );

  if( system->my_rank == MASTER_NODE ) {
    mat[0][0] = tot_mat[3] + tot_mat[5];  // yy + zz;
    mat[0][1] = mat[1][0] = -tot_mat[1];  // -xy;
    mat[0][2] = mat[2][0] = -tot_mat[2];  // -xz;
    mat[1][1] = tot_mat[0] + tot_mat[5];  // xx + zz;
    mat[2][1] = mat[1][2] = -tot_mat[4];  // -yz;
    mat[2][2] = tot_mat[0] + tot_mat[3];  // xx + yy;

    /* invert the inertial tensor */
    det = ( mat[0][0] * mat[1][1] * mat[2][2] +
            mat[0][1] * mat[1][2] * mat[2][0] +
            mat[0][2] * mat[1][0] * mat[2][1] ) -
      ( mat[0][0] * mat[1][2] * mat[2][1] +
        mat[0][1] * mat[1][0] * mat[2][2] +
        mat[0][2] * mat[1][1] * mat[2][0] );

    inv[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    inv[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
    inv[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    inv[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    inv[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
    inv[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
    inv[2][0] = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1];
    inv[2][1] = mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1];
    inv[2][2] = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];

    if( det > ALMOST_ZERO )
      rtensor_Scale( inv, 1./det, inv );
    else rtensor_MakeZero( inv );

    /* Compute the angular velocity about the centre of mass */
    rtensor_MatVec( data->avcm, inv, data->amcm );
  }

  MPI_Bcast( data->avcm, 3, MPI_DOUBLE, MASTER_NODE, comm );

  /* Compute the rotational energy */
  data->erot_cm = 0.5 * E_CONV * rvec_Dot( data->avcm, data->amcm );

}

void Compute_Pressure(reax_system* system, control_params *control,
                      simulation_data* data, mpi_datatypes *mpi_data)
{
  int i;
  reax_atom *p_atom;
  rvec tmp, tx, int_press;
  simulation_box *big_box = &(system->big_box);

  /* Calculate internal pressure */
  rvec_MakeZero( int_press );

  // 0: both int and ext, 1: ext only, 2: int only
  if( control->press_mode == 0 || control->press_mode == 2 ) {
    for( i = 0; i < system->n; ++i ) {
      p_atom = &( system->my_atoms[i] );

      /* transform x into unitbox coordinates */
      Transform_to_UnitBox( p_atom->x, big_box, 1, tx );

      /* this atom's contribution to internal pressure */
      rvec_Multiply( tmp, p_atom->f, tx );
      rvec_Add( int_press, tmp );

    }
  }

  /* sum up internal and external pressure */
  MPI_Allreduce( int_press, data->int_press,
                 3, MPI_DOUBLE, MPI_SUM, mpi_data->comm_mesh3D );
  MPI_Allreduce( data->my_ext_press, data->ext_press,
                 3, MPI_DOUBLE, MPI_SUM, mpi_data->comm_mesh3D );

  /* kinetic contribution */
  data->kin_press = 2.*(E_CONV*data->sys_en.e_kin) / (3.*big_box->V*P_CONV);

  /* Calculate total pressure in each direction */
  data->tot_press[0] = data->kin_press -
    (( data->int_press[0] + data->ext_press[0] ) /
     ( big_box->box_norms[1] * big_box->box_norms[2] * P_CONV ));

  data->tot_press[1] = data->kin_press -
    (( data->int_press[1] + data->ext_press[1] ) /
     ( big_box->box_norms[0] * big_box->box_norms[2] * P_CONV ));

  data->tot_press[2] = data->kin_press -
    (( data->int_press[2] + data->ext_press[2] ) /
     ( big_box->box_norms[0] * big_box->box_norms[1] * P_CONV ));

  data->iso_bar.P =
    ( data->tot_press[0] + data->tot_press[1] + data->tot_press[2] ) / 3.;
}
