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

#include "reaxc_valence_angles.h"
#include <cmath>
#include "pair.h"
#include "reaxc_defs.h"
#include "reaxc_list.h"
#include "reaxc_vector.h"

#include "error.h"

static double Dot( double* v1, double* v2, int k )
{
  double ret = 0.0;

  for( int i=0; i < k; ++i )
    ret +=  v1[i] * v2[i];

  return ret;
}

void Calculate_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                      double *theta, double *cos_theta )
{
  (*cos_theta) = Dot( dvec_ji, dvec_jk, 3 ) / ( d_ji * d_jk );
  if (*cos_theta > 1.) *cos_theta  = 1.0;
  if (*cos_theta < -1.) *cos_theta  = -1.0;

  (*theta) = acos( *cos_theta );
}

void Calculate_dCos_Theta( rvec dvec_ji, double d_ji, rvec dvec_jk, double d_jk,
                           rvec* dcos_theta_di,
                           rvec* dcos_theta_dj,
                           rvec* dcos_theta_dk )
{
  int t;
  double sqr_d_ji = SQR(d_ji);
  double sqr_d_jk = SQR(d_jk);
  double inv_dists = 1.0 / (d_ji * d_jk);
  double inv_dists3 = pow( inv_dists, 3.0 );
  double dot_dvecs = Dot( dvec_ji, dvec_jk, 3 );
  double Cdot_inv3 = dot_dvecs * inv_dists3;

  for( t = 0; t < 3; ++t ) {
    (*dcos_theta_di)[t] = dvec_jk[t] * inv_dists -
      Cdot_inv3 * sqr_d_jk * dvec_ji[t];
    (*dcos_theta_dj)[t] = -(dvec_jk[t] + dvec_ji[t]) * inv_dists +
      Cdot_inv3 * ( sqr_d_jk * dvec_ji[t] + sqr_d_ji * dvec_jk[t] );
    (*dcos_theta_dk)[t] = dvec_ji[t] * inv_dists -
      Cdot_inv3 * sqr_d_ji * dvec_jk[t];
  }
}


void Valence_Angles( reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, output_controls * /*out_control*/ )
{
  int i, j, pi, k, pk, t;
  int type_i, type_j, type_k;
  int start_j, end_j, start_pk, end_pk;
  int cnt, num_thb_intrs;

  double temp, temp_bo_jt, pBOjt7;
  double p_val1, p_val2, p_val3, p_val4, p_val5;
  double p_val6, p_val7, p_val8, p_val9, p_val10;
  double p_pen1, p_pen2, p_pen3, p_pen4;
  double p_coa1, p_coa2, p_coa3, p_coa4;
  double trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
  double exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
  double dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
  double CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
  double CEpen1, CEpen2, CEpen3;
  double e_ang, e_coa, e_pen;
  double CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
  double Cf7ij, Cf7jk, Cf8j, Cf9j;
  double f7_ij, f7_jk, f8_Dj, f9_Dj;
  double Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta;
  double BOA_ij, BOA_jk;
  rvec force, ext_press;

  // Tallying variables
  double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
  double delij[3], delkj[3];

  three_body_header *thbh;
  three_body_parameters *thbp;
  three_body_interaction_data *p_ijk, *p_kji;
  bond_data *pbond_ij, *pbond_jk, *pbond_jt;
  bond_order_data *bo_ij, *bo_jk, *bo_jt;
  reax_list *bonds = (*lists) + BONDS;
  reax_list *thb_intrs =  (*lists) + THREE_BODIES;

  /* global parameters used in these calculations */
  p_val6 = system->reax_param.gp.l[14];
  p_val8 = system->reax_param.gp.l[33];
  p_val9 = system->reax_param.gp.l[16];
  p_val10 = system->reax_param.gp.l[17];
  num_thb_intrs = 0;


  for( j = 0; j < system->N; ++j ) {         // Ray: the first one with system->N
    type_j = system->my_atoms[j].type;
    if (type_j < 0) continue;
    start_j = Start_Index(j, bonds);
    end_j = End_Index(j, bonds);

    p_val3 = system->reax_param.sbp[ type_j ].p_val3;
    p_val5 = system->reax_param.sbp[ type_j ].p_val5;

    SBOp = 0, prod_SBO = 1;
    for( t = start_j; t < end_j; ++t ) {
      bo_jt = &(bonds->select.bond_list[t].bo_data);
      SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
      temp = SQR( bo_jt->BO );
      temp *= temp;
      temp *= temp;
      prod_SBO *= exp( -temp );
    }

    if (workspace->vlpex[j] >= 0) {
      vlpadj = 0;
      dSBO2 = prod_SBO - 1;
    } else {
      vlpadj = workspace->nlp[j];
      dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
    }

    SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
    dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );

    if (SBO <= 0)
      SBO2 = 0, CSBO2 = 0;
    else if (SBO > 0 && SBO <= 1) {
        SBO2 = pow( SBO, p_val9 );
        CSBO2 = p_val9 * pow( SBO, p_val9 - 1 );
    }
    else if (SBO > 1 && SBO < 2) {
      SBO2 = 2 - pow( 2-SBO, p_val9 );
      CSBO2 = p_val9 * pow( 2 - SBO, p_val9 - 1 );
    }
    else
      SBO2 = 2, CSBO2 = 0;

    expval6 = exp( p_val6 * workspace->Delta_boc[j] );

    for( pi = start_j; pi < end_j; ++pi ) {
      Set_Start_Index( pi, num_thb_intrs, thb_intrs );
      pbond_ij = &(bonds->select.bond_list[pi]);
      bo_ij = &(pbond_ij->bo_data);
      BOA_ij = bo_ij->BO - control->thb_cut;


      if( BOA_ij/*bo_ij->BO*/ > 0.0 &&
          ( j < system->n || pbond_ij->nbr < system->n ) ) {
        i = pbond_ij->nbr;
        type_i = system->my_atoms[i].type;

        for( pk = start_j; pk < pi; ++pk ) {
          start_pk = Start_Index( pk, thb_intrs );
          end_pk = End_Index( pk, thb_intrs );

          for( t = start_pk; t < end_pk; ++t )
            if (thb_intrs->select.three_body_list[t].thb == i) {
              p_ijk = &(thb_intrs->select.three_body_list[num_thb_intrs] );
              p_kji = &(thb_intrs->select.three_body_list[t]);

              p_ijk->thb = bonds->select.bond_list[pk].nbr;
              p_ijk->pthb  = pk;
              p_ijk->theta = p_kji->theta;
              rvec_Copy( p_ijk->dcos_di, p_kji->dcos_dk );
              rvec_Copy( p_ijk->dcos_dj, p_kji->dcos_dj );
              rvec_Copy( p_ijk->dcos_dk, p_kji->dcos_di );

              ++num_thb_intrs;
              break;
            }
        }

        for( pk = pi+1; pk < end_j; ++pk ) {
          pbond_jk = &(bonds->select.bond_list[pk]);
          bo_jk    = &(pbond_jk->bo_data);
          BOA_jk   = bo_jk->BO - control->thb_cut;
          k        = pbond_jk->nbr;
          type_k   = system->my_atoms[k].type;
          p_ijk    = &( thb_intrs->select.three_body_list[num_thb_intrs] );

          Calculate_Theta( pbond_ij->dvec, pbond_ij->d,
                           pbond_jk->dvec, pbond_jk->d,
                           &theta, &cos_theta );

          Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d,
                                pbond_jk->dvec, pbond_jk->d,
                                &(p_ijk->dcos_di), &(p_ijk->dcos_dj),
                                &(p_ijk->dcos_dk) );
          p_ijk->thb = k;
          p_ijk->pthb = pk;
          p_ijk->theta = theta;

          sin_theta = sin( theta );
          if (sin_theta < 1.0e-5)
            sin_theta = 1.0e-5;

          ++num_thb_intrs;


          if ((j < system->n) && (BOA_jk > 0.0) &&
              (bo_ij->BO > control->thb_cut) &&
              (bo_jk->BO > control->thb_cut) &&
              (bo_ij->BO * bo_jk->BO > control->thb_cutsq) ) {
            thbh = &( system->reax_param.thbp[ type_i ][ type_j ][ type_k ] );

            for( cnt = 0; cnt < thbh->cnt; ++cnt ) {
              if (fabs(thbh->prm[cnt].p_val1) > 0.001) {
                thbp = &( thbh->prm[cnt] );

                /* ANGLE ENERGY */
                p_val1 = thbp->p_val1;
                p_val2 = thbp->p_val2;
                p_val4 = thbp->p_val4;
                p_val7 = thbp->p_val7;
                theta_00 = thbp->theta_00;

                exp3ij = exp( -p_val3 * pow( BOA_ij, p_val4 ) );
                f7_ij = 1.0 - exp3ij;
                Cf7ij = p_val3 * p_val4 * pow( BOA_ij, p_val4 - 1.0 ) * exp3ij;

                exp3jk = exp( -p_val3 * pow( BOA_jk, p_val4 ) );
                f7_jk = 1.0 - exp3jk;
                Cf7jk = p_val3 * p_val4 * pow( BOA_jk, p_val4 - 1.0 ) * exp3jk;

                expval7 = exp( -p_val7 * workspace->Delta_boc[j] );
                trm8 = 1.0 + expval6 + expval7;
                f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
                Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *
                  ( p_val6 * expval6 * trm8 -
                    (2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );

                theta_0 = 180.0 - theta_00 * (1.0 -
                                              exp(-p_val10 * (2.0 - SBO2)));
                theta_0 = DEG2RAD( theta_0 );

                expval2theta  = exp( -p_val2 * SQR(theta_0 - theta) );
                if (p_val1 >= 0)
                  expval12theta = p_val1 * (1.0 - expval2theta);
                else // To avoid linear Me-H-Me angles (6/6/06)
                  expval12theta = p_val1 * -expval2theta;

                CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
                CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
                CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
                CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj *
                  expval2theta * (theta_0 - theta);

                Ctheta_0 = p_val10 * DEG2RAD(theta_00) *
                  exp( -p_val10 * (2.0 - SBO2) );

                CEval5 = -CEval4 * Ctheta_0 * CSBO2;
                CEval6 = CEval5 * dSBO1;
                CEval7 = CEval5 * dSBO2;
                CEval8 = -CEval4 / sin_theta;

                data->my_en.e_ang += e_ang =
                  f7_ij * f7_jk * f8_Dj * expval12theta;
                /* END ANGLE ENERGY*/

                /* PENALTY ENERGY */
                p_pen1 = thbp->p_pen1;
                p_pen2 = system->reax_param.gp.l[19];
                p_pen3 = system->reax_param.gp.l[20];
                p_pen4 = system->reax_param.gp.l[21];

                exp_pen2ij = exp( -p_pen2 * SQR( BOA_ij - 2.0 ) );
                exp_pen2jk = exp( -p_pen2 * SQR( BOA_jk - 2.0 ) );
                exp_pen3 = exp( -p_pen3 * workspace->Delta[j] );
                exp_pen4 = exp(  p_pen4 * workspace->Delta[j] );
                trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
                f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
                Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 -
                         (2.0 + exp_pen3) * ( -p_pen3 * exp_pen3 +
                                              p_pen4 * exp_pen4 ) ) /
                  SQR( trm_pen34 );

                data->my_en.e_pen += e_pen =
                  p_pen1 * f9_Dj * exp_pen2ij * exp_pen2jk;

                CEpen1 = e_pen * Cf9j / f9_Dj;
                temp   = -2.0 * p_pen2 * e_pen;
                CEpen2 = temp * (BOA_ij - 2.0);
                CEpen3 = temp * (BOA_jk - 2.0);
                /* END PENALTY ENERGY */

                /* COALITION ENERGY */
                p_coa1 = thbp->p_coa1;
                p_coa2 = system->reax_param.gp.l[2];
                p_coa3 = system->reax_param.gp.l[38];
                p_coa4 = system->reax_param.gp.l[30];

                exp_coa2 = exp( p_coa2 * workspace->Delta_val[j] );
                data->my_en.e_coa += e_coa =
                  p_coa1 / (1. + exp_coa2) *
                  exp( -p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij) ) *
                  exp( -p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk) ) *
                  exp( -p_coa4 * SQR(BOA_ij - 1.5) ) *
                  exp( -p_coa4 * SQR(BOA_jk - 1.5) );

                CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
                CEcoa2 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
                CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
                CEcoa4 = -2 * p_coa3 *
                  (workspace->total_bond_order[i]-BOA_ij) * e_coa;
                CEcoa5 = -2 * p_coa3 *
                  (workspace->total_bond_order[k]-BOA_jk) * e_coa;
                /* END COALITION ENERGY */

                /* FORCES */
                bo_ij->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
                bo_jk->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));
                workspace->CdDelta[j] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
                workspace->CdDelta[i] += CEcoa4;
                workspace->CdDelta[k] += CEcoa5;

                for( t = start_j; t < end_j; ++t ) {
                    pbond_jt = &( bonds->select.bond_list[t] );
                    bo_jt = &(pbond_jt->bo_data);
                    temp_bo_jt = bo_jt->BO;
                    temp = CUBE( temp_bo_jt );
                    pBOjt7 = temp * temp * temp_bo_jt;

                    bo_jt->Cdbo += (CEval6 * pBOjt7);
                    bo_jt->Cdbopi += CEval5;
                    bo_jt->Cdbopi2 += CEval5;
                }

                if (control->virial == 0) {
                  rvec_ScaledAdd( workspace->f[i], CEval8, p_ijk->dcos_di );
                  rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );
                  rvec_ScaledAdd( workspace->f[k], CEval8, p_ijk->dcos_dk );
                } else {
                  rvec_Scale( force, CEval8, p_ijk->dcos_di );
                  rvec_Add( workspace->f[i], force );
                  rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
                  rvec_Add( data->my_ext_press, ext_press );

                  rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );

                  rvec_Scale( force, CEval8, p_ijk->dcos_dk );
                  rvec_Add( workspace->f[k], force );
                  rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
                  rvec_Add( data->my_ext_press, ext_press );
                }

                /* tally into per-atom virials */
                if (system->pair_ptr->vflag_atom || system->pair_ptr->evflag) {

                  /* Acquire vectors */
                  rvec_ScaledSum( delij, 1., system->my_atoms[i].x,
                                        -1., system->my_atoms[j].x );
                  rvec_ScaledSum( delkj, 1., system->my_atoms[k].x,
                                        -1., system->my_atoms[j].x );

                  rvec_Scale( fi_tmp, -CEval8, p_ijk->dcos_di );
                  rvec_Scale( fj_tmp, -CEval8, p_ijk->dcos_dj );
                  rvec_Scale( fk_tmp, -CEval8, p_ijk->dcos_dk );

                  eng_tmp = e_ang + e_pen + e_coa;

                  if (system->pair_ptr->evflag)
                          system->pair_ptr->ev_tally(j,j,system->N,1,eng_tmp,0.0,0.0,0.0,0.0,0.0);
                  if (system->pair_ptr->vflag_atom)
                          system->pair_ptr->v_tally3(i,j,k,fi_tmp,fk_tmp,delij,delkj);
                }
              }
            }
          }
        }
      }

      Set_End_Index(pi, num_thb_intrs, thb_intrs );
    }
  }

  if (num_thb_intrs >= thb_intrs->num_intrs * DANGER_ZONE) {
    workspace->realloc.num_3body = num_thb_intrs;
    if (num_thb_intrs > thb_intrs->num_intrs) {
      char errmsg[128];
      snprintf(errmsg, 128, "step%d-ran out of space on angle_list: top=%d, max=%d",
               data->step, num_thb_intrs, thb_intrs->num_intrs);
      control->error_ptr->one(FLERR, errmsg);
    }
  }

}
