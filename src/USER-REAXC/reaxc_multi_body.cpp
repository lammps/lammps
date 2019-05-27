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

#include "pair_reaxc.h"
#include "reaxc_multi_body.h"
#include "reaxc_bond_orders.h"
#include "reaxc_list.h"
#include "reaxc_vector.h"

void Atom_Energy( reax_system *system, control_params *control,
                  simulation_data *data, storage *workspace, reax_list **lists,
                  output_controls * /*out_control*/ )
{
  int i, j, pj, type_i, type_j;
  double Delta_lpcorr, dfvl;
  double e_lp, expvd2, inv_expvd2, dElp, CElp, DlpVi;
  double e_lph, Di, vov3, deahu2dbo, deahu2dsbo;
  double e_ov, CEover1, CEover2, CEover3, CEover4;
  double exp_ovun1, exp_ovun2, sum_ovun1, sum_ovun2;
  double exp_ovun2n, exp_ovun6, exp_ovun8;
  double inv_exp_ovun1, inv_exp_ovun2, inv_exp_ovun2n, inv_exp_ovun8;
  double e_un, CEunder1, CEunder2, CEunder3, CEunder4;
  double p_lp2, p_lp3;
  double p_ovun2, p_ovun3, p_ovun4, p_ovun5, p_ovun6, p_ovun7, p_ovun8;
  double eng_tmp;
  int numbonds;

  single_body_parameters *sbp_i;
  two_body_parameters *twbp;
  bond_data *pbond;
  bond_order_data *bo_ij;
  reax_list *bonds = (*lists) + BONDS;

  /* Initialize parameters */
  p_lp3 = system->reax_param.gp.l[5];
  p_ovun3 = system->reax_param.gp.l[32];
  p_ovun4 = system->reax_param.gp.l[31];
  p_ovun6 = system->reax_param.gp.l[6];
  p_ovun7 = system->reax_param.gp.l[8];
  p_ovun8 = system->reax_param.gp.l[9];

  for( i = 0; i < system->n; ++i ) {
    /* set the parameter pointer */
    type_i = system->my_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[ type_i ]);

    /* lone-pair Energy */
    p_lp2 = sbp_i->p_lp2;
    expvd2 = exp( -75 * workspace->Delta_lp[i] );
    inv_expvd2 = 1. / (1. + expvd2 );

    numbonds = 0;
    e_lp = 0.0;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
      numbonds ++;

    /* calculate the energy */
    if (numbonds > 0 || control->enobondsflag)
      data->my_en.e_lp += e_lp =
        p_lp2 * workspace->Delta_lp[i] * inv_expvd2;

    dElp = p_lp2 * inv_expvd2 +
      75 * p_lp2 * workspace->Delta_lp[i] * expvd2 * SQR(inv_expvd2);
    CElp = dElp * workspace->dDelta_lp[i];

    if (numbonds > 0 || control->enobondsflag)
      workspace->CdDelta[i] += CElp;  // lp - 1st term

    /* tally into per-atom energy */
    if (system->pair_ptr->evflag)
      system->pair_ptr->ev_tally(i,i,system->n,1,e_lp,0.0,0.0,0.0,0.0,0.0);

    /* correction for C2 */
    if (p_lp3 > 0.001 && !strcmp(system->reax_param.sbp[type_i].name, "C"))
      for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->my_atoms[j].type;
        if (type_j < 0) continue;

        if (!strcmp( system->reax_param.sbp[type_j].name, "C" )) {
          twbp = &( system->reax_param.tbp[type_i][type_j]);
          bo_ij = &( bonds->select.bond_list[pj].bo_data );
          Di = workspace->Delta[i];
          vov3 = bo_ij->BO - Di - 0.040*pow(Di, 4.);

          if (vov3 > 3.) {
            data->my_en.e_lp += e_lph = p_lp3 * SQR(vov3-3.0);

            deahu2dbo = 2.*p_lp3*(vov3 - 3.);
            deahu2dsbo = 2.*p_lp3*(vov3 - 3.)*(-1. - 0.16*pow(Di, 3.));

            bo_ij->Cdbo += deahu2dbo;
            workspace->CdDelta[i] += deahu2dsbo;

            /* tally into per-atom energy */
            if (system->pair_ptr->evflag)
              system->pair_ptr->ev_tally(i,j,system->n,1,e_lph,0.0,0.0,0.0,0.0,0.0);

          }
        }
      }
  }


  for( i = 0; i < system->n; ++i ) {
    type_i = system->my_atoms[i].type;
    if (type_i < 0) continue;
    sbp_i = &(system->reax_param.sbp[ type_i ]);

    /* over-coordination energy */
    if (sbp_i->mass > 21.0)
      dfvl = 0.0;
    else dfvl = 1.0; // only for 1st-row elements

    p_ovun2 = sbp_i->p_ovun2;
    sum_ovun1 = sum_ovun2 = 0;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->my_atoms[j].type;
        if (type_j < 0) continue;
        bo_ij = &(bonds->select.bond_list[pj].bo_data);
        twbp = &(system->reax_param.tbp[ type_i ][ type_j ]);

        sum_ovun1 += twbp->p_ovun1 * twbp->De_s * bo_ij->BO;
        sum_ovun2 += (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j])*
          ( bo_ij->BO_pi + bo_ij->BO_pi2 );

      }

    exp_ovun1 = p_ovun3 * exp( p_ovun4 * sum_ovun2 );
    inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
    Delta_lpcorr  = workspace->Delta[i] -
      (dfvl * workspace->Delta_lp_temp[i]) * inv_exp_ovun1;

    exp_ovun2 = exp( p_ovun2 * Delta_lpcorr );
    inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);

    DlpVi = 1.0 / (Delta_lpcorr + sbp_i->valency + 1e-8);
    CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;

    data->my_en.e_ov += e_ov = sum_ovun1 * CEover1;

    CEover2 = sum_ovun1 * DlpVi * inv_exp_ovun2 *
      (1.0 - Delta_lpcorr * ( DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2 ));

    CEover3 = CEover2 * (1.0 - dfvl * workspace->dDelta_lp[i] * inv_exp_ovun1 );

    CEover4 = CEover2 * (dfvl * workspace->Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);


    /* under-coordination potential */
    p_ovun2 = sbp_i->p_ovun2;
    p_ovun5 = sbp_i->p_ovun5;

    exp_ovun2n = 1.0 / exp_ovun2;
    exp_ovun6 = exp( p_ovun6 * Delta_lpcorr );
    exp_ovun8 = p_ovun7 * exp(p_ovun8 * sum_ovun2);
    inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
    inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

    numbonds = 0;
    e_un = 0.0;
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
      numbonds ++;

    if (numbonds > 0 || control->enobondsflag)
      data->my_en.e_un += e_un =
        -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;

    CEunder1 = inv_exp_ovun2n *
      ( p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
        p_ovun2 * e_un * exp_ovun2n );
    CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
    CEunder3 = CEunder1 * (1.0 - dfvl*workspace->dDelta_lp[i]*inv_exp_ovun1);
    CEunder4 = CEunder1 * (dfvl*workspace->Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1) + CEunder2;

    /* tally into per-atom energy */
    if (system->pair_ptr->evflag) {
      eng_tmp = e_ov;
      if (numbonds > 0 || control->enobondsflag)
        eng_tmp += e_un;
      system->pair_ptr->ev_tally(i,i,system->n,1,eng_tmp,0.0,0.0,0.0,0.0,0.0);
    }

    /* forces */
    workspace->CdDelta[i] += CEover3;   // OvCoor - 2nd term
    if (numbonds > 0 || control->enobondsflag)
      workspace->CdDelta[i] += CEunder3;  // UnCoor - 1st term

    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      pbond = &(bonds->select.bond_list[pj]);
      j = pbond->nbr;
      bo_ij = &(pbond->bo_data);
      twbp  = &(system->reax_param.tbp[ system->my_atoms[i].type ]
                [system->my_atoms[pbond->nbr].type]);


      bo_ij->Cdbo += CEover1 * twbp->p_ovun1 * twbp->De_s;// OvCoor-1st
      workspace->CdDelta[j] += CEover4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
        (bo_ij->BO_pi + bo_ij->BO_pi2); // OvCoor-3a
      bo_ij->Cdbopi += CEover4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]); // OvCoor-3b
      bo_ij->Cdbopi2 += CEover4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // OvCoor-3b


      workspace->CdDelta[j] += CEunder4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
        (bo_ij->BO_pi + bo_ij->BO_pi2);   // UnCoor - 2a
      bo_ij->Cdbopi += CEunder4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b
      bo_ij->Cdbopi2 += CEunder4 *
        (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b

    }

  }
}
