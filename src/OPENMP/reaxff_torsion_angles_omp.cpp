// clang-format off
/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program
  Website: https://www.cs.purdue.edu/puremd

  Copyright (2010) Purdue University

  Contributing authors:
  H. M. Aktulga, J. Fogarty, S. Pandit, A. Grama
  Corresponding author:
  Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, 38 (4-5), 245-259

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxff_omp.h"

#include "fix_omp.h"
#include "pair_reaxff_omp.h"
#include "reaxff_api.h"

#include <cmath>

#define MIN_SINE 1e-10

using namespace LAMMPS_NS;

namespace ReaxFF {

  /* ---------------------------------------------------------------------- */
  void Torsion_AnglesOMP(reax_system *system, control_params *control,
                         simulation_data *data, storage *workspace,
                         reax_list **lists)
  {
    int natoms = system->n;
    reax_list *bonds = (*lists) + BONDS;
    reax_list *thb_intrs = (*lists) + THREE_BODIES;
    double p_tor2 = system->reax_param.gp.l[23];
    double p_tor3 = system->reax_param.gp.l[24];
    double p_tor4 = system->reax_param.gp.l[25];
    double p_cot2 = system->reax_param.gp.l[27];
    double total_Etor = 0;
    double total_Econ = 0;
    int  nthreads = control->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(shared) reduction(+: total_Etor, total_Econ)
#endif
    {
      int i, j, k, l, pi, pj, pk, pl, pij, plk;
      int type_i, type_j, type_k, type_l;
      int start_j, end_j;
      int start_pj, end_pj, start_pk, end_pk;
      int num_frb_intrs = 0;

      double Delta_j, Delta_k;
      double r_ij, r_jk, r_kl, r_li;
      double BOA_ij, BOA_jk, BOA_kl;

      double exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
      double exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
      double exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
      double fn10, f11_DjDk, dfn11, fn12;
      double theta_ijk, theta_jkl;
      double sin_ijk, sin_jkl;
      double cos_ijk, cos_jkl;
      double tan_ijk_i, tan_jkl_i;
      double omega, cos_omega, cos2omega, cos3omega;
      rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
      double CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
      double CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
      double Cconj, CEconj1, CEconj2, CEconj3;
      double CEconj4, CEconj5, CEconj6;
      double e_tor, e_con;
      rvec dvec_li;
      four_body_header *fbh;
      four_body_parameters *fbp;
      bond_data *pbond_ij, *pbond_jk, *pbond_kl;
      bond_order_data *bo_ij, *bo_jk, *bo_kl;
      three_body_interaction_data *p_ijk, *p_jkl;

      // Virial tallying variables
      double delil[3], deljl[3], delkl[3];
      double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];

      int tid = get_tid();

      long reductionOffset = (system->N * tid);
      class PairReaxFFOMP *pair_reax_ptr;
      pair_reax_ptr = static_cast<class PairReaxFFOMP*>(system->pair_ptr);
      class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
      for (j = 0; j < system->N; ++j) {
        start_j = Start_Index(j, bonds);
        end_j = End_Index(j, bonds);

        for (pk = start_j; pk < end_j; ++pk) {
          bo_jk = &(bonds->select.bond_list[pk].bo_data);
          for (k = 0; k < nthreads; ++k)
            bo_jk->CdboReduction[k] = 0.;
        }
      }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
      for (j = 0; j < natoms; ++j) {
        type_j = system->my_atoms[j].type;
        Delta_j = workspace->Delta_boc[j];
        start_j = Start_Index(j, bonds);
        end_j = End_Index(j, bonds);

        for (pk = start_j; pk < end_j; ++pk) {
          pbond_jk = &(bonds->select.bond_list[pk]);
          k = pbond_jk->nbr;
          bo_jk = &(pbond_jk->bo_data);
          BOA_jk = bo_jk->BO - control->thb_cut;

          /* see if there are any 3-body interactions involving j&k
             where j is the central atom. Otherwise there is no point in
             trying to form a 4-body interaction out of this neighborhood */
          if (system->my_atoms[j].orig_id < system->my_atoms[k].orig_id &&
              bo_jk->BO > control->thb_cut/*0*/ && Num_Entries(pk, thb_intrs)) {
            pj = pbond_jk->sym_index; // pj points to j on k's list

            /* do the same check as above:
               are there any 3-body interactions involving k&j
               where k is the central atom */
            if (Num_Entries(pj, thb_intrs)) {
              type_k = system->my_atoms[k].type;
              Delta_k = workspace->Delta_boc[k];
              r_jk = pbond_jk->d;

              start_pk = Start_Index(pk, thb_intrs);
              end_pk = End_Index(pk, thb_intrs);
              start_pj = Start_Index(pj, thb_intrs);
              end_pj = End_Index(pj, thb_intrs);

              exp_tor2_jk = exp(-p_tor2 * BOA_jk);
              exp_cot2_jk = exp(-p_cot2 * SQR(BOA_jk - 1.5));
              exp_tor3_DjDk = exp(-p_tor3 * (Delta_j + Delta_k));
              exp_tor4_DjDk = exp(p_tor4  * (Delta_j + Delta_k));
              exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
              f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;


              /* pick i up from j-k interaction where j is the central atom */
              for (pi = start_pk; pi < end_pk; ++pi) {
                p_ijk = &(thb_intrs->select.three_body_list[pi]);
                pij = p_ijk->pthb; // pij is pointer to i on j's bond_list
                pbond_ij = &(bonds->select.bond_list[pij]);
                bo_ij = &(pbond_ij->bo_data);

                if (bo_ij->BO > control->thb_cut/*0*/) {
                  i = p_ijk->thb;
                  type_i = system->my_atoms[i].type;
                  r_ij = pbond_ij->d;
                  BOA_ij = bo_ij->BO - control->thb_cut;

                  theta_ijk = p_ijk->theta;
                  sin_ijk = sin(theta_ijk);
                  cos_ijk = cos(theta_ijk);
                  //tan_ijk_i = 1. / tan(theta_ijk);
                  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE)
                    tan_ijk_i = cos_ijk / MIN_SINE;
                  else if (sin_ijk <= 0 && sin_ijk >= -MIN_SINE)
                    tan_ijk_i = cos_ijk / -MIN_SINE;
                  else tan_ijk_i = cos_ijk / sin_ijk;

                  exp_tor2_ij = exp(-p_tor2 * BOA_ij);
                  exp_cot2_ij = exp(-p_cot2 * SQR(BOA_ij -1.5));


                  /* pick l up from j-k interaction where k is the central atom */
                  for (pl = start_pj; pl < end_pj; ++pl) {
                    p_jkl = &(thb_intrs->select.three_body_list[pl]);
                    l = p_jkl->thb;
                    plk = p_jkl->pthb; //pointer to l on k's bond_list!
                    pbond_kl = &(bonds->select.bond_list[plk]);
                    bo_kl = &(pbond_kl->bo_data);
                    type_l = system->my_atoms[l].type;
                    fbh = &(system->reax_param.fbp[type_i][type_j]
                            [type_k][type_l]);
                    fbp = &(system->reax_param.fbp[type_i][type_j]
                            [type_k][type_l].prm[0]);

                    if (i != l && fbh->cnt &&
                        bo_kl->BO > control->thb_cut/*0*/ &&
                        bo_ij->BO * bo_jk->BO * bo_kl->BO > control->thb_cut/*0*/) {
                      ++num_frb_intrs;
                      //fprintf(stderr,
                      //      "%5d: %6d %6d %6d %6d\n", num_frb_intrs,
                      //      system->my_atoms[i].orig_id,system->my_atoms[j].orig_id,
                      //      system->my_atoms[k].orig_id,system->my_atoms[l].orig_id);

                      r_kl = pbond_kl->d;
                      BOA_kl = bo_kl->BO - control->thb_cut;

                      theta_jkl = p_jkl->theta;
                      sin_jkl = sin(theta_jkl);
                      cos_jkl = cos(theta_jkl);
                      //tan_jkl_i = 1. / tan(theta_jkl);
                      if (sin_jkl >= 0 && sin_jkl <= MIN_SINE)
                        tan_jkl_i = cos_jkl / MIN_SINE;
                      else if (sin_jkl <= 0 && sin_jkl >= -MIN_SINE)
                        tan_jkl_i = cos_jkl / -MIN_SINE;
                      else tan_jkl_i = cos_jkl /sin_jkl;

                      rvec_ScaledSum(dvec_li, 1., system->my_atoms[i].x,
                                     -1., system->my_atoms[l].x);
                      r_li = rvec_Norm(dvec_li);


                      /* omega and its derivative */
                      omega = Calculate_Omega(pbond_ij->dvec, r_ij,
                                              pbond_jk->dvec, r_jk,
                                              pbond_kl->dvec, r_kl,
                                              dvec_li, r_li,
                                              p_ijk, p_jkl,
                                              dcos_omega_di, dcos_omega_dj,
                                              dcos_omega_dk, dcos_omega_dl);

                      cos_omega = cos(omega);
                      cos2omega = cos(2. * omega);
                      cos3omega = cos(3. * omega);
                      /* end omega calculations */

                      /* torsion energy */
                      exp_tor1 = exp(fbp->p_tor1 *
                                     SQR(2.0 - bo_jk->BO_pi - f11_DjDk));
                      exp_tor2_kl = exp(-p_tor2 * BOA_kl);
                      exp_cot2_kl = exp(-p_cot2 * SQR(BOA_kl - 1.5));
                      fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) *
                        (1.0 - exp_tor2_kl);

                      CV = 0.5 * (fbp->V1 * (1.0 + cos_omega) +
                                  fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
                                  fbp->V3 * (1.0 + cos3omega));

                      total_Etor += e_tor = fn10 * sin_ijk * sin_jkl * CV;

                      dfn11 = (-p_tor3 * exp_tor3_DjDk +
                               (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
                               (2.0 + exp_tor3_DjDk) * exp_tor34_inv) *
                        exp_tor34_inv;

                      CEtors1 = sin_ijk * sin_jkl * CV;

                      CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                        (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0 - SQR(cos_omega)) *
                        sin_ijk * sin_jkl;
                      CEtors3 = CEtors2 * dfn11;

                      CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                        (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
                      CEtors5 = CEtors1 * p_tor2 *
                        (1.0 - exp_tor2_ij) * exp_tor2_jk * (1.0 - exp_tor2_kl);
                      CEtors6 = CEtors1 * p_tor2 *
                        (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * exp_tor2_kl;

                      cmn = -fn10 * CV;
                      CEtors7 = cmn * sin_jkl * tan_ijk_i;
                      CEtors8 = cmn * sin_ijk * tan_jkl_i;

                      CEtors9 = fn10 * sin_ijk * sin_jkl *
                        (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
                         1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
                      /* end  of torsion energy */


                      /* 4-body conjugation energy */
                      fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
                      //data->my_en.e_con += e_con =
                      total_Econ += e_con =
                        fbp->p_cot1 * fn12 *
                        (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                      Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                        (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                      CEconj1 = Cconj * (BOA_ij - 1.5e0);
                      CEconj2 = Cconj * (BOA_jk - 1.5e0);
                      CEconj3 = Cconj * (BOA_kl - 1.5e0);

                      CEconj4 = -fbp->p_cot1 * fn12 *
                        (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
                      CEconj5 = -fbp->p_cot1 * fn12 *
                        (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
                      CEconj6 = 2.0 * fbp->p_cot1 * fn12 *
                        cos_omega * sin_ijk * sin_jkl;
                      /* end 4-body conjugation energy */

                      /* FORCES */
                      bo_jk->Cdbopi += CEtors2;
                      workspace->CdDelta[j] += CEtors3;
                      //workspace->CdDelta[k] += CEtors3;
                      workspace->CdDeltaReduction[reductionOffset+k] += CEtors3;
                      bo_ij->Cdbo += (CEtors4 + CEconj1);
                      bo_jk->Cdbo += (CEtors5 + CEconj2);
                      //bo_kl->Cdbo += (CEtors6 + CEconj3);
                      bo_kl->CdboReduction[tid] += (CEtors6 + CEconj3);

                      /* dcos_theta_ijk */
                      rvec_ScaledAdd(workspace->f[j],
                                     CEtors7 + CEconj4, p_ijk->dcos_dj);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+i],
                                     CEtors7 + CEconj4, p_ijk->dcos_dk);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k],
                                     CEtors7 + CEconj4, p_ijk->dcos_di);

                      /* dcos_theta_jkl */
                      rvec_ScaledAdd(workspace->f[j],
                                     CEtors8 + CEconj5, p_jkl->dcos_di);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k],
                                     CEtors8 + CEconj5, p_jkl->dcos_dj);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+l],
                                     CEtors8 + CEconj5, p_jkl->dcos_dk);

                      /* dcos_omega */
                      rvec_ScaledAdd(workspace->f[j],
                                     CEtors9 + CEconj6, dcos_omega_dj);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+i],
                                     CEtors9 + CEconj6, dcos_omega_di);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k],
                                     CEtors9 + CEconj6, dcos_omega_dk);
                      rvec_ScaledAdd(workspace->forceReduction[reductionOffset+l],
                                     CEtors9 + CEconj6, dcos_omega_dl);

                      /* tally into per-atom virials */
                      if (system->pair_ptr->evflag) {

                        // acquire vectors
                        rvec_ScaledSum(delil, 1., system->my_atoms[l].x,
                                       -1., system->my_atoms[i].x);
                        rvec_ScaledSum(deljl, 1., system->my_atoms[l].x,
                                       -1., system->my_atoms[j].x);
                        rvec_ScaledSum(delkl, 1., system->my_atoms[l].x,
                                       -1., system->my_atoms[k].x);
                        // dcos_theta_ijk
                        rvec_Scale(fi_tmp, CEtors7 + CEconj4, p_ijk->dcos_dk);
                        rvec_Scale(fj_tmp, CEtors7 + CEconj4, p_ijk->dcos_dj);
                        rvec_Scale(fk_tmp, CEtors7 + CEconj4, p_ijk->dcos_di);

                        // dcos_theta_jkl
                        rvec_ScaledAdd(fj_tmp, CEtors8 + CEconj5, p_jkl->dcos_di);
                        rvec_ScaledAdd(fk_tmp, CEtors8 + CEconj5, p_jkl->dcos_dj);

                        // dcos_omega
                        rvec_ScaledAdd(fi_tmp, CEtors9 + CEconj6, dcos_omega_di);
                        rvec_ScaledAdd(fj_tmp, CEtors9 + CEconj6, dcos_omega_dj);
                        rvec_ScaledAdd(fk_tmp, CEtors9 + CEconj6, dcos_omega_dk);

                        // tally
                        eng_tmp = e_tor + e_con;

                        if (system->pair_ptr->eflag_either)
                          pair_reax_ptr->ev_tally_thr_proxy( j, k, system->n, 1,
                                                            eng_tmp, 0.0, 0.0, 0.0, 0.0, 0.0, thr);

                        if (system->pair_ptr->vflag_either)
                          pair_reax_ptr->v_tally4_thr_proxy(i, j, k, l, fi_tmp, fj_tmp, fk_tmp,
                                                               delil, deljl, delkl, thr);
                      }

                    } // pl check ends
                  } // pl loop ends
                } // pi check ends
              } // pi loop ends
            } // k-j neighbor check ends
          } // j<k && j-k neighbor check ends
        } // pk loop ends
      } // j loop

    } // end omp parallel

    data->my_en.e_tor = total_Etor;
    data->my_en.e_con = total_Econ;
  }
}
