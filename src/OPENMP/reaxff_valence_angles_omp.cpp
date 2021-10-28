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

#include "error.h"
#include "fix_omp.h"
#include "pair_reaxff_omp.h"
#include "reaxff_api.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

namespace ReaxFF {
  void Calculate_dCos_ThetaOMP(rvec dvec_ji, double d_ji, rvec dvec_jk,
                               double d_jk, rvec *dcos_theta_di,
                               rvec *dcos_theta_dj, rvec *dcos_theta_dk)
  {
    double sqr_d_ji = SQR(d_ji);
    double sqr_d_jk = SQR(d_jk);
    double inv_dists = 1.0 / (d_ji * d_jk);
    double inv_dists3 = inv_dists * inv_dists * inv_dists;
    double dot_dvecs = dvec_ji[0]*dvec_jk[0] + dvec_ji[1]*dvec_jk[1] + dvec_ji[2]*dvec_jk[2];
    double Cdot_inv3 = dot_dvecs * inv_dists3;

    double csqr_jk = Cdot_inv3 * sqr_d_jk;
    double csqr_ji = Cdot_inv3 * sqr_d_ji;

    // Try to help compiler out by unrolling
    // x-component
    double dinv_jk = dvec_jk[0] * inv_dists;
    double dinv_ji = dvec_ji[0] * inv_dists;

    double cdev_ji = csqr_jk * dvec_ji[0];
    double cdev_jk = csqr_ji * dvec_jk[0];

    (*dcos_theta_di)[0] =   dinv_jk            - cdev_ji;
    (*dcos_theta_dj)[0] = -(dinv_jk + dinv_ji) + cdev_ji + cdev_jk;
    (*dcos_theta_dk)[0] =             dinv_ji            - cdev_jk;

    // y-component
    dinv_jk = dvec_jk[1] * inv_dists;
    dinv_ji = dvec_ji[1] * inv_dists;

    cdev_ji = csqr_jk * dvec_ji[1];
    cdev_jk = csqr_ji * dvec_jk[1];

    (*dcos_theta_di)[1] =   dinv_jk            - cdev_ji;
    (*dcos_theta_dj)[1] = -(dinv_jk + dinv_ji) + cdev_ji + cdev_jk;
    (*dcos_theta_dk)[1] =             dinv_ji            - cdev_jk;

    // z-component
    dinv_jk = dvec_jk[2] * inv_dists;
    dinv_ji = dvec_ji[2] * inv_dists;

    cdev_ji = csqr_jk * dvec_ji[2];
    cdev_jk = csqr_ji * dvec_jk[2];

    (*dcos_theta_di)[2] =   dinv_jk            - cdev_ji;
    (*dcos_theta_dj)[2] = -(dinv_jk + dinv_ji) + cdev_ji + cdev_jk;
    (*dcos_theta_dk)[2] =             dinv_ji            - cdev_jk;
  }

  /* ---------------------------------------------------------------------- */

  /* this is a 3-body interaction in which the main role is
     played by j which sits in the middle of the other two. */
  void Valence_AnglesOMP(reax_system *system, control_params *control,
                         simulation_data *data, storage *workspace,
                         reax_list **lists)
  {
    reax_list *bonds = (*lists) + BONDS;
    reax_list *thb_intrs =  (*lists) + THREE_BODIES;

    // Precompute and store valence_angle offsets for OpenMP code.
    int * _my_offset = workspace->valence_angle_atom_myoffset;

    /* global parameters used in these calculations */
    double p_val6 = system->reax_param.gp.l[14];
    double p_val8 = system->reax_param.gp.l[33];
    double p_val9 = system->reax_param.gp.l[16];
    double p_val10 = system->reax_param.gp.l[17];
    double total_Eang = 0;
    double total_Epen = 0;
    double total_Ecoa = 0;

    int  nthreads = control->nthreads;
    int  num_thb_intrs = 0;
    int  TWICE = 2;
#if defined(_OPENMP)
#pragma omp parallel default(shared) reduction(+:total_Eang, total_Epen, total_Ecoa, num_thb_intrs)
#endif
    {
      int i, j, pi, k, pk, t;
      int type_i, type_j, type_k;
      int start_j, end_j, start_pk, end_pk;
      int cnt, my_offset;

      double temp, temp_bo_jt, pBOjt7;
      double p_val1, p_val2, p_val3, p_val4, p_val5, p_val7;
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

      // Tallying variables
      double eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
      double delij[3], delkj[3];

      three_body_header *thbh;
      three_body_parameters *thbp;
      three_body_interaction_data *p_ijk, *p_kji;
      bond_data *pbond_ij, *pbond_jk, *pbond_jt;
      bond_order_data *bo_ij, *bo_jk, *bo_jt;

      int tid = get_tid();

      long reductionOffset = (system->N * tid);
      class PairReaxFFOMP *pair_reax_ptr;
      pair_reax_ptr = static_cast<class PairReaxFFOMP*>(system->pair_ptr);
      class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

      // Run through a minimal for (j<N) loop once to precompute offsets with safe number of threads

      const int per_thread = thb_intrs->num_intrs / nthreads;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
      for (j = 0; j < system->N; ++j) {
        type_j = system->my_atoms[j].type;
        _my_offset[j] = 0;
        if (type_j < 0) continue;

        start_j = Start_Index(j, bonds);
        end_j = End_Index(j, bonds);

        // Always point to start of workspace to count angles
        my_offset = tid * per_thread;

        for (pi = start_j; pi < end_j; ++pi) {
          Set_Start_Index(pi, my_offset, thb_intrs);
          pbond_ij = &(bonds->select.bond_list[pi]);
          bo_ij = &(pbond_ij->bo_data);
          BOA_ij = bo_ij->BO - control->thb_cut;

          if (BOA_ij > 0.0) {
            i = pbond_ij->nbr;

            /* first copy 3-body intrs from previously computed ones where i>k.
               in the second for-loop below,
               we compute only new 3-body intrs where i < k */
            for (pk = start_j; pk < pi; ++pk) {
              start_pk = Start_Index(pk, thb_intrs);
              end_pk = End_Index(pk, thb_intrs);

              for (t = start_pk; t < end_pk; ++t)
                if (thb_intrs->select.three_body_list[t].thb == i) {

                  p_ijk = &(thb_intrs->select.three_body_list[my_offset]);
                  p_ijk->thb = bonds->select.bond_list[pk].nbr;

                  ++my_offset;
                  break;
                }
            } // for (pk)

            /* and this is the second for loop mentioned above */
            for (pk = pi+1; pk < end_j; ++pk) {
              pbond_jk = &(bonds->select.bond_list[pk]);
              k        = pbond_jk->nbr;

              if (j >= system->n && i >= system->n && k >= system->n) continue;

              p_ijk    = &(thb_intrs->select.three_body_list[my_offset]);
              p_ijk->thb = k;

              ++my_offset; // add this  to the list of 3-body interactions
            } // for (pk)
          } // if ()

          Set_End_Index(pi, my_offset, thb_intrs);
        } // for (pi)

        // Confirm that thb_intrs->num_intrs / nthreads is enough to hold all angles from a single atom
        if (my_offset >= (tid+1)*per_thread)
          control->error_ptr->one(FLERR, fmt::format("step {}: ran out of space on "
                                                     "angle_list for atom {}:\n"
                                                     " nthreads={} tid={} my_offset={} per_thread={}\n"
                                                     " num_intrs={} N={}",data->step,j,nthreads,tid,
                                                     my_offset,per_thread,thb_intrs->num_intrs,system->N));
        // Number of angles owned by this atom
        _my_offset[j] = my_offset - tid * per_thread;
      } // for (j)

      // Wait for all threads to finish counting angles
#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp barrier
#endif
      // Master thread uses angle counts to compute offsets
      // This can be threaded
#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp master
#endif
      {
        int current_count = 0;
        int m = _my_offset[0];
        _my_offset[0] = current_count;
        for (j=1; j<system->N; j++) {
          current_count+= m;
          m = _my_offset[j];
          _my_offset[j] = current_count;
        }
        _my_offset[system->N] = current_count + m; // Used to test if last particle has any angles
      }

      // All threads wait till master thread finished computing offsets
#if defined(_OPENMP) && !defined(__NVCC__)
#pragma omp barrier
#endif
      // Original loop, but now using precomputed offsets
      // Safe to use all threads available, regardless of threads tasked above
      // We also now skip over atoms that have no angles assigned
#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)//(dynamic,chunksize)//(guided)
#endif
      for (j = 0; j < system->N; ++j) {         // Ray: the first one with system->N
        type_j = system->my_atoms[j].type;
        if (type_j < 0) continue;

        // Skip if no angles for this atom
        if (_my_offset[j] == _my_offset[j+1]) continue;

        start_j = Start_Index(j, bonds);
        end_j = End_Index(j, bonds);

        type_j = system->my_atoms[j].type;

        my_offset = _my_offset[j];

        p_val3 = system->reax_param.sbp[type_j].p_val3;
        p_val5 = system->reax_param.sbp[type_j].p_val5;

        SBOp = 0, prod_SBO = 1;
        for (t = start_j; t < end_j; ++t) {
          bo_jt = &(bonds->select.bond_list[t].bo_data);
          SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
          temp = SQR(bo_jt->BO);
          temp *= temp;
          temp *= temp;
          prod_SBO *= exp(-temp);
        }

        // modifications to match Adri's code - 09/01/09
        if (workspace->vlpex[j] >= 0) {
          vlpadj = 0;
          dSBO2 = prod_SBO - 1;
        } else {
          vlpadj = workspace->nlp[j];
          dSBO2 = (prod_SBO - 1) * (1 - p_val8 * workspace->dDelta_lp[j]);
        }

        SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
        dSBO1 = -8 * prod_SBO * (workspace->Delta_boc[j] + p_val8 * vlpadj);

        if (SBO <= 0)
          SBO2 = 0, CSBO2 = 0;
        else if (SBO > 0 && SBO <= 1) {
          SBO2 = pow(SBO, p_val9);
          CSBO2 = p_val9 * pow(SBO, p_val9 - 1);
        }
        else if (SBO > 1 && SBO < 2) {
          SBO2 = 2 - pow(2-SBO, p_val9);
          CSBO2 = p_val9 * pow(2 - SBO, p_val9 - 1);
        }
        else
          SBO2 = 2, CSBO2 = 0;

        expval6 = exp(p_val6 * workspace->Delta_boc[j]);

        for (pi = start_j; pi < end_j; ++pi) {
          Set_Start_Index(pi, my_offset, thb_intrs);
          pbond_ij = &(bonds->select.bond_list[pi]);
          bo_ij = &(pbond_ij->bo_data);
          BOA_ij = bo_ij->BO - control->thb_cut;


          if (BOA_ij > 0.0) {
            i = pbond_ij->nbr;
            type_i = system->my_atoms[i].type;


            /* first copy 3-body intrs from previously computed ones where i>k.
               in the second for-loop below,
               we compute only new 3-body intrs where i < k */
            for (pk = start_j; pk < pi; ++pk) {
              start_pk = Start_Index(pk, thb_intrs);
              end_pk = End_Index(pk, thb_intrs);

              for (t = start_pk; t < end_pk; ++t)
                if (thb_intrs->select.three_body_list[t].thb == i) {
                  p_ijk = &(thb_intrs->select.three_body_list[my_offset]);
                  p_kji = &(thb_intrs->select.three_body_list[t]);

                  p_ijk->thb = bonds->select.bond_list[pk].nbr;
                  p_ijk->pthb  = pk;
                  p_ijk->theta = p_kji->theta;
                  rvec_Copy(p_ijk->dcos_di, p_kji->dcos_dk);
                  rvec_Copy(p_ijk->dcos_dj, p_kji->dcos_dj);
                  rvec_Copy(p_ijk->dcos_dk, p_kji->dcos_di);

                  ++my_offset;
                  ++num_thb_intrs;
                  break;
                }
            } // for (pk)


            /* and this is the second for loop mentioned above */
            for (pk = pi+1; pk < end_j; ++pk) {
              pbond_jk = &(bonds->select.bond_list[pk]);
              bo_jk    = &(pbond_jk->bo_data);
              BOA_jk   = bo_jk->BO - control->thb_cut;
              k        = pbond_jk->nbr;
              type_k   = system->my_atoms[k].type;
              p_ijk    = &(thb_intrs->select.three_body_list[my_offset]);

              // Fix by Sudhir
              // if (BOA_jk <= 0) continue;
              if (j >= system->n && i >= system->n && k >= system->n) continue;

              Calculate_Theta(pbond_ij->dvec, pbond_ij->d,
                              pbond_jk->dvec, pbond_jk->d,
                              &theta, &cos_theta);

              Calculate_dCos_ThetaOMP(pbond_ij->dvec, pbond_ij->d,
                                      pbond_jk->dvec, pbond_jk->d,
                                      &(p_ijk->dcos_di), &(p_ijk->dcos_dj),
                                      &(p_ijk->dcos_dk));
              p_ijk->thb = k;
              p_ijk->pthb = pk;
              p_ijk->theta = theta;

              sin_theta = sin(theta);
              if (sin_theta < 1.0e-5)
                sin_theta = 1.0e-5;

              ++my_offset; // add this  to the list of 3-body interactions
              ++num_thb_intrs;

              if ((j < system->n) && (BOA_jk > 0.0) &&
                  (bo_ij->BO > control->thb_cut) &&
                  (bo_jk->BO > control->thb_cut) &&
                  (bo_ij->BO * bo_jk->BO > control->thb_cutsq)) {
                thbh = &(system->reax_param.thbp[type_i][type_j][type_k]);

                for (cnt = 0; cnt < thbh->cnt; ++cnt) {

                  if (fabs(thbh->prm[cnt].p_val1) > 0.001) {
                    thbp = &(thbh->prm[cnt]);

                    /* ANGLE ENERGY */
                    p_val1 = thbp->p_val1;
                    p_val2 = thbp->p_val2;
                    p_val4 = thbp->p_val4;
                    p_val7 = thbp->p_val7;
                    theta_00 = thbp->theta_00;

                    exp3ij = exp(-p_val3 * pow(BOA_ij, p_val4));
                    f7_ij = 1.0 - exp3ij;
                    Cf7ij = p_val3 * p_val4 * pow(BOA_ij, p_val4 - 1.0) * exp3ij;

                    exp3jk = exp(-p_val3 * pow(BOA_jk, p_val4));
                    f7_jk = 1.0 - exp3jk;
                    Cf7jk = p_val3 * p_val4 * pow(BOA_jk, p_val4 - 1.0) * exp3jk;

                    expval7 = exp(-p_val7 * workspace->Delta_boc[j]);
                    trm8 = 1.0 + expval6 + expval7;
                    f8_Dj = p_val5 - ((p_val5 - 1.0) * (2.0 + expval6) / trm8);
                    Cf8j = ((1.0 - p_val5) / SQR(trm8)) *
                      (p_val6 * expval6 * trm8 -
                       (2.0 + expval6) * (p_val6*expval6 - p_val7*expval7));

                    theta_0 = 180.0 - theta_00 * (1.0 -
                                                  exp(-p_val10 * (2.0 - SBO2)));
                    theta_0 = DEG2RAD(theta_0);

                    expval2theta  = exp(-p_val2 * SQR(theta_0 - theta));
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
                      exp(-p_val10 * (2.0 - SBO2));

                    CEval5 = -CEval4 * Ctheta_0 * CSBO2;
                    CEval6 = CEval5 * dSBO1;
                    CEval7 = CEval5 * dSBO2;
                    CEval8 = -CEval4 / sin_theta;

                    total_Eang += e_ang =
                      f7_ij * f7_jk * f8_Dj * expval12theta;
                    /* END ANGLE ENERGY*/


                    /* PENALTY ENERGY */
                    p_pen1 = thbp->p_pen1;
                    p_pen2 = system->reax_param.gp.l[19];
                    p_pen3 = system->reax_param.gp.l[20];
                    p_pen4 = system->reax_param.gp.l[21];

                    exp_pen2ij = exp(-p_pen2 * SQR(BOA_ij - 2.0));
                    exp_pen2jk = exp(-p_pen2 * SQR(BOA_jk - 2.0));
                    exp_pen3 = exp(-p_pen3 * workspace->Delta[j]);
                    exp_pen4 = exp(p_pen4 * workspace->Delta[j]);
                    trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
                    f9_Dj = (2.0 + exp_pen3) / trm_pen34;
                    Cf9j = (-p_pen3 * exp_pen3 * trm_pen34 -
                            (2.0 + exp_pen3) * (-p_pen3 * exp_pen3 +
                                                p_pen4 * exp_pen4)) /
                      SQR(trm_pen34);

                    total_Epen += e_pen =
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

                    exp_coa2 = exp(p_coa2 * workspace->Delta_val[j]);
                    total_Ecoa += e_coa =
                      p_coa1 / (1. + exp_coa2) *
                      exp(-p_coa3 * SQR(workspace->total_bond_order[i]-BOA_ij)) *
                      exp(-p_coa3 * SQR(workspace->total_bond_order[k]-BOA_jk)) *
                      exp(-p_coa4 * SQR(BOA_ij - 1.5)) *
                      exp(-p_coa4 * SQR(BOA_jk - 1.5));

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
                    workspace->CdDeltaReduction[reductionOffset+i] += CEcoa4;
                    workspace->CdDeltaReduction[reductionOffset+k] += CEcoa5;

                    for (t = start_j; t < end_j; ++t) {
                      pbond_jt = &(bonds->select.bond_list[t]);
                      bo_jt = &(pbond_jt->bo_data);
                      temp_bo_jt = bo_jt->BO;
                      temp = CUBE(temp_bo_jt);
                      pBOjt7 = temp * temp * temp_bo_jt;

                      bo_jt->Cdbo += (CEval6 * pBOjt7);
                      bo_jt->Cdbopi += CEval5;
                      bo_jt->Cdbopi2 += CEval5;
                    }

                    rvec_ScaledAdd(workspace->f[j], CEval8, p_ijk->dcos_dj);
                    rvec_ScaledAdd(workspace->forceReduction[reductionOffset+i], CEval8, p_ijk->dcos_di);
                    rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k], CEval8, p_ijk->dcos_dk);

                    /* tally into per-atom virials */
                    if (system->pair_ptr->evflag) {

                      /* Acquire vectors */
                      rvec_ScaledSum(delij, 1., system->my_atoms[i].x,
                                     -1., system->my_atoms[j].x);
                      rvec_ScaledSum(delkj, 1., system->my_atoms[k].x,
                                     -1., system->my_atoms[j].x);

                      rvec_Scale(fi_tmp, -CEval8, p_ijk->dcos_di);
                      rvec_Scale(fj_tmp, -CEval8, p_ijk->dcos_dj);
                      rvec_Scale(fk_tmp, -CEval8, p_ijk->dcos_dk);

                      eng_tmp = e_ang + e_pen + e_coa;

                      if (system->pair_ptr->eflag_either)
                        pair_reax_ptr->ev_tally_thr_proxy( j, j, system->N, 1,
                                                          eng_tmp, 0.0, 0.0, 0.0, 0.0, 0.0, thr);
                      if (system->pair_ptr->vflag_either)
                        pair_reax_ptr->v_tally3_thr_proxy(i, j, k, fi_tmp, fk_tmp, delij, delkj, thr);
                    }

                  } // if (p_val1>0.001)
                } // for (cnt)
              } // if (j<n && BOA_jk>0)
            } // for (pk)
          } // if (BOA_ij>0)

          Set_End_Index(pi, my_offset, thb_intrs);
        } // for (pi)
      } // for (j)
    } // end omp parallel

    data->my_en.e_ang = total_Eang;
    data->my_en.e_pen = total_Epen;
    data->my_en.e_coa = total_Ecoa;

    if (num_thb_intrs >= thb_intrs->num_intrs * DANGER_ZONE) {
      workspace->realloc.num_3body = num_thb_intrs * TWICE;
      if (num_thb_intrs > thb_intrs->num_intrs)
        control->error_ptr->one(FLERR, fmt::format("step {}: ran out of space on "
                                                   "angle_list: top={}, max={}",
                                                   data->step, num_thb_intrs,
                                                   thb_intrs->num_intrs));
    }
  }
}
