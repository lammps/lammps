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

#include "pair_reaxff_omp.h"
#include "reaxff_api.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

namespace ReaxFF {
  void vdW_Coulomb_Energy_OMP(reax_system *system, control_params *control,
                              simulation_data *data, storage *workspace,
                              reax_list **lists)
  {
    int natoms = system->n;
    reax_list *far_nbrs = (*lists) + FAR_NBRS;
    double p_vdW1 = system->reax_param.gp.l[28];
    double p_vdW1i = 1.0 / p_vdW1;
    double total_EvdW = 0.;
    double total_Eele = 0.;

#if defined(_OPENMP)
#pragma omp parallel default(shared) reduction(+: total_EvdW, total_Eele)
#endif
    {
      int tid = get_tid();
      int i, j, pj;
      int start_i, end_i, orig_i, orig_j, flag;
      double powr_vdW1, powgi_vdW1;
      double tmp, r_ij, fn13, exp1, exp2;
      double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
      double dr3gamij_1, dr3gamij_3;
      double e_ele, e_vdW, e_core;
      const double SMALL = 0.0001;
      double e_lg, de_lg, r_ij5, r_ij6, re6;
      two_body_parameters *twbp;
      far_neighbor_data *nbr_pj;

      // Tallying variables:
      double pe_vdw, f_tmp, delij[3];

      long reductionOffset = (system->N * tid);

      class PairReaxFFOMP *pair_reax_ptr;
      pair_reax_ptr = static_cast<class PairReaxFFOMP*>(system->pair_ptr);
      class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

      e_core = 0;
      e_vdW = 0;
      e_lg = 0;
      de_lg = 0.0;

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (i = 0; i < natoms; ++i) {
        if (system->my_atoms[i].type < 0) continue;
        start_i = Start_Index(i, far_nbrs);
        end_i   = End_Index(i, far_nbrs);
        orig_i  = system->my_atoms[i].orig_id;

        for (pj = start_i; pj < end_i; ++pj) {
          nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
          j = nbr_pj->nbr;
          orig_j  = system->my_atoms[j].orig_id;

          flag = 0;
          if (nbr_pj->d <= control->nonb_cut) {
            if (j < natoms) flag = 1;
            else if (orig_i < orig_j) flag = 1;
            else if (orig_i == orig_j) {
              if (nbr_pj->dvec[2] > SMALL) flag = 1;
              else if (fabs(nbr_pj->dvec[2]) < SMALL) {
                if (nbr_pj->dvec[1] > SMALL) flag = 1;
                else if (fabs(nbr_pj->dvec[1]) < SMALL && nbr_pj->dvec[0] > SMALL)
                  flag = 1;
              }
            }
          }

          if (flag) {

            r_ij = nbr_pj->d;
            twbp = &(system->reax_param.tbp[system->my_atoms[i].type]
                     [system->my_atoms[j].type]);

            /* Calculate Taper and its derivative */
            // Tap = nbr_pj->Tap;   -- precomputed during compte_H
            Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
            Tap = Tap * r_ij + workspace->Tap[5];
            Tap = Tap * r_ij + workspace->Tap[4];
            Tap = Tap * r_ij + workspace->Tap[3];
            Tap = Tap * r_ij + workspace->Tap[2];
            Tap = Tap * r_ij + workspace->Tap[1];
            Tap = Tap * r_ij + workspace->Tap[0];

            dTap = 7*workspace->Tap[7] * r_ij + 6*workspace->Tap[6];
            dTap = dTap * r_ij + 5*workspace->Tap[5];
            dTap = dTap * r_ij + 4*workspace->Tap[4];
            dTap = dTap * r_ij + 3*workspace->Tap[3];
            dTap = dTap * r_ij + 2*workspace->Tap[2];
            dTap += workspace->Tap[1]/r_ij;

            /*vdWaals Calculations*/
            if (system->reax_param.gp.vdw_type==1 || system->reax_param.gp.vdw_type==3)
              { // shielding
                powr_vdW1 = pow(r_ij, p_vdW1);
                powgi_vdW1 = pow(1.0 / twbp->gamma_w, p_vdW1);

                fn13 = pow(powr_vdW1 + powgi_vdW1, p_vdW1i);
                exp1 = exp(twbp->alpha * (1.0 - fn13 / twbp->r_vdW));
                exp2 = exp(0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW));

                e_vdW = twbp->D * (exp1 - 2.0 * exp2);
                total_EvdW += Tap * e_vdW;

                dfn13 = pow(powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0) *
                  pow(r_ij, p_vdW1 - 2.0);

                CEvd = dTap * e_vdW -
                  Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
              }
            else { // no shielding
              exp1 = exp(twbp->alpha * (1.0 - r_ij / twbp->r_vdW));
              exp2 = exp(0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW));

              e_vdW = twbp->D * (exp1 - 2.0 * exp2);
              total_EvdW += Tap * e_vdW;

              CEvd = dTap * e_vdW -
                Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) / r_ij;
            }

            if (system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
              { // innner wall
                e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij/twbp->rcore)));
                total_EvdW += Tap * e_core;

                de_core = -(twbp->acore/twbp->rcore) * e_core;
                CEvd += dTap * e_core + Tap * de_core / r_ij;

                //  lg correction, only if lgvdw is yes
                if (control->lgflag) {
                  r_ij5 = pow(r_ij, 5.0);
                  r_ij6 = pow(r_ij, 6.0);
                  re6 = pow(twbp->lgre, 6.0);

                  e_lg = -(twbp->lgcij/(r_ij6 + re6));
                  total_EvdW += Tap * e_lg;

                  de_lg = -6.0 * e_lg *  r_ij5 / (r_ij6 + re6) ;
                  CEvd += dTap * e_lg + Tap * de_lg / r_ij;
                }

              }

            /*Coulomb Calculations*/
            dr3gamij_1 = (r_ij * r_ij * r_ij + twbp->gamma);
            dr3gamij_3 = pow(dr3gamij_1 , 0.33333333333333);

            tmp = Tap / dr3gamij_3;
            total_Eele += e_ele =
              C_ele * system->my_atoms[i].q * system->my_atoms[j].q * tmp;

            CEclmb = C_ele * system->my_atoms[i].q * system->my_atoms[j].q *
              (dTap -  Tap * r_ij / dr3gamij_1) / dr3gamij_3;

            /* tally into per-atom energy */
            if (system->pair_ptr->evflag || system->pair_ptr->vflag_atom) {
              pe_vdw = Tap * (e_vdW + e_core + e_lg);
              rvec_ScaledSum(delij, 1., system->my_atoms[i].x,
                              -1., system->my_atoms[j].x);
              f_tmp = -(CEvd + CEclmb);
              pair_reax_ptr->ev_tally_thr_proxy( i, j, natoms,
                                                1, pe_vdw, e_ele, f_tmp,
                                                delij[0], delij[1], delij[2], thr);
            }

            rvec_ScaledAdd(workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec);
            rvec_ScaledAdd(workspace->forceReduction[reductionOffset+j],
                           +(CEvd + CEclmb), nbr_pj->dvec);
          }
        }
      }

      pair_reax_ptr->reduce_thr_proxy(system->pair_ptr, system->pair_ptr->eflag_either,
                                      system->pair_ptr->vflag_either, thr);
    } // parallel region

    data->my_en.e_vdW = total_EvdW;
    data->my_en.e_ele = total_Eele;

    Compute_Polarization_Energy(system, data, workspace);
  }

/* ---------------------------------------------------------------------- */

  void Tabulated_vdW_Coulomb_Energy_OMP(reax_system *system,control_params *control,
                                        simulation_data *data, storage *workspace,
                                        reax_list **lists) {

    double SMALL = 0.0001;
    int  natoms = system->n;
    reax_list *far_nbrs = (*lists) + FAR_NBRS;
    double total_EvdW = 0.;
    double total_Eele = 0.;

#if defined(_OPENMP)
#pragma omp parallel default(shared) reduction(+:total_EvdW, total_Eele)
#endif
    {
      int i, j, pj, r;
      int type_i, type_j, tmin, tmax;
      int start_i, end_i, orig_i, orig_j, flag;
      double r_ij, base, dif;
      double e_vdW, e_ele;
      double CEvd, CEclmb;
      double f_tmp, delij[3];
      far_neighbor_data *nbr_pj;
      LR_lookup_table *t;

      int tid = get_tid();
      long froffset = (system->N * tid);
      LR_lookup_table ** & LR = system->LR;

      class PairReaxFFOMP *pair_reax_ptr;
      pair_reax_ptr = static_cast<class PairReaxFFOMP*>(system->pair_ptr);
      class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (i = 0; i < natoms; ++i) {
        type_i  = system->my_atoms[i].type;
        if (type_i < 0) continue;
        start_i = Start_Index(i,far_nbrs);
        end_i   = End_Index(i,far_nbrs);
        orig_i  = system->my_atoms[i].orig_id;

        for (pj = start_i; pj < end_i; ++pj) {
          nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
          j = nbr_pj->nbr;
          type_j = system->my_atoms[j].type;
          if (type_j < 0) continue;
          orig_j  = system->my_atoms[j].orig_id;

          flag = 0;
          if (nbr_pj->d <= control->nonb_cut) {
            if (j < natoms) flag = 1;
            else if (orig_i < orig_j) flag = 1;
            else if (orig_i == orig_j) {
              if (nbr_pj->dvec[2] > SMALL) flag = 1;
              else if (fabs(nbr_pj->dvec[2]) < SMALL) {
                if (nbr_pj->dvec[1] > SMALL) flag = 1;
                else if (fabs(nbr_pj->dvec[1]) < SMALL && nbr_pj->dvec[0] > SMALL)
                  flag = 1;
              }
            }

          }

          if (flag) {

            r_ij   = nbr_pj->d;
            tmin  = MIN(type_i, type_j);
            tmax  = MAX(type_i, type_j);
            t = &(LR[tmin][tmax]);

            /* Cubic Spline Interpolation */
            r = (int)(r_ij * t->inv_dx);
            if (r == 0)  ++r;
            base = (double)(r+1) * t->dx;
            dif = r_ij - base;

            e_vdW = ((t->vdW[r].d*dif + t->vdW[r].c)*dif + t->vdW[r].b)*dif +
              t->vdW[r].a;

            e_ele = ((t->ele[r].d*dif + t->ele[r].c)*dif + t->ele[r].b)*dif +
              t->ele[r].a;
            e_ele *= system->my_atoms[i].q * system->my_atoms[j].q;

            total_EvdW += e_vdW;
            total_Eele += e_ele;

            CEvd = ((t->CEvd[r].d*dif + t->CEvd[r].c)*dif + t->CEvd[r].b)*dif +
              t->CEvd[r].a;

            CEclmb = ((t->CEclmb[r].d*dif+t->CEclmb[r].c)*dif+t->CEclmb[r].b)*dif +
              t->CEclmb[r].a;
            CEclmb *= system->my_atoms[i].q * system->my_atoms[j].q;

            /* tally into per-atom energy */
            if (system->pair_ptr->evflag) {
              rvec_ScaledSum(delij, 1., system->my_atoms[i].x,
                              -1., system->my_atoms[j].x);
              f_tmp = -(CEvd + CEclmb);
              pair_reax_ptr->ev_tally_thr_proxy( i, j, natoms, 1, e_vdW, e_ele,
                                                f_tmp, delij[0], delij[1], delij[2], thr);
            }

            rvec_ScaledAdd(workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec);
            rvec_ScaledAdd(workspace->forceReduction[froffset+j],
                           +(CEvd + CEclmb), nbr_pj->dvec);
          }
        }
      }

      pair_reax_ptr->reduce_thr_proxy(system->pair_ptr, system->pair_ptr->eflag_either,
                                      system->pair_ptr->vflag_either, thr);
    } // end omp parallel

    data->my_en.e_vdW = total_EvdW;
    data->my_en.e_ele = total_Eele;

    Compute_Polarization_Energy(system, data, workspace);
  }
}
