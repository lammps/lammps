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
#include "pair_reaxff_omp.h"

#include "reaxff_api.h"

#include <cmath>

using namespace LAMMPS_NS;

namespace ReaxFF {
/* ---------------------------------------------------------------------- */

  static void Compute_Bonded_ForcesOMP(reax_system *system, control_params *control,
                                       simulation_data *data, storage *workspace,
                                       reax_list **lists)
  {

    BOOMP(system, workspace, lists);
    BondsOMP(system, data, workspace, lists);
    Atom_EnergyOMP(system, data, workspace, lists);
    Valence_AnglesOMP(system, control, data, workspace, lists);
    Torsion_AnglesOMP(system, control, data, workspace, lists);
    if (control->hbond_cut > 0)
      Hydrogen_BondsOMP(system, control, data, workspace, lists);
  }

  static void Compute_NonBonded_ForcesOMP(reax_system *system, control_params *control,
                                          simulation_data *data, storage *workspace,
                                          reax_list **lists)
  {
    /* van der Waals and Coulomb interactions */

    if (control->tabulate == 0)
      vdW_Coulomb_Energy_OMP(system, control, data, workspace, lists);
    else
      Tabulated_vdW_Coulomb_Energy_OMP(system, control, data, workspace, lists);
  }

/* ---------------------------------------------------------------------- */

/* this version of Compute_Total_Force computes forces from
   coefficients accumulated by all interaction functions.
   Saves enormous time & space! */
  static void Compute_Total_ForceOMP(reax_system *system, control_params *control,
                                     storage *workspace, reax_list **lists)
  {
    int natoms = system->N;
    int nthreads = control->nthreads;
    long totalReductionSize = (bigint)system->N * nthreads;
    reax_list *bonds = (*lists) + BONDS;

#if defined(_OPENMP)
#pragma omp parallel default(shared) //LMP_DEFAULT_NONE
#endif
    {
      int i, j, k, pj, pk, start_j, end_j;

      int tid = get_tid();
      bond_order_data *bo_jk;

      class PairReaxFFOMP *pair_reax_ptr;
      pair_reax_ptr = static_cast<class PairReaxFFOMP*>(system->pair_ptr);
      class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

      pair_reax_ptr->ev_setup_thr_proxy(0, 1, natoms, system->pair_ptr->eatom,
                                        system->pair_ptr->vatom, nullptr, thr);

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (i = 0; i < system->N; ++i) {
        for (j = 0; j < nthreads; ++j)
          workspace->CdDelta[i] += workspace->CdDeltaReduction[system->N*j+i];
      }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
      for (j = 0; j < system->N; ++j) {
        start_j = Start_Index(j, bonds);
        end_j = End_Index(j, bonds);

        for (pk = start_j; pk < end_j; ++pk) {
          bo_jk = &(bonds->select.bond_list[pk].bo_data);
          for (k = 0; k < nthreads; ++k)
            bo_jk->Cdbo += bo_jk->CdboReduction[k];
        }
      }

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
      for (i = 0; i < system->N; ++i) {
        const int startj = Start_Index(i, bonds);
        const int endj  = End_Index(i, bonds);
        for (pj = startj; pj < endj; ++pj)
          if (i < bonds->select.bond_list[pj].nbr)
            Add_dBond_to_ForcesOMP(system, i, pj, workspace, lists);
      }

      pair_reax_ptr->reduce_thr_proxy(system->pair_ptr, 0, 1, thr);

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (i = 0; i < system->N; ++i) {
        for (j = 0; j < nthreads; ++j)
          rvec_Add(workspace->f[i], workspace->forceReduction[system->N*j+i]);
      }


#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (i = 0; i < totalReductionSize; i++) {
        workspace->forceReduction[i][0] = 0;
        workspace->forceReduction[i][1] = 0;
        workspace->forceReduction[i][2] = 0;
        workspace->CdDeltaReduction[i] = 0;
      }
    } // parallel region
  }

/* ---------------------------------------------------------------------- */

  static void Validate_ListsOMP(reax_system *system, reax_list **lists,
                         int step, int n, int N, int numH)
  {
    int comp, Hindex;
    reax_list *bonds, *hbonds;
    double saferzone = system->saferzone;

#if defined(_OPENMP)
#pragma omp parallel default(shared) private(comp,Hindex)
#endif
    {

      /* bond list */
      if (N > 0) {
        bonds = *lists + BONDS;

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
        for (int i = 0; i < N; ++i) {
          system->my_atoms[i].num_bonds = MAX(Num_Entries(i,bonds)*2, MIN_BONDS);

          if (i < N-1)
            comp = Start_Index(i+1, bonds);
          else comp = bonds->num_intrs;

          if (End_Index(i, bonds) > comp)
            system->error_ptr->one(FLERR, fmt::format("step {}: bondchk failed: "
                                                      "i={} end(i)={} str(i+1)={}\n",
                                                      step,i,End_Index(i,bonds),comp));
        }
      }


      /* hbonds list */
      if (numH > 0) {
        hbonds = *lists + HBONDS;

#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
        for (int i = 0; i < n; ++i) {
          Hindex = system->my_atoms[i].Hindex;
          if (Hindex > -1) {
            system->my_atoms[i].num_hbonds =
              (int)(MAX(Num_Entries(Hindex,hbonds)*saferzone,system->minhbonds));

            if (Hindex < numH-1)
              comp = Start_Index(Hindex+1, hbonds);
            else comp = hbonds->num_intrs;

            if (End_Index(Hindex, hbonds) > comp)
              system->error_ptr->one(FLERR, fmt::format("step {}: hbondchk failed: "
                                                        "H={} end(H)={} str(H+1)={}\n",
                                                        step, Hindex,End_Index(Hindex,hbonds),comp));
          }
        }
      }

    } // omp parallel
  }


  void Init_Forces_noQEq_OMP(reax_system *system, control_params *control,
                              simulation_data *data, storage *workspace,
                              reax_list **lists) {
    int j, pj;
    int start_i, end_i;
    int type_i, type_j;
    int ihb, jhb, ihb_top, jhb_top;
    double cutoff;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data *nbr_pj;
    reax_atom *atom_i, *atom_j;
    reax_list *far_nbrs = *lists + FAR_NBRS;
    reax_list *bonds = *lists + BONDS;
    reax_list *hbonds = *lists + HBONDS;
    int num_bonds = 0;
    int num_hbonds = 0;
    int btop_i = 0;

    // We will use CdDeltaReduction as a temporary (double) buffer to accumulate total_bond_order
    // This is safe because CdDeltaReduction is currently zeroed and its accumulation doesn't start until BondsOMP()
    double * tmp_bond_order = workspace->CdDeltaReduction;

    // We do the same with forceReduction as a temporary (rvec) buffer to accumulate dDeltap_self
    // This is safe because forceReduction is currently zeroed and its accumulation does start until Hydrogen_BondsOMP()
    rvec * tmp_ddelta = workspace->forceReduction;

    /* uncorrected bond orders */
    cutoff = control->bond_cut;

#if defined(_OPENMP)
#pragma omp parallel default(shared)                                    \
  private(atom_i, type_i, start_i, end_i, sbp_i, btop_i, ihb, ihb_top,  \
          atom_j, type_j, pj, sbp_j, nbr_pj, jhb, twbp)
#endif
    {

      int nthreads = control->nthreads;
      int tid = get_tid();
      long reductionOffset = (bigint)system->N * tid;
      long totalReductionSize = (bigint)system->N * nthreads;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) reduction(+:num_bonds)
#endif
      for (int i = 0; i < system->N; ++i) {
        atom_i = &(system->my_atoms[i]);
        type_i  = atom_i->type;
        if (type_i < 0) continue;
        sbp_i = &(system->reax_param.sbp[type_i]);

        start_i = Start_Index(i, far_nbrs);
        end_i   = End_Index(i, far_nbrs);

        for (pj = start_i; pj < end_i; ++pj) {
          nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
          if (nbr_pj->d <= cutoff) {
            int j = nbr_pj->nbr;
            atom_j = &(system->my_atoms[j]);
            type_j = atom_j->type;
            if (type_j < 0) continue;
            sbp_j = &(system->reax_param.sbp[type_j]);
            twbp = &(system->reax_param.tbp[type_i][type_j]);

            // Trying to minimize time spent in critical section by moving initial part of BOp()
            // outside of critical section.

            // Start top portion of BOp()
            double C12, C34, C56;
            double BO, BO_s, BO_pi, BO_pi2;
            double bo_cut = control->bo_cut;

            if (sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
              C12 = twbp->p_bo1 * pow(nbr_pj->d / twbp->r_s, twbp->p_bo2);
              BO_s = (1.0 + bo_cut) * exp(C12);
            }
            else BO_s = C12 = 0.0;

            if (sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
              C34 = twbp->p_bo3 * pow(nbr_pj->d / twbp->r_p, twbp->p_bo4);
              BO_pi = exp(C34);
            }
            else BO_pi = C34 = 0.0;

            if (sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
              C56 = twbp->p_bo5 * pow(nbr_pj->d / twbp->r_pp, twbp->p_bo6);
              BO_pi2= exp(C56);
            }
            else BO_pi2 = C56 = 0.0;

            /* Initially BO values are the uncorrected ones, page 1 */
            BO = BO_s + BO_pi + BO_pi2;
            // End top portion of BOp()

            if (BO >= bo_cut) {
              int btop_j;

              // Update indices in critical section
#if defined(_OPENMP)
#pragma omp critical
#endif
              {
                btop_i = End_Index(i, bonds);
                btop_j = End_Index(j, bonds);
                Set_End_Index(j, btop_j+1, bonds);
                Set_End_Index(i, btop_i+1, bonds);
              } // omp critical

              // Finish remaining BOp() work
              BOp_OMP(workspace, bonds, bo_cut,
                      i , btop_i, nbr_pj, sbp_i, sbp_j, twbp, btop_j,
                      C12, C34, C56, BO, BO_s, BO_pi, BO_pi2);

              bond_data * ibond = &(bonds->select.bond_list[btop_i]);
              bond_order_data * bo_ij = &(ibond->bo_data);

              bond_data * jbond = &(bonds->select.bond_list[btop_j]);
              bond_order_data * bo_ji = &(jbond->bo_data);

              workspace->total_bond_order[i]      += bo_ij->BO;
              tmp_bond_order[reductionOffset + j] += bo_ji->BO;

              rvec_Add(workspace->dDeltap_self[i],      bo_ij->dBOp);
              rvec_Add(tmp_ddelta[reductionOffset + j], bo_ji->dBOp);

              btop_i++;
              num_bonds++;
            } // if (BO>=bo_cut)

          } // if (cutoff)

        } // for (pj)
      } // for (i)

      // Need to wait for all indices and tmp arrays accumulated.
#if defined(_OPENMP)
#pragma omp barrier
#endif

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
      for (int i=0; i<system->N; i++)
        for (int t=0; t<nthreads; t++) {
          const int indx = t*system->N + i;
          workspace->dDeltap_self[i][0]  += tmp_ddelta[indx][0];
          workspace->dDeltap_self[i][1]  += tmp_ddelta[indx][1];
          workspace->dDeltap_self[i][2]  += tmp_ddelta[indx][2];
          workspace->total_bond_order[i] += tmp_bond_order[indx];
        }

      /* hydrogen bond list */
      if (control->hbond_cut > 0) {
        cutoff = control->hbond_cut;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50) reduction(+ : num_hbonds)
#endif
        for (int i = 0; i < system->n; ++i) {
          atom_i = &(system->my_atoms[i]);
          type_i  = atom_i->type;
          if (type_i < 0) continue;
          sbp_i = &(system->reax_param.sbp[type_i]);
          ihb = sbp_i->p_hbond;

#if defined(_OPENMP)
#pragma omp critical
#endif
          {

            if (ihb == 1 || ihb == 2) {
              start_i = Start_Index(i, far_nbrs);
              end_i   = End_Index(i, far_nbrs);

              for (pj = start_i; pj < end_i; ++pj) {
                nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
                j = nbr_pj->nbr;
                atom_j = &(system->my_atoms[j]);
                type_j = atom_j->type;
                if (type_j < 0) continue;
                sbp_j = &(system->reax_param.sbp[type_j]);
                jhb = sbp_j->p_hbond;

                if (nbr_pj->d <= control->hbond_cut) {
                  int iflag = 0;
                  int jflag = 0;

                  if (ihb==1 && jhb==2) iflag = 1;
                  else if (j<system->n && ihb == 2 && jhb == 1) jflag = 1;

                  if (iflag || jflag) {
                    if (iflag) {
                      ihb_top = End_Index(atom_i->Hindex, hbonds);
                      Set_End_Index(atom_i->Hindex, ihb_top+1, hbonds);
                    } else if (jflag) {
                      jhb_top = End_Index(atom_j->Hindex, hbonds);
                      Set_End_Index(atom_j->Hindex, jhb_top+1, hbonds);
                    }

                    if (iflag) {
                      hbonds->select.hbond_list[ihb_top].nbr = j;
                      hbonds->select.hbond_list[ihb_top].scl = 1;
                      hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
                    } else if (jflag) {
                      hbonds->select.hbond_list[jhb_top].nbr = i;
                      hbonds->select.hbond_list[jhb_top].scl = -1;
                      hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
                    }

                    num_hbonds++;
                  } // if (iflag || jflag)

                }
              }
            }

          } // omp critical
        }

      } // if (control->hbond > 0)

      // Zero buffers for others to use as intended.
#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
      for (int i=0; i<totalReductionSize; i++) {
        tmp_ddelta[i][0]  = 0.0;
        tmp_ddelta[i][1]  = 0.0;
        tmp_ddelta[i][2]  = 0.0;
        tmp_bond_order[i] = 0.0;
      }

    } // omp

    workspace->realloc.num_bonds = num_bonds;
    workspace->realloc.num_hbonds = num_hbonds;

    Validate_ListsOMP(system, lists, data->step,
                      system->n, system->N, system->numH);
  }

/* ---------------------------------------------------------------------- */

  void Compute_ForcesOMP(reax_system *system, control_params *control,
                          simulation_data *data, storage *workspace,
                          reax_list **lists)
  {
    // Init Forces
    Init_Forces_noQEq_OMP(system, control, data, workspace, lists);

    // Bonded Interactions
    Compute_Bonded_ForcesOMP(system, control, data, workspace, lists);

    // Nonbonded Interactions
    Compute_NonBonded_ForcesOMP(system, control, data, workspace, lists);

    // Total Force
    Compute_Total_ForceOMP(system, control, workspace, lists);
  }
}
