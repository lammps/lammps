// clang-format off
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
  <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxff_api.h"

#include "error.h"
#include "memory.h"
#include "pair.h"

namespace ReaxFF {

  /* allocate space for my_atoms
     important: we cannot know the exact number of atoms that will fall into a
     process's box throughout the whole simulation. therefore
     we need to make upper bound estimates for various data structures */
  void PreAllocate_Space(reax_system *system, storage * workspace)
  {
    const int mincap = system->mincap;
    const double safezone = system->safezone;

    // determine the local and total capacity

    system->local_cap = MAX((int)(system->n * safezone), mincap);
    system->total_cap = MAX((int)(system->N * safezone), mincap);

    system->my_atoms = (reax_atom*) scalloc(system->error_ptr,
      system->total_cap, sizeof(reax_atom), "my_atoms");

    // Nullify some arrays only used in omp styles
    // Should be safe to do here since called in pair->setup();

    workspace->CdDeltaReduction = nullptr;
    workspace->forceReduction = nullptr;
    workspace->valence_angle_atom_myoffset = nullptr;
  }

  /*************       system        *************/

  void DeAllocate_System(reax_system *system)
  {
    auto error = system->error_ptr;
    auto memory = system->mem_ptr;

    // deallocate the atom list
    sfree(error, system->my_atoms, "system->my_atoms");

    // deallocate the ffield parameters storage
    memory->destroy(system->reax_param.gp.l);
    memory->destroy(system->reax_param.sbp);
    memory->destroy(system->reax_param.tbp);
    memory->destroy(system->reax_param.thbp);
    memory->destroy(system->reax_param.hbp);
    memory->destroy(system->reax_param.fbp);
  }

  /*************       workspace        *************/
  void DeAllocate_Workspace(control_params *control, storage *workspace)
  {
    if (!workspace->allocated)
      return;

    workspace->allocated = 0;
    auto error = control->error_ptr;

    /* bond order storage */
    sfree(error, workspace->total_bond_order, "total_bo");
    sfree(error, workspace->Deltap, "Deltap");
    sfree(error, workspace->Deltap_boc, "Deltap_boc");
    sfree(error, workspace->dDeltap_self, "dDeltap_self");
    sfree(error, workspace->Delta, "Delta");
    sfree(error, workspace->Delta_lp, "Delta_lp");
    sfree(error, workspace->Delta_lp_temp, "Delta_lp_temp");
    sfree(error, workspace->dDelta_lp, "dDelta_lp");
    sfree(error, workspace->dDelta_lp_temp, "dDelta_lp_temp");
    sfree(error, workspace->Delta_e, "Delta_e");
    sfree(error, workspace->Delta_boc, "Delta_boc");
    sfree(error, workspace->Delta_val, "Delta_val");
    sfree(error, workspace->nlp, "nlp");
    sfree(error, workspace->nlp_temp, "nlp_temp");
    sfree(error, workspace->Clp, "Clp");
    sfree(error, workspace->vlpex, "vlpex");
    sfree(error, workspace->bond_mark, "bond_mark");

    /* force related storage */
    sfree(error, workspace->f, "f");
    sfree(error, workspace->CdDelta, "CdDelta");

    /* reductions */

    if (workspace->CdDeltaReduction)
      sfree(error, workspace->CdDeltaReduction, "cddelta_reduce");
    if (workspace->forceReduction)
      sfree(error, workspace->forceReduction, "f_reduce");
    if (workspace->valence_angle_atom_myoffset)
      sfree(error, workspace->valence_angle_atom_myoffset, "valence_angle_atom_myoffset");
  }

  void Allocate_Workspace(control_params *control, storage *workspace, int total_cap)
  {
    int total_real, total_rvec;
    auto error = control->error_ptr;

    workspace->allocated = 1;
    total_real = total_cap * sizeof(double);
    total_rvec = total_cap * sizeof(rvec);

    /* bond order related storage  */
    workspace->total_bond_order = (double*) smalloc(error, total_real, "total_bo");
    workspace->Deltap = (double*) smalloc(error, total_real, "Deltap");
    workspace->Deltap_boc = (double*) smalloc(error, total_real, "Deltap_boc");
    workspace->dDeltap_self = (rvec*) smalloc(error, total_rvec, "dDeltap_self");
    workspace->Delta = (double*) smalloc(error, total_real, "Delta");
    workspace->Delta_lp = (double*) smalloc(error, total_real, "Delta_lp");
    workspace->Delta_lp_temp = (double*) smalloc(error, total_real, "Delta_lp_temp");
    workspace->dDelta_lp = (double*) smalloc(error, total_real, "dDelta_lp");
    workspace->dDelta_lp_temp = (double*) smalloc(error, total_real, "dDelta_lp_temp");
    workspace->Delta_e = (double*) smalloc(error, total_real, "Delta_e");
    workspace->Delta_boc = (double*) smalloc(error, total_real, "Delta_boc");
    workspace->Delta_val = (double*) smalloc(error, total_real, "Delta_val");
    workspace->nlp = (double*) smalloc(error, total_real, "nlp");
    workspace->nlp_temp = (double*) smalloc(error, total_real, "nlp_temp");
    workspace->Clp = (double*) smalloc(error, total_real, "Clp");
    workspace->vlpex = (double*) smalloc(error, total_real, "vlpex");
    workspace->bond_mark = (int*) scalloc(error, total_cap, sizeof(int), "bond_mark");

    /* force related storage */
    workspace->f = (rvec*) scalloc(error, total_cap, sizeof(rvec), "f");
    workspace->CdDelta = (double*) scalloc(error, total_cap, sizeof(double), "CdDelta");

    // storage for reductions with multiple threads

    workspace->CdDeltaReduction = (double *) scalloc(error,
      sizeof(double), (rc_bigint)total_cap*control->nthreads, "cddelta_reduce");
    workspace->forceReduction = (rvec *) scalloc(error,
      sizeof(rvec), (rc_bigint)total_cap*control->nthreads, "forceReduction");
    workspace->valence_angle_atom_myoffset = (int *) scalloc(error,
     sizeof(int), total_cap, "valence_angle_atom_myoffset");
  }


  static void Reallocate_Neighbor_List(reax_list *far_nbrs, int n, int num_intrs)
  {
    Delete_List(far_nbrs);
    Make_List(n, num_intrs, TYP_FAR_NEIGHBOR, far_nbrs);
  }

  static int Reallocate_HBonds_List(reax_system *system, reax_list *hbonds)
  {
    int i, total_hbonds;

    int mincap = system->mincap;
    double saferzone = system->saferzone;

    total_hbonds = 0;
    for (i = 0; i < system->n; ++i)
      if ((system->my_atoms[i].Hindex) >= 0) {
        total_hbonds += system->my_atoms[i].num_hbonds;
      }
    total_hbonds = (int)(MAX(total_hbonds*saferzone, mincap*system->minhbonds));

    Delete_List(hbonds);
    Make_List(system->Hcap, total_hbonds, TYP_HBOND, hbonds);

    return total_hbonds;
  }

  static void Reallocate_Bonds_List(control_params *control, reax_system *system,
                                    reax_list *bonds, int *total_bonds, int *est_3body)
  {
    int i;

    int mincap = system->mincap;
    double safezone = system->safezone;

    *total_bonds = 0;
    *est_3body = 0;
    for (i = 0; i < system->N; ++i) {
      *est_3body += SQR(system->my_atoms[i].num_bonds);
      *total_bonds += system->my_atoms[i].num_bonds;
    }
    *total_bonds = (int)(MAX(*total_bonds * safezone, mincap*MIN_BONDS));

    if (system->omp_active)
      for (i = 0; i < bonds->num_intrs; ++i)
        sfree(system->error_ptr, bonds->select.bond_list[i].bo_data.CdboReduction, "CdboReduction");

    Delete_List(bonds);
    Make_List(system->total_cap, *total_bonds, TYP_BOND, bonds);

    if (system->omp_active)
      for (i = 0; i < bonds->num_intrs; ++i)
        bonds->select.bond_list[i].bo_data.CdboReduction =
          (double*) smalloc(system->error_ptr, sizeof(double)*control->nthreads, "CdboReduction");
  }

  void ReAllocate(reax_system *system, control_params *control,
                  simulation_data *data, storage *workspace, reax_list **lists)
  {
    int num_bonds, est_3body, Hflag;
    int newsize;
    reax_list *far_nbrs;

    int mincap = system->mincap;
    double safezone = system->safezone;
    double saferzone = system->saferzone;

    auto error = system->error_ptr;
    reallocate_data *wsr = &(workspace->realloc);

    if (system->n >= DANGER_ZONE * system->local_cap ||
        (0 && system->n <= LOOSE_ZONE * system->local_cap)) {
      system->local_cap = MAX((int)(system->n * safezone), mincap);
    }

    int Nflag = 0;
    if (system->N >= DANGER_ZONE * system->total_cap ||
        (0 && system->N <= LOOSE_ZONE * system->total_cap)) {
      Nflag = 1;
      system->total_cap = MAX((int)(system->N * safezone), mincap);
    }

    if (Nflag) {
      /* system */
      system->my_atoms = (reax_atom *)::realloc(system->my_atoms,
        system->total_cap*sizeof(reax_atom));
      /* workspace */
      DeAllocate_Workspace(control, workspace);
      Allocate_Workspace(control, workspace, system->total_cap);
    }

    /* far neighbors */

    far_nbrs = *lists + FAR_NBRS;

    if (Nflag || wsr->num_far >= far_nbrs->num_intrs * DANGER_ZONE) {
      if (wsr->num_far > far_nbrs->num_intrs)
        error->one(FLERR,fmt::format("step{}: ran out of space on far_nbrs: top={}, max={}",
                                   data->step, wsr->num_far, far_nbrs->num_intrs));

      newsize = static_cast<int>
        (MAX(wsr->num_far*safezone, mincap*REAX_MIN_NBRS));

      Reallocate_Neighbor_List(far_nbrs, system->total_cap, newsize);
      wsr->num_far = 0;
    }

    /* hydrogen bonds list */
    if (control->hbond_cut > 0) {
      Hflag = 0;
      if (system->numH >= DANGER_ZONE * system->Hcap ||
          (0 && system->numH <= LOOSE_ZONE * system->Hcap)) {
        Hflag = 1;
        system->Hcap = int(MAX(system->numH * saferzone, mincap));
      }

      if (Hflag || wsr->hbonds) {
        Reallocate_HBonds_List(system, (*lists)+HBONDS);
        wsr->hbonds = 0;
      }
    }

    /* bonds list */
    num_bonds = est_3body = -1;
    if (Nflag || wsr->bonds) {
      Reallocate_Bonds_List(control, system, (*lists)+BONDS, &num_bonds, &est_3body);
      wsr->bonds = 0;
      wsr->num_3body = MAX(wsr->num_3body, est_3body) * 2;


      if (system->omp_active) {
        int nthreads = control->nthreads;
        reax_list *bonds = (*lists)+BONDS;

        for (int i = 0; i < bonds->num_intrs; ++i) {
          sfree(error, bonds->select.bond_list[i].bo_data.CdboReduction,
                "CdboReduction");

          bonds->select.bond_list[i].bo_data.CdboReduction =
            (double*) smalloc(error, sizeof(double)*nthreads, "CdboReduction");
        }
      }
    }

    /* 3-body list */
    if (wsr->num_3body > 0) {
      Delete_List((*lists)+THREE_BODIES);

      if (num_bonds == -1)
        num_bonds = ((*lists)+BONDS)->num_intrs;

      wsr->num_3body = (int)(MAX(wsr->num_3body*safezone, MIN_3BODIES));

      Make_List(num_bonds, wsr->num_3body, TYP_THREE_BODY,
                (*lists)+THREE_BODIES);
      wsr->num_3body = -1;
    }
  }
}
