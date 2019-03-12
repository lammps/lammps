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
#include "reaxc_allocate.h"
#include "reaxc_list.h"
#include "reaxc_reset_tools.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "lammps.h"
#include "error.h"
using namespace LAMMPS_NS;

/* allocate space for my_atoms
   important: we cannot know the exact number of atoms that will fall into a
   process's box throughout the whole simulation. therefore
   we need to make upper bound estimates for various data structures */
int PreAllocate_Space( LAMMPS_NS::LAMMPS* lmp, reax_system *system, control_params * /*control*/,
                       storage * workspace, MPI_Comm comm )
{
  int mincap = system->mincap;
  double safezone = system->safezone;

  // determine the local and total capacity

  system->local_cap = MAX( (int)(system->n * safezone), mincap );
  system->total_cap = MAX( (int)(system->N * safezone), mincap );

  system->my_atoms = (reax_atom*)
    scalloc(lmp,  system->total_cap, sizeof(reax_atom), "my_atoms", comm );

  // Nullify some arrays only used in omp styles
  // Should be safe to do here since called in pair->setup();
#ifdef LMP_USER_OMP
  workspace->CdDeltaReduction = NULL;
  workspace->forceReduction = NULL;
  workspace->valence_angle_atom_myoffset = NULL;
  workspace->my_ext_pressReduction = NULL;
#else
  LMP_UNUSED_PARAM(workspace);
#endif

  return SUCCESS;
}


/*************       system        *************/

int Allocate_System( reax_system *system, int /*local_cap*/, int total_cap,
                     char * /*msg*/ )
{
  system->my_atoms = (reax_atom*)
    realloc( system->my_atoms, total_cap*sizeof(reax_atom) );

  return SUCCESS;
}


void DeAllocate_System( LAMMPS_NS::LAMMPS *lmp, reax_system *system )
{
  int i, j, k;
  int ntypes;
  reax_interaction *ff_params;

  // dealloocate the atom list
  sfree(lmp,  system->my_atoms, "system->my_atoms" );

  // deallocate the ffield parameters storage
  ff_params = &(system->reax_param);
  ntypes = ff_params->num_atom_types;

  sfree(lmp,  ff_params->gp.l, "ff:globals" );

  for( i = 0; i < ntypes; ++i ) {
    for( j = 0; j < ntypes; ++j ) {
      for( k = 0; k < ntypes; ++k ) {
        sfree(lmp,  ff_params->fbp[i][j][k], "ff:fbp[i,j,k]" );
      }
      sfree(lmp,  ff_params->fbp[i][j], "ff:fbp[i,j]" );
      sfree(lmp,  ff_params->thbp[i][j], "ff:thbp[i,j]" );
      sfree(lmp,  ff_params->hbp[i][j], "ff:hbp[i,j]" );
    }
    sfree(lmp,  ff_params->fbp[i], "ff:fbp[i]" );
    sfree(lmp,  ff_params->thbp[i], "ff:thbp[i]" );
    sfree(lmp,  ff_params->hbp[i], "ff:hbp[i]" );
    sfree(lmp,  ff_params->tbp[i], "ff:tbp[i]" );
  }
  sfree(lmp,  ff_params->fbp, "ff:fbp" );
  sfree(lmp,  ff_params->thbp, "ff:thbp" );
  sfree(lmp,  ff_params->hbp, "ff:hbp" );
  sfree(lmp,  ff_params->tbp, "ff:tbp" );
  sfree(lmp,  ff_params->sbp, "ff:sbp" );
}


/*************       workspace        *************/
void DeAllocate_Workspace( LAMMPS_NS::LAMMPS* lmp, control_params * /*control*/, storage *workspace )
{
  int i;

  if (!workspace->allocated)
    return;

  workspace->allocated = 0;

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    sfree(lmp,  workspace->tmp_dbl[i], "tmp_dbl[i]" );
    sfree(lmp,  workspace->tmp_rvec[i], "tmp_rvec[i]" );
    sfree(lmp,  workspace->tmp_rvec2[i], "tmp_rvec2[i]" );
  }

  /* bond order storage */
  sfree(lmp,  workspace->within_bond_box, "skin" );
  sfree(lmp,  workspace->total_bond_order, "total_bo" );
  sfree(lmp,  workspace->Deltap, "Deltap" );
  sfree(lmp,  workspace->Deltap_boc, "Deltap_boc" );
  sfree(lmp,  workspace->dDeltap_self, "dDeltap_self" );
  sfree(lmp,  workspace->Delta, "Delta" );
  sfree(lmp,  workspace->Delta_lp, "Delta_lp" );
  sfree(lmp,  workspace->Delta_lp_temp, "Delta_lp_temp" );
  sfree(lmp,  workspace->dDelta_lp, "dDelta_lp" );
  sfree(lmp,  workspace->dDelta_lp_temp, "dDelta_lp_temp" );
  sfree(lmp,  workspace->Delta_e, "Delta_e" );
  sfree(lmp,  workspace->Delta_boc, "Delta_boc" );
  sfree(lmp,  workspace->Delta_val, "Delta_val" );
  sfree(lmp,  workspace->nlp, "nlp" );
  sfree(lmp,  workspace->nlp_temp, "nlp_temp" );
  sfree(lmp,  workspace->Clp, "Clp" );
  sfree(lmp,  workspace->vlpex, "vlpex" );
  sfree(lmp,  workspace->bond_mark, "bond_mark" );
  sfree(lmp,  workspace->done_after, "done_after" );

  /* QEq storage */
  sfree(lmp,  workspace->Hdia_inv, "Hdia_inv" );
  sfree(lmp,  workspace->b_s, "b_s" );
  sfree(lmp,  workspace->b_t, "b_t" );
  sfree(lmp,  workspace->b_prc, "b_prc" );
  sfree(lmp,  workspace->b_prm, "b_prm" );
  sfree(lmp,  workspace->s, "s" );
  sfree(lmp,  workspace->t, "t" );
  sfree(lmp,  workspace->droptol, "droptol" );
  sfree(lmp,  workspace->b, "b" );
  sfree(lmp,  workspace->x, "x" );

  /* GMRES storage */
  for( i = 0; i < RESTART+1; ++i ) {
    sfree(lmp,  workspace->h[i], "h[i]" );
    sfree(lmp,  workspace->v[i], "v[i]" );
  }
  sfree(lmp,  workspace->h, "h" );
  sfree(lmp,  workspace->v, "v" );
  sfree(lmp,  workspace->y, "y" );
  sfree(lmp,  workspace->z, "z" );
  sfree(lmp,  workspace->g, "g" );
  sfree(lmp,  workspace->hs, "hs" );
  sfree(lmp,  workspace->hc, "hc" );
  /* CG storage */
  sfree(lmp,  workspace->r, "r" );
  sfree(lmp,  workspace->d, "d" );
  sfree(lmp,  workspace->q, "q" );
  sfree(lmp,  workspace->p, "p" );
  sfree(lmp,  workspace->r2, "r2" );
  sfree(lmp,  workspace->d2, "d2" );
  sfree(lmp,  workspace->q2, "q2" );
  sfree(lmp,  workspace->p2, "p2" );

  /* integrator storage */
  sfree(lmp,  workspace->v_const, "v_const" );

  /* force related storage */
  sfree(lmp,  workspace->f, "f" );
  sfree(lmp,  workspace->CdDelta, "CdDelta" );

  /* reductions */
#ifdef LMP_USER_OMP
  if (workspace->CdDeltaReduction) sfree(lmp,  workspace->CdDeltaReduction, "cddelta_reduce" );
  if (workspace->forceReduction) sfree(lmp,  workspace->forceReduction, "f_reduce" );
  if (workspace->valence_angle_atom_myoffset) sfree(lmp,  workspace->valence_angle_atom_myoffset, "valence_angle_atom_myoffset");
  if (workspace->my_ext_pressReduction) sfree(lmp,  workspace->my_ext_pressReduction, "ext_press_reduce");
#endif
}


int Allocate_Workspace( LAMMPS_NS::LAMMPS* lmp, reax_system * /*system*/, control_params * control,
                        storage *workspace, int local_cap, int total_cap,
                        MPI_Comm comm, char * /*msg*/ )
{
  int i, total_real, total_rvec, local_rvec;

  workspace->allocated = 1;
  total_real = total_cap * sizeof(double);
  total_rvec = total_cap * sizeof(rvec);
  local_rvec = local_cap * sizeof(rvec);

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    workspace->tmp_dbl[i] = (double*)
      scalloc(lmp,  total_cap, sizeof(double), "tmp_dbl", comm );
    workspace->tmp_rvec[i] = (rvec*)
      scalloc(lmp,  total_cap, sizeof(rvec), "tmp_rvec", comm );
    workspace->tmp_rvec2[i] = (rvec2*)
      scalloc(lmp,  total_cap, sizeof(rvec2), "tmp_rvec2", comm );
  }

  /* bond order related storage  */
  workspace->within_bond_box = (int*)
    scalloc(lmp,  total_cap, sizeof(int), "skin", comm );
  workspace->total_bond_order = (double*) smalloc(lmp,  total_real, "total_bo", comm );
  workspace->Deltap = (double*) smalloc(lmp,  total_real, "Deltap", comm );
  workspace->Deltap_boc = (double*) smalloc(lmp,  total_real, "Deltap_boc", comm );
  workspace->dDeltap_self = (rvec*) smalloc(lmp,  total_rvec, "dDeltap_self", comm );
  workspace->Delta = (double*) smalloc(lmp,  total_real, "Delta", comm );
  workspace->Delta_lp = (double*) smalloc(lmp,  total_real, "Delta_lp", comm );
  workspace->Delta_lp_temp = (double*)
    smalloc(lmp,  total_real, "Delta_lp_temp", comm );
  workspace->dDelta_lp = (double*) smalloc(lmp,  total_real, "dDelta_lp", comm );
  workspace->dDelta_lp_temp = (double*)
    smalloc(lmp,  total_real, "dDelta_lp_temp", comm );
  workspace->Delta_e = (double*) smalloc(lmp,  total_real, "Delta_e", comm );
  workspace->Delta_boc = (double*) smalloc(lmp,  total_real, "Delta_boc", comm );
  workspace->Delta_val = (double*) smalloc(lmp,  total_real, "Delta_val", comm );
  workspace->nlp = (double*) smalloc(lmp,  total_real, "nlp", comm );
  workspace->nlp_temp = (double*) smalloc(lmp,  total_real, "nlp_temp", comm );
  workspace->Clp = (double*) smalloc(lmp,  total_real, "Clp", comm );
  workspace->vlpex = (double*) smalloc(lmp,  total_real, "vlpex", comm );
  workspace->bond_mark = (int*)
    scalloc(lmp,  total_cap, sizeof(int), "bond_mark", comm );
  workspace->done_after = (int*)
    scalloc(lmp,  total_cap, sizeof(int), "done_after", comm );

  /* QEq storage */
  workspace->Hdia_inv = (double*)
    scalloc(lmp,  total_cap, sizeof(double), "Hdia_inv", comm );
  workspace->b_s = (double*) scalloc(lmp,  total_cap, sizeof(double), "b_s", comm );
  workspace->b_t = (double*) scalloc(lmp,  total_cap, sizeof(double), "b_t", comm );
  workspace->b_prc = (double*) scalloc(lmp,  total_cap, sizeof(double), "b_prc", comm );
  workspace->b_prm = (double*) scalloc(lmp,  total_cap, sizeof(double), "b_prm", comm );
  workspace->s = (double*) scalloc(lmp,  total_cap, sizeof(double), "s", comm );
  workspace->t = (double*) scalloc(lmp,  total_cap, sizeof(double), "t", comm );
  workspace->droptol = (double*)
    scalloc(lmp,  total_cap, sizeof(double), "droptol", comm );
  workspace->b = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "b", comm );
  workspace->x = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "x", comm );

  /* GMRES storage */
  workspace->y = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "y", comm );
  workspace->z = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "z", comm );
  workspace->g = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "g", comm );
  workspace->h = (double**) scalloc(lmp,  RESTART+1, sizeof(double*), "h", comm );
  workspace->hs = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "hs", comm );
  workspace->hc = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "hc", comm );
  workspace->v = (double**) scalloc(lmp,  RESTART+1, sizeof(double*), "v", comm );

  for( i = 0; i < RESTART+1; ++i ) {
    workspace->h[i] = (double*) scalloc(lmp,  RESTART+1, sizeof(double), "h[i]", comm );
    workspace->v[i] = (double*) scalloc(lmp,  total_cap, sizeof(double), "v[i]", comm );
  }

  /* CG storage */
  workspace->r = (double*) scalloc(lmp,  total_cap, sizeof(double), "r", comm );
  workspace->d = (double*) scalloc(lmp,  total_cap, sizeof(double), "d", comm );
  workspace->q = (double*) scalloc(lmp,  total_cap, sizeof(double), "q", comm );
  workspace->p = (double*) scalloc(lmp,  total_cap, sizeof(double), "p", comm );
  workspace->r2 = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "r2", comm );
  workspace->d2 = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "d2", comm );
  workspace->q2 = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "q2", comm );
  workspace->p2 = (rvec2*) scalloc(lmp,  total_cap, sizeof(rvec2), "p2", comm );

  /* integrator storage */
  workspace->v_const = (rvec*) smalloc(lmp,  local_rvec, "v_const", comm );

  /* force related storage */
  workspace->f = (rvec*) scalloc(lmp,  total_cap, sizeof(rvec), "f", comm );
  workspace->CdDelta = (double*)
    scalloc(lmp,  total_cap, sizeof(double), "CdDelta", comm );

  // storage for reductions with multiple threads
#ifdef LMP_USER_OMP
  workspace->CdDeltaReduction = (double *) scalloc(lmp, sizeof(double), total_cap*control->nthreads,
                                                 "cddelta_reduce", comm);

  workspace->forceReduction = (rvec *) scalloc(lmp, sizeof(rvec), total_cap*control->nthreads,
                                               "forceReduction", comm);

  workspace->valence_angle_atom_myoffset = (int *) scalloc(lmp, sizeof(int), total_cap, "valence_angle_atom_myoffset", comm);
  workspace->my_ext_pressReduction = (rvec *) calloc(sizeof(rvec), control->nthreads);
#else
  LMP_UNUSED_PARAM(control);
#endif

  return SUCCESS;
}


static void Reallocate_Neighbor_List( LAMMPS* lmp, reax_list *far_nbrs, int n,
                                      int num_intrs, MPI_Comm comm )
{
  Delete_List( lmp, far_nbrs, comm );
  if(!Make_List( lmp, n, num_intrs, TYP_FAR_NEIGHBOR, far_nbrs, comm )){
    lmp->error->one(FLERR,"Problem in initializing far neighbors list");
  }
}


static int Reallocate_HBonds_List( LAMMPS *lmp, reax_system *system, reax_list *hbonds,
                                   MPI_Comm comm )
{
  int i, total_hbonds;

  int mincap = system->mincap;
  double saferzone = system->saferzone;

  total_hbonds = 0;
  for( i = 0; i < system->n; ++i )
    if ((system->my_atoms[i].Hindex) >= 0) {
      total_hbonds += system->my_atoms[i].num_hbonds;
    }
  total_hbonds = (int)(MAX( total_hbonds*saferzone, mincap*MIN_HBONDS ));

  Delete_List( lmp, hbonds, comm );
  if (!Make_List( lmp, system->Hcap, total_hbonds, TYP_HBOND, hbonds, comm )) {
    lmp->error->one(FLERR, "Not enough space for hydrogen bonds list");
  }

  return total_hbonds;
}


static int Reallocate_Bonds_List( LAMMPS *lmp, reax_system *system, reax_list *bonds,
                                  int *total_bonds, int *est_3body,
                                  MPI_Comm comm )
{
  int i;

  int mincap = system->mincap;
  double safezone = system->safezone;

  *total_bonds = 0;
  *est_3body = 0;
  for( i = 0; i < system->N; ++i ){
    *est_3body += SQR(system->my_atoms[i].num_bonds);
    *total_bonds += system->my_atoms[i].num_bonds;
  }
  *total_bonds = (int)(MAX( *total_bonds * safezone, mincap*MIN_BONDS ));

#ifdef LMP_USER_OMP
  if (system->omp_active)
    for (i = 0; i < bonds->num_intrs; ++i)
      sfree(lmp, bonds->select.bond_list[i].bo_data.CdboReduction, "CdboReduction");
#endif

  Delete_List( lmp, bonds, comm );
  if(!Make_List(lmp, system->total_cap, *total_bonds, TYP_BOND, bonds, comm)) {
    lmp->error->one(FLERR, "Not enough space for bonds list");
  }

#ifdef LMP_USER_OMP
#if defined(_OPENMP)
  int nthreads = omp_get_num_threads();
#else
  int nthreads = 1;
#endif

  if (system->omp_active)
    for (i = 0; i < bonds->num_intrs; ++i)
      bonds->select.bond_list[i].bo_data.CdboReduction =
        (double*) smalloc(lmp, sizeof(double)*nthreads, "CdboReduction", comm);
#endif

  return SUCCESS;
}


void ReAllocate( LAMMPS *lmp, reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists,
                 mpi_datatypes *mpi_data )
{
  int num_bonds, est_3body, Hflag, ret;
  int renbr, newsize;
  reallocate_data *realloc;
  reax_list *far_nbrs;
  MPI_Comm comm;
  char msg[200];

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  realloc = &(workspace->realloc);
  comm = mpi_data->world;

  if( system->n >= DANGER_ZONE * system->local_cap ||
      (0 && system->n <= LOOSE_ZONE * system->local_cap) ) {
    system->local_cap = MAX( (int)(system->n * safezone), mincap );
  }

  int Nflag = 0;
  if( system->N >= DANGER_ZONE * system->total_cap ||
      (0 && system->N <= LOOSE_ZONE * system->total_cap) ) {
    Nflag = 1;
    system->total_cap = MAX( (int)(system->N * safezone), mincap );
  }

  if (Nflag) {
    /* system */
    ret = Allocate_System( system, system->local_cap, system->total_cap, msg );
    if (ret != SUCCESS) {
      char errmsg[128];
      snprintf(errmsg, 128, "Not enough space for atom_list: total_cap=%d", system->total_cap);
      lmp->error->one(FLERR, errmsg);
    }

    /* workspace */
    DeAllocate_Workspace( lmp, control, workspace );
    ret = Allocate_Workspace( lmp, system, control, workspace, system->local_cap,
                              system->total_cap, comm, msg );
    if (ret != SUCCESS) {
      char errmsg[128];
      snprintf(errmsg, 128, "Not enough space for workspace: local_cap=%d total_cap=%d", system->local_cap, system->total_cap);
      lmp->error->one(FLERR, errmsg);
    }
  }


  renbr = (data->step - data->prev_steps) % control->reneighbor == 0;
  /* far neighbors */
  if (renbr) {
    far_nbrs = *lists + FAR_NBRS;

    if (Nflag || realloc->num_far >= far_nbrs->num_intrs * DANGER_ZONE) {
      if (realloc->num_far > far_nbrs->num_intrs) {
        char errmsg[128];
        snprintf(errmsg, 128, "step%d-ran out of space on far_nbrs: top=%d, max=%d", data->step, realloc->num_far, far_nbrs->num_intrs);
        lmp->error->one(FLERR, errmsg);
      }

      newsize = static_cast<int>
        (MAX( realloc->num_far*safezone, mincap*MIN_NBRS ));

      Reallocate_Neighbor_List( lmp, far_nbrs, system->total_cap, newsize, comm );
      realloc->num_far = 0;
    }
  }

  /* hydrogen bonds list */
  if (control->hbond_cut > 0) {
    Hflag = 0;
    if( system->numH >= DANGER_ZONE * system->Hcap ||
        (0 && system->numH <= LOOSE_ZONE * system->Hcap) ) {
      Hflag = 1;
      system->Hcap = int(MAX( system->numH * saferzone, mincap ));
    }

    if (Hflag || realloc->hbonds) {
      ret = Reallocate_HBonds_List( lmp, system, (*lists)+HBONDS, comm );
      realloc->hbonds = 0;
    }
  }

  /* bonds list */
  num_bonds = est_3body = -1;
  if (Nflag || realloc->bonds) {
    Reallocate_Bonds_List( lmp, system, (*lists)+BONDS, &num_bonds,
                           &est_3body, comm );
    realloc->bonds = 0;
    realloc->num_3body = MAX( realloc->num_3body, est_3body ) * 2;
  }

  /* 3-body list */
  if (realloc->num_3body > 0) {
    Delete_List( lmp, (*lists)+THREE_BODIES, comm );

    if (num_bonds == -1)
      num_bonds = ((*lists)+BONDS)->num_intrs;

    realloc->num_3body = (int)(MAX(realloc->num_3body*safezone, MIN_3BODIES));

    if( !Make_List( lmp, num_bonds, realloc->num_3body, TYP_THREE_BODY,
                    (*lists)+THREE_BODIES, comm ) ) {
      lmp->error->one(FLERR, "Problem in initializing angles list");
    }
    realloc->num_3body = -1;
  }

}
