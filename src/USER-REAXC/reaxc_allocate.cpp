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

#include "reaxc_allocate.h"
#include <cstdlib>
#include "reaxc_defs.h"
#include "reaxc_list.h"
#include "reaxc_tool_box.h"

#if defined(LMP_USER_OMP) && defined(_OPENMP)
#include <omp.h>
#endif

#include "error.h"

/* allocate space for my_atoms
   important: we cannot know the exact number of atoms that will fall into a
   process's box throughout the whole simulation. therefore
   we need to make upper bound estimates for various data structures */
int PreAllocate_Space( reax_system *system, control_params * /*control*/,
                       storage * workspace )
{
  int mincap = system->mincap;
  double safezone = system->safezone;

  // determine the local and total capacity

  system->local_cap = MAX( (int)(system->n * safezone), mincap );
  system->total_cap = MAX( (int)(system->N * safezone), mincap );

  system->my_atoms = (reax_atom*)
    scalloc(system->error_ptr,  system->total_cap, sizeof(reax_atom), "my_atoms");

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


void DeAllocate_System( reax_system *system )
{
  int i, j, k;
  int ntypes;
  reax_interaction *ff_params;

  // deallocate the atom list
  sfree(system->error_ptr,  system->my_atoms, "system->my_atoms" );

  // deallocate the ffield parameters storage
  ff_params = &(system->reax_param);
  ntypes = ff_params->num_atom_types;

  sfree(system->error_ptr,  ff_params->gp.l, "ff:globals" );

  for( i = 0; i < ntypes; ++i ) {
    for( j = 0; j < ntypes; ++j ) {
      for( k = 0; k < ntypes; ++k ) {
        sfree(system->error_ptr,  ff_params->fbp[i][j][k], "ff:fbp[i,j,k]" );
      }
      sfree(system->error_ptr,  ff_params->fbp[i][j], "ff:fbp[i,j]" );
      sfree(system->error_ptr,  ff_params->thbp[i][j], "ff:thbp[i,j]" );
      sfree(system->error_ptr,  ff_params->hbp[i][j], "ff:hbp[i,j]" );
    }
    sfree(system->error_ptr,  ff_params->fbp[i], "ff:fbp[i]" );
    sfree(system->error_ptr,  ff_params->thbp[i], "ff:thbp[i]" );
    sfree(system->error_ptr,  ff_params->hbp[i], "ff:hbp[i]" );
    sfree(system->error_ptr,  ff_params->tbp[i], "ff:tbp[i]" );
  }
  sfree(system->error_ptr,  ff_params->fbp, "ff:fbp" );
  sfree(system->error_ptr,  ff_params->thbp, "ff:thbp" );
  sfree(system->error_ptr,  ff_params->hbp, "ff:hbp" );
  sfree(system->error_ptr,  ff_params->tbp, "ff:tbp" );
  sfree(system->error_ptr,  ff_params->sbp, "ff:sbp" );
}


/*************       workspace        *************/
void DeAllocate_Workspace( control_params * control, storage *workspace )
{
  int i;

  if (!workspace->allocated)
    return;

  workspace->allocated = 0;

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    sfree(control->error_ptr,  workspace->tmp_dbl[i], "tmp_dbl[i]" );
    sfree(control->error_ptr,  workspace->tmp_rvec[i], "tmp_rvec[i]" );
    sfree(control->error_ptr,  workspace->tmp_rvec2[i], "tmp_rvec2[i]" );
  }

  /* bond order storage */
  sfree(control->error_ptr,  workspace->within_bond_box, "skin" );
  sfree(control->error_ptr,  workspace->total_bond_order, "total_bo" );
  sfree(control->error_ptr,  workspace->Deltap, "Deltap" );
  sfree(control->error_ptr,  workspace->Deltap_boc, "Deltap_boc" );
  sfree(control->error_ptr,  workspace->dDeltap_self, "dDeltap_self" );
  sfree(control->error_ptr,  workspace->Delta, "Delta" );
  sfree(control->error_ptr,  workspace->Delta_lp, "Delta_lp" );
  sfree(control->error_ptr,  workspace->Delta_lp_temp, "Delta_lp_temp" );
  sfree(control->error_ptr,  workspace->dDelta_lp, "dDelta_lp" );
  sfree(control->error_ptr,  workspace->dDelta_lp_temp, "dDelta_lp_temp" );
  sfree(control->error_ptr,  workspace->Delta_e, "Delta_e" );
  sfree(control->error_ptr,  workspace->Delta_boc, "Delta_boc" );
  sfree(control->error_ptr,  workspace->Delta_val, "Delta_val" );
  sfree(control->error_ptr,  workspace->nlp, "nlp" );
  sfree(control->error_ptr,  workspace->nlp_temp, "nlp_temp" );
  sfree(control->error_ptr,  workspace->Clp, "Clp" );
  sfree(control->error_ptr,  workspace->vlpex, "vlpex" );
  sfree(control->error_ptr,  workspace->bond_mark, "bond_mark" );
  sfree(control->error_ptr,  workspace->done_after, "done_after" );

  /* QEq storage */
  sfree(control->error_ptr,  workspace->Hdia_inv, "Hdia_inv" );
  sfree(control->error_ptr,  workspace->b_s, "b_s" );
  sfree(control->error_ptr,  workspace->b_t, "b_t" );
  sfree(control->error_ptr,  workspace->b_prc, "b_prc" );
  sfree(control->error_ptr,  workspace->b_prm, "b_prm" );
  sfree(control->error_ptr,  workspace->s, "s" );
  sfree(control->error_ptr,  workspace->t, "t" );
  sfree(control->error_ptr,  workspace->droptol, "droptol" );
  sfree(control->error_ptr,  workspace->b, "b" );
  sfree(control->error_ptr,  workspace->x, "x" );

  /* GMRES storage */
  for( i = 0; i < RESTART+1; ++i ) {
    sfree(control->error_ptr,  workspace->h[i], "h[i]" );
    sfree(control->error_ptr,  workspace->v[i], "v[i]" );
  }
  sfree(control->error_ptr,  workspace->h, "h" );
  sfree(control->error_ptr,  workspace->v, "v" );
  sfree(control->error_ptr,  workspace->y, "y" );
  sfree(control->error_ptr,  workspace->z, "z" );
  sfree(control->error_ptr,  workspace->g, "g" );
  sfree(control->error_ptr,  workspace->hs, "hs" );
  sfree(control->error_ptr,  workspace->hc, "hc" );
  /* CG storage */
  sfree(control->error_ptr,  workspace->r, "r" );
  sfree(control->error_ptr,  workspace->d, "d" );
  sfree(control->error_ptr,  workspace->q, "q" );
  sfree(control->error_ptr,  workspace->p, "p" );
  sfree(control->error_ptr,  workspace->r2, "r2" );
  sfree(control->error_ptr,  workspace->d2, "d2" );
  sfree(control->error_ptr,  workspace->q2, "q2" );
  sfree(control->error_ptr,  workspace->p2, "p2" );

  /* integrator storage */
  sfree(control->error_ptr,  workspace->v_const, "v_const" );

  /* force related storage */
  sfree(control->error_ptr,  workspace->f, "f" );
  sfree(control->error_ptr,  workspace->CdDelta, "CdDelta" );

  /* reductions */
#ifdef LMP_USER_OMP
  if (workspace->CdDeltaReduction) sfree(control->error_ptr,  workspace->CdDeltaReduction, "cddelta_reduce" );
  if (workspace->forceReduction) sfree(control->error_ptr,  workspace->forceReduction, "f_reduce" );
  if (workspace->valence_angle_atom_myoffset) sfree(control->error_ptr,  workspace->valence_angle_atom_myoffset, "valence_angle_atom_myoffset");
  if (workspace->my_ext_pressReduction) sfree(control->error_ptr,  workspace->my_ext_pressReduction, "ext_press_reduce");
#endif
}


int Allocate_Workspace( reax_system * /*system*/, control_params * control,
                        storage *workspace, int local_cap, int total_cap,
                        char * /*msg*/ )
{
  int i, total_real, total_rvec, local_rvec;

  workspace->allocated = 1;
  total_real = total_cap * sizeof(double);
  total_rvec = total_cap * sizeof(rvec);
  local_rvec = local_cap * sizeof(rvec);

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    workspace->tmp_dbl[i] = (double*)
      scalloc(control->error_ptr,  total_cap, sizeof(double), "tmp_dbl");
    workspace->tmp_rvec[i] = (rvec*)
      scalloc(control->error_ptr,  total_cap, sizeof(rvec), "tmp_rvec");
    workspace->tmp_rvec2[i] = (rvec2*)
      scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "tmp_rvec2");
  }

  /* bond order related storage  */
  workspace->within_bond_box = (int*)
    scalloc(control->error_ptr,  total_cap, sizeof(int), "skin");
  workspace->total_bond_order = (double*) smalloc(control->error_ptr,  total_real, "total_bo");
  workspace->Deltap = (double*) smalloc(control->error_ptr,  total_real, "Deltap");
  workspace->Deltap_boc = (double*) smalloc(control->error_ptr,  total_real, "Deltap_boc");
  workspace->dDeltap_self = (rvec*) smalloc(control->error_ptr,  total_rvec, "dDeltap_self");
  workspace->Delta = (double*) smalloc(control->error_ptr,  total_real, "Delta");
  workspace->Delta_lp = (double*) smalloc(control->error_ptr,  total_real, "Delta_lp");
  workspace->Delta_lp_temp = (double*)
    smalloc(control->error_ptr,  total_real, "Delta_lp_temp");
  workspace->dDelta_lp = (double*) smalloc(control->error_ptr,  total_real, "dDelta_lp");
  workspace->dDelta_lp_temp = (double*)
    smalloc(control->error_ptr,  total_real, "dDelta_lp_temp");
  workspace->Delta_e = (double*) smalloc(control->error_ptr,  total_real, "Delta_e");
  workspace->Delta_boc = (double*) smalloc(control->error_ptr,  total_real, "Delta_boc");
  workspace->Delta_val = (double*) smalloc(control->error_ptr,  total_real, "Delta_val");
  workspace->nlp = (double*) smalloc(control->error_ptr,  total_real, "nlp");
  workspace->nlp_temp = (double*) smalloc(control->error_ptr,  total_real, "nlp_temp");
  workspace->Clp = (double*) smalloc(control->error_ptr,  total_real, "Clp");
  workspace->vlpex = (double*) smalloc(control->error_ptr,  total_real, "vlpex");
  workspace->bond_mark = (int*)
    scalloc(control->error_ptr,  total_cap, sizeof(int), "bond_mark");
  workspace->done_after = (int*)
    scalloc(control->error_ptr,  total_cap, sizeof(int), "done_after");

  /* QEq storage */
  workspace->Hdia_inv = (double*)
    scalloc(control->error_ptr,  total_cap, sizeof(double), "Hdia_inv");
  workspace->b_s = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "b_s");
  workspace->b_t = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "b_t");
  workspace->b_prc = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "b_prc");
  workspace->b_prm = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "b_prm");
  workspace->s = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "s");
  workspace->t = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "t");
  workspace->droptol = (double*)
    scalloc(control->error_ptr,  total_cap, sizeof(double), "droptol");
  workspace->b = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "b");
  workspace->x = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "x");

  /* GMRES storage */
  workspace->y = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "y");
  workspace->z = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "z");
  workspace->g = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "g");
  workspace->h = (double**) scalloc(control->error_ptr,  RESTART+1, sizeof(double*), "h");
  workspace->hs = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "hs");
  workspace->hc = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "hc");
  workspace->v = (double**) scalloc(control->error_ptr,  RESTART+1, sizeof(double*), "v");

  for( i = 0; i < RESTART+1; ++i ) {
    workspace->h[i] = (double*) scalloc(control->error_ptr,  RESTART+1, sizeof(double), "h[i]");
    workspace->v[i] = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "v[i]");
  }

  /* CG storage */
  workspace->r = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "r");
  workspace->d = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "d");
  workspace->q = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "q");
  workspace->p = (double*) scalloc(control->error_ptr,  total_cap, sizeof(double), "p");
  workspace->r2 = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "r2");
  workspace->d2 = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "d2");
  workspace->q2 = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "q2");
  workspace->p2 = (rvec2*) scalloc(control->error_ptr,  total_cap, sizeof(rvec2), "p2");

  /* integrator storage */
  workspace->v_const = (rvec*) smalloc(control->error_ptr,  local_rvec, "v_const");

  /* force related storage */
  workspace->f = (rvec*) scalloc(control->error_ptr,  total_cap, sizeof(rvec), "f");
  workspace->CdDelta = (double*)
    scalloc(control->error_ptr,  total_cap, sizeof(double), "CdDelta");

  // storage for reductions with multiple threads
#ifdef LMP_USER_OMP
  workspace->CdDeltaReduction = (double *) scalloc(control->error_ptr, sizeof(double), total_cap*control->nthreads,
                                                 "cddelta_reduce");

  workspace->forceReduction = (rvec *) scalloc(control->error_ptr, sizeof(rvec), total_cap*control->nthreads,
                                               "forceReduction");

  workspace->valence_angle_atom_myoffset = (int *) scalloc(control->error_ptr, sizeof(int), total_cap, "valence_angle_atom_myoffset");
  workspace->my_ext_pressReduction = (rvec *) calloc(sizeof(rvec), control->nthreads);
#else
  LMP_UNUSED_PARAM(control);
#endif

  return SUCCESS;
}


static void Reallocate_Neighbor_List( reax_list *far_nbrs, int n,
                                      int num_intrs )
{
  Delete_List( far_nbrs);
  if(!Make_List( n, num_intrs, TYP_FAR_NEIGHBOR, far_nbrs )){
    far_nbrs->error_ptr->one(FLERR,"Problem in initializing far neighbors list");
  }
}


static int Reallocate_HBonds_List( reax_system *system, reax_list *hbonds )
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

  Delete_List( hbonds);
  if (!Make_List( system->Hcap, total_hbonds, TYP_HBOND, hbonds )) {
    hbonds->error_ptr->one(FLERR, "Not enough space for hydrogen bonds list");
  }

  return total_hbonds;
}


static int Reallocate_Bonds_List( reax_system *system, reax_list *bonds,
                                  int *total_bonds, int *est_3body )
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
      sfree(system->error_ptr, bonds->select.bond_list[i].bo_data.CdboReduction, "CdboReduction");
#endif

  Delete_List( bonds);
  if(!Make_List(system->total_cap, *total_bonds, TYP_BOND, bonds)) {
    bonds->error_ptr->one(FLERR, "Not enough space for bonds list");
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
        (double*) smalloc(system->error_ptr, sizeof(double)*nthreads, "CdboReduction");
#endif

  return SUCCESS;
}


void ReAllocate( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists )
{
  int num_bonds, est_3body, Hflag, ret;
  int renbr, newsize;
  reallocate_data *realloc;
  reax_list *far_nbrs;
  char msg[200];

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  realloc = &(workspace->realloc);

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
      char errmsg[256];
      snprintf(errmsg, 256, "Not enough space for atom_list: total_cap=%d", system->total_cap);
      system->error_ptr->one(FLERR, errmsg);
    }

    /* workspace */
    DeAllocate_Workspace( control, workspace );
    ret = Allocate_Workspace( system, control, workspace, system->local_cap,
                              system->total_cap, msg );
    if (ret != SUCCESS) {
      char errmsg[256];
      snprintf(errmsg, 256, "Not enough space for workspace: local_cap=%d total_cap=%d", system->local_cap, system->total_cap);
      system->error_ptr->one(FLERR, errmsg);
    }
  }


  renbr = (data->step - data->prev_steps) % control->reneighbor == 0;
  /* far neighbors */
  if (renbr) {
    far_nbrs = *lists + FAR_NBRS;

    if (Nflag || realloc->num_far >= far_nbrs->num_intrs * DANGER_ZONE) {
      if (realloc->num_far > far_nbrs->num_intrs) {
        char errmsg[256];
        snprintf(errmsg, 256, "step%d-ran out of space on far_nbrs: top=%d, max=%d", data->step, realloc->num_far, far_nbrs->num_intrs);
        system->error_ptr->one(FLERR, errmsg);
      }

      newsize = static_cast<int>
        (MAX( realloc->num_far*safezone, mincap*MIN_NBRS ));

      Reallocate_Neighbor_List( far_nbrs, system->total_cap, newsize);
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
      ret = Reallocate_HBonds_List( system, (*lists)+HBONDS);
      realloc->hbonds = 0;
    }
  }

  /* bonds list */
  num_bonds = est_3body = -1;
  if (Nflag || realloc->bonds) {
    Reallocate_Bonds_List( system, (*lists)+BONDS, &num_bonds,
                           &est_3body);
    realloc->bonds = 0;
    realloc->num_3body = MAX( realloc->num_3body, est_3body ) * 2;
  }

  /* 3-body list */
  if (realloc->num_3body > 0) {
    Delete_List( (*lists)+THREE_BODIES);

    if (num_bonds == -1)
      num_bonds = ((*lists)+BONDS)->num_intrs;

    realloc->num_3body = (int)(MAX(realloc->num_3body*safezone, MIN_3BODIES));

    if( !Make_List( num_bonds, realloc->num_3body, TYP_THREE_BODY,
                    (*lists)+THREE_BODIES ) ) {
      system->error_ptr->one(FLERR, "Problem in initializing angles list");
    }
    realloc->num_3body = -1;
  }

}
