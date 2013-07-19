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
#if defined(PURE_REAX)
#include "allocate.h"
#include "list.h"
#include "reset_tools.h"
#include "tool_box.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_allocate.h"
#include "reaxc_list.h"
#include "reaxc_reset_tools.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"
#endif


/* allocate space for my_atoms
   important: we cannot know the exact number of atoms that will fall into a
   process's box throughout the whole simulation. therefore
   we need to make upper bound estimates for various data structures */
int PreAllocate_Space( reax_system *system, control_params *control,
                       storage *workspace, MPI_Comm comm )
{
  int mincap = system->mincap;
  double safezone = system->safezone;

  // determine the local and total capacity

  system->local_cap = MAX( (int)(system->n * safezone), mincap );
  system->total_cap = MAX( (int)(system->N * safezone), mincap );

#if defined(DEBUG)
  fprintf( stderr, "p%d: local_cap=%d total_cap=%d\n",
           system->my_rank, system->local_cap, system->total_cap );
#endif

  system->my_atoms = (reax_atom*)
    scalloc( system->total_cap, sizeof(reax_atom), "my_atoms", comm );

  /* space for keeping restriction info, if any */
  // not yet implemented in the parallel version!!!
  // if( control->restrict_bonds ) {
  //   workspace->restricted  = (int*)
  //     scalloc( system->local_cap, sizeof(int), "restricted_atoms", comm );

  //   workspace->restricted_list = (int**)
  //     scalloc( system->local_cap, sizeof(int*), "restricted_list", comm );

  //   for( i = 0; i < system->local_cap; ++i )
  //     workspace->restricted_list[i] = (int*)
  //         scalloc( MAX_RESTRICT, sizeof(int), "restricted_list[i]", comm );
  // }

  return SUCCESS;
}


/*************       system        *************/
inline void reax_atom_Copy( reax_atom *dest, reax_atom *src )
{
  dest->orig_id = src->orig_id;
  dest->type = src->type;
  strcpy( dest->name, src->name );
  rvec_Copy( dest->x, src->x );
  rvec_Copy( dest->v, src->v );
  rvec_Copy( dest->f_old, src->f_old );
  rvec_Copy( dest->s, src->s );
  rvec_Copy( dest->t, src->t );
  dest->Hindex = src->Hindex;
  dest->num_bonds = src->num_bonds;
  dest->num_hbonds = src->num_hbonds;
  dest->numbonds = src->numbonds;
}


void Copy_Atom_List( reax_atom *dest, reax_atom *src, int n )
{
  int i;

  for( i = 0; i < n; ++i )
    memcpy( dest+i, src+i, sizeof(reax_atom) );
}


int Allocate_System( reax_system *system, int local_cap, int total_cap,
                     char *msg )
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

  // dealloocate the atom list
  sfree( system->my_atoms, "system->my_atoms" );

  // deallocate the ffield parameters storage
  ff_params = &(system->reax_param);
  ntypes = ff_params->num_atom_types;

  sfree( ff_params->gp.l, "ff:globals" );

  for( i = 0; i < ntypes; ++i ) {
    for( j = 0; j < ntypes; ++j ) {
      for( k = 0; k < ntypes; ++k ) {
        sfree( ff_params->fbp[i][j][k], "ff:fbp[i,j,k]" );
      }
      sfree( ff_params->fbp[i][j], "ff:fbp[i,j]" );
      sfree( ff_params->thbp[i][j], "ff:thbp[i,j]" );
      sfree( ff_params->hbp[i][j], "ff:hbp[i,j]" );
    }
    sfree( ff_params->fbp[i], "ff:fbp[i]" );
    sfree( ff_params->thbp[i], "ff:thbp[i]" );
    sfree( ff_params->hbp[i], "ff:hbp[i]" );
    sfree( ff_params->tbp[i], "ff:tbp[i]" );
  }
  sfree( ff_params->fbp, "ff:fbp" );
  sfree( ff_params->thbp, "ff:thbp" );
  sfree( ff_params->hbp, "ff:hbp" );
  sfree( ff_params->tbp, "ff:tbp" );
  sfree( ff_params->sbp, "ff:sbp" );
}


/*************       workspace        *************/
void DeAllocate_Workspace( control_params *control, storage *workspace )
{
  int i;

  if( !workspace->allocated )
    return;

  workspace->allocated = 0;

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    sfree( workspace->tmp_dbl[i], "tmp_dbl[i]" );
    sfree( workspace->tmp_rvec[i], "tmp_rvec[i]" );
    sfree( workspace->tmp_rvec2[i], "tmp_rvec2[i]" );
  }

  /* bond order storage */
  sfree( workspace->within_bond_box, "skin" );
  sfree( workspace->total_bond_order, "total_bo" );
  sfree( workspace->Deltap, "Deltap" );
  sfree( workspace->Deltap_boc, "Deltap_boc" );
  sfree( workspace->dDeltap_self, "dDeltap_self" );
  sfree( workspace->Delta, "Delta" );
  sfree( workspace->Delta_lp, "Delta_lp" );
  sfree( workspace->Delta_lp_temp, "Delta_lp_temp" );
  sfree( workspace->dDelta_lp, "dDelta_lp" );
  sfree( workspace->dDelta_lp_temp, "dDelta_lp_temp" );
  sfree( workspace->Delta_e, "Delta_e" );
  sfree( workspace->Delta_boc, "Delta_boc" );
  sfree( workspace->Delta_val, "Delta_val" );
  sfree( workspace->nlp, "nlp" );
  sfree( workspace->nlp_temp, "nlp_temp" );
  sfree( workspace->Clp, "Clp" );
  sfree( workspace->vlpex, "vlpex" );
  sfree( workspace->bond_mark, "bond_mark" );
  sfree( workspace->done_after, "done_after" );

  /* QEq storage */
  sfree( workspace->Hdia_inv, "Hdia_inv" );
  sfree( workspace->b_s, "b_s" );
  sfree( workspace->b_t, "b_t" );
  sfree( workspace->b_prc, "b_prc" );
  sfree( workspace->b_prm, "b_prm" );
  sfree( workspace->s, "s" );
  sfree( workspace->t, "t" );
  sfree( workspace->droptol, "droptol" );
  sfree( workspace->b, "b" );
  sfree( workspace->x, "x" );

  /* GMRES storage */
  for( i = 0; i < RESTART+1; ++i ) {
    sfree( workspace->h[i], "h[i]" );
    sfree( workspace->v[i], "v[i]" );
  }
  sfree( workspace->h, "h" );
  sfree( workspace->v, "v" );
  sfree( workspace->y, "y" );
  sfree( workspace->z, "z" );
  sfree( workspace->g, "g" );
  sfree( workspace->hs, "hs" );
  sfree( workspace->hc, "hc" );
  /* CG storage */
  sfree( workspace->r, "r" );
  sfree( workspace->d, "d" );
  sfree( workspace->q, "q" );
  sfree( workspace->p, "p" );
  sfree( workspace->r2, "r2" );
  sfree( workspace->d2, "d2" );
  sfree( workspace->q2, "q2" );
  sfree( workspace->p2, "p2" );

  /* integrator */
  // sfree( workspace->f_old );
  sfree( workspace->v_const, "v_const" );

  /*workspace->realloc.num_far = -1;
    workspace->realloc.Htop = -1;
    workspace->realloc.hbonds = -1;
    workspace->realloc.bonds = -1;
    workspace->realloc.num_3body = -1;
    workspace->realloc.gcell_atoms = -1;*/

  /* storage for analysis */
  // if( control->molecular_analysis || control->diffusion_coef ) {
  //   sfree( workspace->mark, "mark" );
  //   sfree( workspace->old_mark, "old_mark" );
  // }

  // if( control->diffusion_coef )
  //   sfree( workspace->x_old, "x_old" );

  /* force related storage */
  sfree( workspace->f, "f" );
  sfree( workspace->CdDelta, "CdDelta" );
#ifdef TEST_FORCES
  sfree(workspace->dDelta, "dDelta" );
  sfree( workspace->f_ele, "f_ele" );
  sfree( workspace->f_vdw, "f_vdw" );
  sfree( workspace->f_bo, "f_bo" );
  sfree( workspace->f_be, "f_be" );
  sfree( workspace->f_lp, "f_lp" );
  sfree( workspace->f_ov, "f_ov" );
  sfree( workspace->f_un, "f_un" );
  sfree( workspace->f_ang, "f_ang" );
  sfree( workspace->f_coa, "f_coa" );
  sfree( workspace->f_pen, "f_pen" );
  sfree( workspace->f_hb, "f_hb" );
  sfree( workspace->f_tor, "f_tor" );
  sfree( workspace->f_con, "f_con" );
  sfree( workspace->f_tot, "f_tot" );

  sfree( workspace->rcounts, "rcounts" );
  sfree( workspace->displs, "displs" );
  sfree( workspace->id_all, "id_all" );
  sfree( workspace->f_all, "f_all" );
#endif

  /* hbond storage */
  //sfree( workspace->Hindex, "Hindex" );
  //sfree( workspace->num_bonds );
  //sfree( workspace->num_hbonds );
  //sfree( workspace->hash, "hash" );
  //sfree( workspace->rev_hash, "rev_hash" );
}


int Allocate_Workspace( reax_system *system, control_params *control,
                        storage *workspace, int local_cap, int total_cap,
                        MPI_Comm comm, char *msg )
{
  int i, total_real, total_rvec, local_int, local_real, local_rvec;

  workspace->allocated = 1;
  total_real = total_cap * sizeof(real);
  total_rvec = total_cap * sizeof(rvec);
  local_int = local_cap * sizeof(int);
  local_real = local_cap * sizeof(real);
  local_rvec = local_cap * sizeof(rvec);

  /* communication storage */
  for( i = 0; i < MAX_NBRS; ++i ) {
    workspace->tmp_dbl[i] = (real*)
      scalloc( total_cap, sizeof(real), "tmp_dbl", comm );
    workspace->tmp_rvec[i] = (rvec*)
      scalloc( total_cap, sizeof(rvec), "tmp_rvec", comm );
    workspace->tmp_rvec2[i] = (rvec2*)
      scalloc( total_cap, sizeof(rvec2), "tmp_rvec2", comm );
  }

  /* bond order related storage  */
  workspace->within_bond_box = (int*)
    scalloc( total_cap, sizeof(int), "skin", comm );
  workspace->total_bond_order = (real*) smalloc( total_real, "total_bo", comm );
  workspace->Deltap = (real*) smalloc( total_real, "Deltap", comm );
  workspace->Deltap_boc = (real*) smalloc( total_real, "Deltap_boc", comm );
  workspace->dDeltap_self = (rvec*) smalloc( total_rvec, "dDeltap_self", comm );
  workspace->Delta = (real*) smalloc( total_real, "Delta", comm );
  workspace->Delta_lp = (real*) smalloc( total_real, "Delta_lp", comm );
  workspace->Delta_lp_temp = (real*)
    smalloc( total_real, "Delta_lp_temp", comm );
  workspace->dDelta_lp = (real*) smalloc( total_real, "dDelta_lp", comm );
  workspace->dDelta_lp_temp = (real*)
    smalloc( total_real, "dDelta_lp_temp", comm );
  workspace->Delta_e = (real*) smalloc( total_real, "Delta_e", comm );
  workspace->Delta_boc = (real*) smalloc( total_real, "Delta_boc", comm );
  workspace->Delta_val = (real*) smalloc( total_real, "Delta_val", comm );
  workspace->nlp = (real*) smalloc( total_real, "nlp", comm );
  workspace->nlp_temp = (real*) smalloc( total_real, "nlp_temp", comm );
  workspace->Clp = (real*) smalloc( total_real, "Clp", comm );
  workspace->vlpex = (real*) smalloc( total_real, "vlpex", comm );
  workspace->bond_mark = (int*)
    scalloc( total_cap, sizeof(int), "bond_mark", comm );
  workspace->done_after = (int*)
    scalloc( total_cap, sizeof(int), "done_after", comm );
  // fprintf( stderr, "p%d: bond order storage\n", system->my_rank );

  /* QEq storage */
  workspace->Hdia_inv = (real*)
    scalloc( total_cap, sizeof(real), "Hdia_inv", comm );
  workspace->b_s = (real*) scalloc( total_cap, sizeof(real), "b_s", comm );
  workspace->b_t = (real*) scalloc( total_cap, sizeof(real), "b_t", comm );
  workspace->b_prc = (real*) scalloc( total_cap, sizeof(real), "b_prc", comm );
  workspace->b_prm = (real*) scalloc( total_cap, sizeof(real), "b_prm", comm );
  workspace->s = (real*) scalloc( total_cap, sizeof(real), "s", comm );
  workspace->t = (real*) scalloc( total_cap, sizeof(real), "t", comm );
  workspace->droptol = (real*)
    scalloc( total_cap, sizeof(real), "droptol", comm );
  workspace->b = (rvec2*) scalloc( total_cap, sizeof(rvec2), "b", comm );
  workspace->x = (rvec2*) scalloc( total_cap, sizeof(rvec2), "x", comm );

  /* GMRES storage */
  workspace->y = (real*) scalloc( RESTART+1, sizeof(real), "y", comm );
  workspace->z = (real*) scalloc( RESTART+1, sizeof(real), "z", comm );
  workspace->g = (real*) scalloc( RESTART+1, sizeof(real), "g", comm );
  workspace->h = (real**) scalloc( RESTART+1, sizeof(real*), "h", comm );
  workspace->hs = (real*) scalloc( RESTART+1, sizeof(real), "hs", comm );
  workspace->hc = (real*) scalloc( RESTART+1, sizeof(real), "hc", comm );
  workspace->v = (real**) scalloc( RESTART+1, sizeof(real*), "v", comm );

  for( i = 0; i < RESTART+1; ++i ) {
    workspace->h[i] = (real*) scalloc( RESTART+1, sizeof(real), "h[i]", comm );
    workspace->v[i] = (real*) scalloc( total_cap, sizeof(real), "v[i]", comm );
  }

  /* CG storage */
  workspace->r = (real*) scalloc( total_cap, sizeof(real), "r", comm );
  workspace->d = (real*) scalloc( total_cap, sizeof(real), "d", comm );
  workspace->q = (real*) scalloc( total_cap, sizeof(real), "q", comm );
  workspace->p = (real*) scalloc( total_cap, sizeof(real), "p", comm );
  workspace->r2 = (rvec2*) scalloc( total_cap, sizeof(rvec2), "r2", comm );
  workspace->d2 = (rvec2*) scalloc( total_cap, sizeof(rvec2), "d2", comm );
  workspace->q2 = (rvec2*) scalloc( total_cap, sizeof(rvec2), "q2", comm );
  workspace->p2 = (rvec2*) scalloc( total_cap, sizeof(rvec2), "p2", comm );

  /* integrator storage */
  workspace->v_const = (rvec*) smalloc( local_rvec, "v_const", comm );

  /* storage for analysis */
  // not yet implemented in the parallel version!!!
  // if( control->molecular_analysis || control->diffusion_coef ) {
  //   workspace->mark = (int*) scalloc( local_cap, sizeof(int), "mark", comm );
  //   workspace->old_mark = (int*)
  //     scalloc( local_cap, sizeof(int), "old_mark", comm );
  // }
  // else
  //   workspace->mark = workspace->old_mark = NULL;

  // if( control->diffusion_coef )
  //   workspace->x_old = (rvec*)
  //     scalloc( local_cap, sizeof(rvec), "x_old", comm );
  // else workspace->x_old = NULL;

  // /* force related storage */
  workspace->f = (rvec*) scalloc( total_cap, sizeof(rvec), "f", comm );
  workspace->CdDelta = (real*)
    scalloc( total_cap, sizeof(real), "CdDelta", comm );

#ifdef TEST_FORCES
  workspace->dDelta=(rvec*) smalloc( total_rvec, "dDelta", comm );
  workspace->f_ele =(rvec*) smalloc( total_rvec, "f_ele", comm );
  workspace->f_vdw =(rvec*) smalloc( total_rvec, "f_vdw", comm );
  workspace->f_bo =(rvec*) smalloc( total_rvec, "f_bo", comm );
  workspace->f_be =(rvec*) smalloc( total_rvec, "f_be", comm );
  workspace->f_lp =(rvec*) smalloc( total_rvec, "f_lp", comm );
  workspace->f_ov =(rvec*) smalloc( total_rvec, "f_ov", comm );
  workspace->f_un =(rvec*) smalloc( total_rvec, "f_un", comm );
  workspace->f_ang =(rvec*) smalloc( total_rvec, "f_ang", comm );
  workspace->f_coa =(rvec*) smalloc( total_rvec, "f_coa", comm );
  workspace->f_pen =(rvec*) smalloc( total_rvec, "f_pen", comm );
  workspace->f_hb =(rvec*) smalloc( total_rvec, "f_hb", comm );
  workspace->f_tor =(rvec*) smalloc( total_rvec, "f_tor", comm );
  workspace->f_con =(rvec*) smalloc( total_rvec, "f_con", comm );
  workspace->f_tot =(rvec*) smalloc( total_rvec, "f_tot", comm );

  if( system->my_rank == MASTER_NODE ) {
    workspace->rcounts = (int*)
      smalloc( system->wsize*sizeof(int), "rcount", comm );
    workspace->displs = (int*)
      smalloc( system->wsize*sizeof(int), "displs", comm );
    workspace->id_all = (int*)
      smalloc( system->bigN*sizeof(int), "id_all", comm );
    workspace->f_all = (rvec*)
      smalloc( system->bigN*sizeof(rvec), "f_all", comm );
  }
  else{
    workspace->rcounts = NULL;
    workspace->displs = NULL;
    workspace->id_all = NULL;
    workspace->f_all = NULL;
  }
#endif

  return SUCCESS;
}


void Reallocate_Neighbor_List( reax_list *far_nbrs, int n, int num_intrs,
                               MPI_Comm comm )
{
  Delete_List( far_nbrs, comm );
  if(!Make_List( n, num_intrs, TYP_FAR_NEIGHBOR, far_nbrs, comm )){
    fprintf(stderr, "Problem in initializing far nbrs list. Terminating!\n");
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }
}


int Allocate_Matrix( sparse_matrix **pH, int cap, int m, MPI_Comm comm )
{
  sparse_matrix *H;

  *pH = (sparse_matrix*)
    smalloc( sizeof(sparse_matrix), "sparse_matrix", comm );
  H = *pH;
  H->cap = cap;
  H->m = m;
  H->start = (int*) smalloc( sizeof(int) * cap, "matrix_start", comm );
  H->end = (int*) smalloc( sizeof(int) * cap, "matrix_end", comm );
  H->entries = (sparse_matrix_entry*)
    smalloc( sizeof(sparse_matrix_entry)*m, "matrix_entries", comm );

  return SUCCESS;
}


void Deallocate_Matrix( sparse_matrix *H )
{
  sfree(H->start, "H->start");
  sfree(H->end, "H->end");
  sfree(H->entries, "H->entries");
  sfree(H, "H");
}


int Reallocate_Matrix( sparse_matrix **H, int n, int m, char *name,
                       MPI_Comm comm )
{
  Deallocate_Matrix( *H );
  if( !Allocate_Matrix( H, n, m, comm ) ) {
    fprintf(stderr, "not enough space for %s matrix. terminating!\n", name);
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "reallocating %s matrix, n = %d, m = %d\n", name, n, m );
  fprintf( stderr, "memory allocated: %s = %dMB\n",
           name, (int)(m * sizeof(sparse_matrix_entry) / (1024*1024)) );
#endif
  return SUCCESS;
}


int Reallocate_HBonds_List( reax_system *system, reax_list *hbonds,
                            MPI_Comm comm )
{
  int i, id, total_hbonds;

  int mincap = system->mincap;
  double saferzone = system->saferzone;

  total_hbonds = 0;
  for( i = 0; i < system->n; ++i )
    if( (id = system->my_atoms[i].Hindex) >= 0 ) {
      // commented out - already updated in validate_lists in forces.c
      // system->my_atoms[i].num_hbonds = MAX(Num_Entries(id,hbonds)*SAFER_ZONE,
      //                                   MIN_HBONDS);
      total_hbonds += system->my_atoms[i].num_hbonds;
    }
  total_hbonds = (int)(MAX( total_hbonds*saferzone, mincap*MIN_HBONDS ));

  Delete_List( hbonds, comm );
  if( !Make_List( system->Hcap, total_hbonds, TYP_HBOND, hbonds, comm ) ) {
    fprintf( stderr, "not enough space for hbonds list. terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return total_hbonds;
}


int Reallocate_Bonds_List( reax_system *system, reax_list *bonds,
                           int *total_bonds, int *est_3body, MPI_Comm comm )
{
  int i;

  int mincap = system->mincap;
  double safezone = system->safezone;

  *total_bonds = 0;
  *est_3body = 0;
  for( i = 0; i < system->N; ++i ){
    *est_3body += SQR(system->my_atoms[i].num_bonds);
    // commented out - already updated in validate_lists in forces.c
    // system->my_atoms[i].num_bonds = MAX( Num_Entries(i,bonds)*2, MIN_BONDS );
    *total_bonds += system->my_atoms[i].num_bonds;
  }
  *total_bonds = (int)(MAX( *total_bonds * safezone, mincap*MIN_BONDS ));

  Delete_List( bonds, comm );
  if(!Make_List(system->total_cap, *total_bonds, TYP_BOND, bonds, comm)) {
    fprintf( stderr, "not enough space for bonds list. terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  return SUCCESS;
}


/*************       grid        *************/
int Estimate_GCell_Population( reax_system* system, MPI_Comm comm )
{
  int d, i, j, k, l, max_atoms, my_max, all_max;
  ivec c;
  grid *g;
  grid_cell *gc;
  simulation_box *big_box, *my_ext_box;
  reax_atom *atoms;

  big_box    = &(system->big_box);
  my_ext_box = &(system->my_ext_box);
  g          = &(system->my_grid);
  atoms      = system->my_atoms;
  Reset_Grid( g );

  for( l = 0; l < system->n; l++ ) {
    for( d = 0; d < 3; ++d ) {
      //if( atoms[l].x[d] < big_box->min[d] )
      //  atoms[l].x[d] += big_box->box_norms[d];
      //else if( atoms[l].x[d] >= big_box->max[d] )
      //  atoms[l].x[d] -= big_box->box_norms[d];

      c[d] = (int)((atoms[l].x[d]-my_ext_box->min[d])*g->inv_len[d]);

      if( c[d] >= g->native_end[d] )
        c[d] = g->native_end[d] - 1;
      else if( c[d] < g->native_str[d] )
        c[d] = g->native_str[d];
    }
#if defined(DEBUG)
    fprintf( stderr, "p%d bin_my_atoms: l:%d - atom%d @ %.5f %.5f %.5f"        \
             "--> cell: %d %d %d\n",
             system->my_rank, l, atoms[l].orig_id,
             atoms[l].x[0], atoms[l].x[1], atoms[l].x[2],
             c[0], c[1], c[2] );
#endif
    gc = &( g->cells[c[0]][c[1]][c[2]] );
    gc->top++;
  }

  max_atoms = 0;
  for( i = 0; i < g->ncells[0]; i++ )
    for( j = 0; j < g->ncells[1]; j++ )
      for( k = 0; k < g->ncells[2]; k++ ) {
        gc = &(g->cells[i][j][k]);
        if( max_atoms < gc->top )
          max_atoms = gc->top;
#if defined(DEBUG)
        fprintf( stderr, "p%d gc[%d,%d,%d]->top=%d\n",
                 system->my_rank, i, j, k, gc->top );
#endif
      }

  my_max = (int)(MAX(max_atoms*SAFE_ZONE, MIN_GCELL_POPL));
  MPI_Allreduce( &my_max, &all_max, 1, MPI_INT, MPI_MAX, comm );
#if defined(DEBUG)
  fprintf( stderr, "p%d max_atoms=%d, my_max=%d, all_max=%d\n",
           system->my_rank, max_atoms, my_max, all_max );
#endif

  return all_max;
}


void Allocate_Grid( reax_system *system, MPI_Comm comm )
{
  int i, j, k, l;
  grid *g;
  grid_cell *gc;

  g = &( system->my_grid );

  /* allocate gcell reordering space */
  g->order = (ivec*) scalloc( g->total+1, sizeof(ivec), "g:order", comm );

  /* allocate the gcells for the new grid */
  g->max_nbrs = (2*g->vlist_span[0]+1)*(2*g->vlist_span[1]+1)*
    (2*g->vlist_span[2]+1)+3;

  g->cells = (grid_cell***)
    scalloc( g->ncells[0], sizeof(grid_cell**), "gcells", comm );
  for( i = 0; i < g->ncells[0]; i++ ) {
    g->cells[i] = (grid_cell**)
      scalloc( g->ncells[1], sizeof(grid_cell*),"gcells[i]", comm );

    for( j = 0; j < g->ncells[1]; ++j )        {
      g->cells[i][j] = (grid_cell*)
        scalloc( g->ncells[2], sizeof(grid_cell), "gcells[i][j]", comm );

      for( k = 0; k < g->ncells[2]; k++ ) {
        gc = &(g->cells[i][j][k]);
        gc->top = gc->mark = gc->str = gc->end = 0;
        gc->nbrs = (grid_cell**)
          scalloc( g->max_nbrs, sizeof(grid_cell*), "g:nbrs", comm );
        gc->nbrs_x = (ivec*)
          scalloc( g->max_nbrs, sizeof(ivec), "g:nbrs_x", comm );
        gc->nbrs_cp = (rvec*)
          scalloc( g->max_nbrs, sizeof(rvec), "g:nbrs_cp", comm );
        for( l = 0; l < g->max_nbrs; ++l )
          gc->nbrs[l] = NULL;
      }
    }
  }

  /* allocate atom id storage in gcells */
  g->max_atoms = Estimate_GCell_Population( system, comm );
  /* space for storing atom id's is required only for native cells */
  for( i = g->native_str[0]; i < g->native_end[0]; ++i )
    for( j = g->native_str[1]; j < g->native_end[1]; ++j )
      for( k = g->native_str[2]; k < g->native_end[2]; ++k )
        g->cells[i][j][k].atoms = (int*) scalloc( g->max_atoms, sizeof(int),
                                                  "g:atoms", comm );
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d-allocated %dx%dx%d grid: nbrs=%d atoms=%d space=%dMB\n",
           system->my_rank, g->ncells[0], g->ncells[1], g->ncells[2],
           g->max_nbrs, g->max_atoms,
           (int)
           ((g->total*sizeof(grid_cell)+g->total*g->max_nbrs*sizeof(int*) +
             g->total*g->max_nbrs*sizeof(rvec) +
             (g->native_end[0]-g->native_str[0])*
             (g->native_end[1]-g->native_str[1])*
             (g->native_end[2]-g->native_str[2])*g->max_atoms*sizeof(int))/
            (1024*1024)) );
#endif
}


void Deallocate_Grid( grid *g )
{
  int i, j, k;
  grid_cell *gc;

  sfree( g->order, "g:order" );

  /* deallocate the grid cells */
  for( i = 0; i < g->ncells[0]; i++ ) {
    for( j = 0; j < g->ncells[1]; j++ ) {
      for( k = 0; k < g->ncells[2]; k++ ) {
        gc = &(g->cells[i][j][k]);
        sfree( gc->nbrs, "g:nbrs" );
        sfree( gc->nbrs_x, "g:nbrs_x" );
        sfree( gc->nbrs_cp, "g:nbrs_cp" );
        if(gc->atoms != NULL )
          sfree( gc->atoms, "g:atoms" );
      }
      sfree( g->cells[i][j], "g:cells[i][j]" );
    }
    sfree( g->cells[i], "g:cells[i]" );
  }
  sfree( g->cells, "g:cells" );
}


/************ mpi buffers ********************/
 /* make the allocations based on the largest need by the three comms:
    1- transfer an atom who has moved into other proc's domain (mpi_atom)
    2- exchange boundary atoms (boundary_atom)
    3- update position info for boundary atoms (mpi_rvec)

    the largest space by far is required for the 2nd comm operation above.
    buffers are void*, type cast to the correct pointer type to access
    the allocated buffers */
int  Allocate_MPI_Buffers( mpi_datatypes *mpi_data, int est_recv,
                           neighbor_proc *my_nbrs, char *msg )
{
  int i;
  mpi_out_data  *mpi_buf;
  MPI_Comm comm;

  comm = mpi_data->world;

  /* in buffers */
  mpi_data->in1_buffer = (void*)
    scalloc( est_recv, sizeof(boundary_atom), "in1_buffer", comm );
  mpi_data->in2_buffer = (void*)
    scalloc( est_recv, sizeof(boundary_atom), "in2_buffer", comm );

  /* out buffers */
  for( i = 0; i < MAX_NBRS; ++i ) {
    mpi_buf = &( mpi_data->out_buffers[i] );
    /* allocate storage for the neighbor processor i */
    mpi_buf->index = (int*)
      scalloc( my_nbrs[i].est_send, sizeof(int), "mpibuf:index", comm );
    mpi_buf->out_atoms = (void*)
      scalloc( my_nbrs[i].est_send, sizeof(boundary_atom), "mpibuf:out_atoms",
               comm );
  }

  return SUCCESS;
}


void Deallocate_MPI_Buffers( mpi_datatypes *mpi_data )
{
  int i;
  mpi_out_data  *mpi_buf;

  sfree( mpi_data->in1_buffer, "in1_buffer" );
  sfree( mpi_data->in2_buffer, "in2_buffer" );

  for( i = 0; i < MAX_NBRS; ++i ) {
    mpi_buf = &( mpi_data->out_buffers[i] );
    sfree( mpi_buf->index, "mpibuf:index" );
    sfree( mpi_buf->out_atoms, "mpibuf:out_atoms" );
  }
}


void ReAllocate( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists,
                 mpi_datatypes *mpi_data )
{
  int num_bonds, est_3body, nflag, Nflag, Hflag, ret;
  int renbr, newsize;
  reallocate_data *realloc;
  reax_list *far_nbrs;
  grid *g;
  MPI_Comm comm;
  char msg[200];

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  realloc = &(workspace->realloc);
  g = &(system->my_grid);
  comm = mpi_data->world;

#if defined(DEBUG)
  fprintf( stderr, "p%d@reallocate: n: %d, N: %d, numH: %d\n",
           system->my_rank, system->n, system->N, system->numH );
  fprintf( stderr, "p%d@reallocate: local_cap: %d, total_cap: %d, Hcap: %d\n",
           system->my_rank, system->local_cap, system->total_cap,
           system->Hcap);
  fprintf( stderr, "p%d: realloc.num_far: %d\n",
           system->my_rank, realloc->num_far );
  fprintf( stderr, "p%d: realloc.H: %d, realloc.Htop: %d\n",
           system->my_rank, realloc->H, realloc->Htop );
  fprintf( stderr, "p%d: realloc.Hbonds: %d, realloc.num_hbonds: %d\n",
           system->my_rank, realloc->hbonds, realloc->num_hbonds );
  fprintf( stderr, "p%d: realloc.bonds: %d, num_bonds: %d\n",
           system->my_rank, realloc->bonds, realloc->num_bonds );
  fprintf( stderr, "p%d: realloc.num_3body: %d\n",
           system->my_rank, realloc->num_3body );
#endif

  // IMPORTANT: LOOSE ZONES CHECKS ARE DISABLED FOR NOW BY &&'ing with 0!!!
  nflag = 0;
  if( system->n >= DANGER_ZONE * system->local_cap ||
      (0 && system->n <= LOOSE_ZONE * system->local_cap) ) {
    nflag = 1;
    system->local_cap = MAX( (int)(system->n * safezone), mincap );
  }

  Nflag = 0;
  if( system->N >= DANGER_ZONE * system->total_cap ||
      (0 && system->N <= LOOSE_ZONE * system->total_cap) ) {
    Nflag = 1;
    system->total_cap = MAX( (int)(system->N * safezone), mincap );
  }

  if( Nflag ) {
    /* system */
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating system and workspace -"\
             "n=%d  N=%d  local_cap=%d  total_cap=%d\n",
             system->my_rank, system->n, system->N,
             system->local_cap, system->total_cap );
#endif
    ret = Allocate_System( system, system->local_cap, system->total_cap, msg );
    if( ret != SUCCESS ) {
      fprintf( stderr, "not enough space for atom_list: total_cap=%d",
               system->total_cap );
      fprintf( stderr, "terminating...\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }

    /* workspace */
    DeAllocate_Workspace( control, workspace );
    ret = Allocate_Workspace( system, control, workspace, system->local_cap,
                              system->total_cap, comm, msg );
    if( ret != SUCCESS ) {
      fprintf( stderr, "no space for workspace: local_cap=%d total_cap=%d",
               system->local_cap, system->total_cap );
      fprintf( stderr, "terminating...\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
  }


  renbr = (data->step - data->prev_steps) % control->reneighbor == 0;
  /* far neighbors */
  if( renbr ) {
    far_nbrs = *lists + FAR_NBRS;

    if( Nflag || realloc->num_far >= far_nbrs->num_intrs * DANGER_ZONE ) {
      if( realloc->num_far > far_nbrs->num_intrs ) {
        fprintf( stderr, "step%d-ran out of space on far_nbrs: top=%d, max=%d",
                 data->step, realloc->num_far, far_nbrs->num_intrs );
        MPI_Abort( comm, INSUFFICIENT_MEMORY );
      }

      newsize = static_cast<int>
        (MAX( realloc->num_far*safezone, mincap*MIN_NBRS ));
#if defined(DEBUG_FOCUS)
      fprintf( stderr, "p%d: reallocating far_nbrs: num_fars=%d, space=%dMB\n",
               system->my_rank, (int)(realloc->num_far*SAFE_ZONE),
               (newsize*sizeof(far_neighbor_data)/(1024*1024)) );
#endif

      Reallocate_Neighbor_List( far_nbrs, system->total_cap, newsize, comm );
      realloc->num_far = 0;
    }
  }

#if defined(PURE_REAX)
  /* qeq coef matrix */
  H = workspace->H;
  if( nflag || realloc->Htop >= H->m * DANGER_ZONE ) {
    if( realloc->Htop > H->m ) {
      fprintf( stderr,
               "step%d - ran out of space on H matrix: Htop=%d, max = %d",
               data->step, realloc->Htop, H->m );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating H matrix: Htop=%d, space=%dMB\n",
             system->my_rank, (int)(realloc->Htop*SAFE_ZONE),
             (int)(realloc->Htop * SAFE_ZONE * sizeof(sparse_matrix_entry) /
                   (1024*1024)) );
#endif

    newsize = static_cast<int>
        (MAX( realloc->Htop*safezone, mincap*MIN_NBRS ));
    Reallocate_Matrix( &(workspace->H), system->local_cap,
                       newsize, "H", comm );
    //Deallocate_Matrix( workspace->L );
    //Deallocate_Matrix( workspace->U );
    workspace->L = NULL;
    workspace->U = NULL;
    realloc->Htop = 0;
  }
#endif /*PURE_REAX*/

  /* hydrogen bonds list */
  if( control->hbond_cut > 0 ) {
    Hflag = 0;
    if( system->numH >= DANGER_ZONE * system->Hcap ||
        (0 && system->numH <= LOOSE_ZONE * system->Hcap) ) {
      Hflag = 1;
      system->Hcap = int(MAX( system->numH * saferzone, mincap ));
    }

    if( Hflag || realloc->hbonds ) {
      ret = Reallocate_HBonds_List( system, (*lists)+HBONDS, comm );
      realloc->hbonds = 0;
#if defined(DEBUG_FOCUS)
      fprintf(stderr, "p%d: reallocating hbonds: total_hbonds=%d space=%dMB\n",
              system->my_rank, ret, (int)(ret*sizeof(hbond_data)/(1024*1024)));
#endif
    }
  }

  /* bonds list */
  num_bonds = est_3body = -1;
  if( Nflag || realloc->bonds ){
    Reallocate_Bonds_List( system, (*lists)+BONDS, &num_bonds,
                           &est_3body, comm );
    realloc->bonds = 0;
    realloc->num_3body = MAX( realloc->num_3body, est_3body );
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating bonds: total_bonds=%d, space=%dMB\n",
             system->my_rank, num_bonds,
             (int)(num_bonds*sizeof(bond_data)/(1024*1024)) );
#endif
  }

  /* 3-body list */
  if( realloc->num_3body > 0 ) {
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating 3body list: num_3body=%d, space=%dMB\n",
             system->my_rank, realloc->num_3body,
             (int)(realloc->num_3body * sizeof(three_body_interaction_data) /
                   (1024*1024)) );
#endif
    Delete_List( (*lists)+THREE_BODIES, comm );

    if( num_bonds == -1 )
      num_bonds = ((*lists)+BONDS)->num_intrs;

    realloc->num_3body = (int)(MAX(realloc->num_3body*safezone, MIN_3BODIES));

    if( !Make_List( num_bonds, realloc->num_3body, TYP_THREE_BODY,
                    (*lists)+THREE_BODIES, comm ) ) {
      fprintf( stderr, "Problem in initializing angles list. Terminating!\n" );
      MPI_Abort( comm, CANNOT_INITIALIZE );
    }
    realloc->num_3body = -1;
  }

#if defined(PURE_REAX)
  /* grid */
  if( renbr && realloc->gcell_atoms > -1 ) {
#if defined(DEBUG_FOCUS)
    fprintf(stderr, "reallocating gcell: g->max_atoms: %d\n", g->max_atoms);
#endif
    for( i = g->native_str[0]; i < g->native_end[0]; i++ )
      for( j = g->native_str[1]; j < g->native_end[1]; j++ )
        for( k = g->native_str[2]; k < g->native_end[2]; k++ ) {
          // reallocate g->atoms
          sfree( g->cells[i][j][k].atoms, "g:atoms" );
          g->cells[i][j][k].atoms = (int*)
            scalloc( realloc->gcell_atoms, sizeof(int), "g:atoms", comm);
        }
    realloc->gcell_atoms = -1;
  }

  /* mpi buffers */
  // we have to be at a renbring step -
  // to ensure correct values at mpi_buffers for update_boundary_positions
  if( !renbr )
    mpi_flag = 0;
  // check whether in_buffer capacity is enough
  else if( system->max_recved >= system->est_recv * 0.90 )
    mpi_flag = 1;
  else {
    // otherwise check individual outgoing buffers
    mpi_flag = 0;
    for( p = 0; p < MAX_NBRS; ++p ) {
      nbr_pr   = &( system->my_nbrs[p] );
      nbr_data = &( mpi_data->out_buffers[p] );
      if( nbr_data->cnt >= nbr_pr->est_send * 0.90 ) {
        mpi_flag = 1;
        break;
      }
    }
  }

  if( mpi_flag ) {
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating mpi_buf: old_recv=%d\n",
             system->my_rank, system->est_recv );
    for( p = 0; p < MAX_NBRS; ++p )
      fprintf( stderr, "p%d: nbr%d old_send=%d\n",
               system->my_rank, p, system->my_nbrs[p].est_send );
#endif
    /* update mpi buffer estimates based on last comm */
    system->est_recv = MAX( system->max_recved*SAFER_ZONE, MIN_SEND );
    system->est_trans =
      (system->est_recv * sizeof(boundary_atom)) / sizeof(mpi_atom);
    total_send = 0;
    for( p = 0; p < MAX_NBRS; ++p ) {
      nbr_pr   = &( system->my_nbrs[p] );
      nbr_data = &( mpi_data->out_buffers[p] );
      nbr_pr->est_send = MAX( nbr_data->cnt*SAFER_ZONE, MIN_SEND );
      total_send += nbr_pr->est_send;
    }
#if defined(DEBUG_FOCUS)
    fprintf( stderr, "p%d: reallocating mpi_buf: recv=%d send=%d total=%dMB\n",
             system->my_rank, system->est_recv, total_send,
             (int)((system->est_recv+total_send)*sizeof(boundary_atom)/
                   (1024*1024)));
    for( p = 0; p < MAX_NBRS; ++p )
      fprintf( stderr, "p%d: nbr%d new_send=%d\n",
               system->my_rank, p, system->my_nbrs[p].est_send );
#endif

    /* reallocate mpi buffers */
    Deallocate_MPI_Buffers( mpi_data );
    ret = Allocate_MPI_Buffers( mpi_data, system->est_recv,
                                system->my_nbrs, msg );
    if( ret != SUCCESS ) {
      fprintf( stderr, "%s", msg );
      fprintf( stderr, "terminating...\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
  }
#endif /*PURE_REAX*/

#if defined(DEBUG_FOCUS)
  fprintf( stderr, "p%d @ step%d: reallocate done\n",
           system->my_rank, data->step );
  MPI_Barrier( comm );
#endif
}
