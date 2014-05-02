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
#include "reaxc_allocate.h"
#include "reaxc_list.h"
#include "reaxc_reset_tools.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

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

  system->my_atoms = (reax_atom*)
    scalloc( system->total_cap, sizeof(reax_atom), "my_atoms", comm );

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
  sfree( workspace->v_const, "v_const" );

  /* force related storage */
  sfree( workspace->f, "f" );
  sfree( workspace->CdDelta, "CdDelta" );

}


int Allocate_Workspace( reax_system *system, control_params *control,
                        storage *workspace, int local_cap, int total_cap,
                        MPI_Comm comm, char *msg )
{
  int i, total_real, total_rvec, local_rvec;

  workspace->allocated = 1;
  total_real = total_cap * sizeof(real);
  total_rvec = total_cap * sizeof(rvec);
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

  // /* force related storage */
  workspace->f = (rvec*) scalloc( total_cap, sizeof(rvec), "f", comm );
  workspace->CdDelta = (real*)
    scalloc( total_cap, sizeof(real), "CdDelta", comm );

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
  simulation_box *my_ext_box;
  reax_atom *atoms;

  my_ext_box = &(system->my_ext_box);
  g          = &(system->my_grid);
  atoms      = system->my_atoms;
  Reset_Grid( g );

  for( l = 0; l < system->n; l++ ) {
    for( d = 0; d < 3; ++d ) {

      c[d] = (int)((atoms[l].x[d]-my_ext_box->min[d])*g->inv_len[d]);

      if( c[d] >= g->native_end[d] )
        c[d] = g->native_end[d] - 1;
      else if( c[d] < g->native_str[d] )
        c[d] = g->native_str[d];
    }
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
      }

  my_max = (int)(MAX(max_atoms*SAFE_ZONE, MIN_GCELL_POPL));
  MPI_Allreduce( &my_max, &all_max, 1, MPI_INT, MPI_MAX, comm );

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
  int num_bonds, est_3body, Hflag, ret;
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

  int nflag = 0;
  if( system->n >= DANGER_ZONE * system->local_cap ||
      (0 && system->n <= LOOSE_ZONE * system->local_cap) ) {
    nflag = 1;
    system->local_cap = MAX( (int)(system->n * safezone), mincap );
  }

  int Nflag = 0;
  if( system->N >= DANGER_ZONE * system->total_cap ||
      (0 && system->N <= LOOSE_ZONE * system->total_cap) ) {
    Nflag = 1;
    system->total_cap = MAX( (int)(system->N * safezone), mincap );
  }

  if( Nflag ) {
    /* system */
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

      Reallocate_Neighbor_List( far_nbrs, system->total_cap, newsize, comm );
      realloc->num_far = 0;
    }
  }

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
    }
  }

  /* bonds list */
  num_bonds = est_3body = -1;
  if( Nflag || realloc->bonds ){
    Reallocate_Bonds_List( system, (*lists)+BONDS, &num_bonds,
                           &est_3body, comm );
    realloc->bonds = 0;
    realloc->num_3body = MAX( realloc->num_3body, est_3body );
  }

  /* 3-body list */
  if( realloc->num_3body > 0 ) {
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

}
