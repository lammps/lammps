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
#include "update.h"
#include "reaxc_io_tools.h"
#include "reaxc_basic_comm.h"
#include "reaxc_list.h"
#include "reaxc_reset_tools.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_traj.h"
#include "reaxc_vector.h"

print_interaction Print_Interactions[NUM_INTRS];

int Init_Output_Files( reax_system *system, control_params *control,
                       output_controls *out_control, mpi_datatypes *mpi_data,
                       char *msg )
{
  char temp[MAX_STR];
  int ret;

  if( out_control->write_steps > 0 ){
    ret = Init_Traj( system, control, out_control, mpi_data, msg );
    if( ret == FAILURE )
      return ret;
  }

  if( system->my_rank == MASTER_NODE ) {
    /* These files are written only by the master node */
    if( out_control->energy_update_freq > 0 ) {

      /* init potentials file */
      sprintf( temp, "%s.pot", control->sim_name );
      if( (out_control->pot = fopen( temp, "w" )) != NULL ) {
        fflush( out_control->pot );
      }
      else {
        strcpy( msg, "init_out_controls: .pot file could not be opened\n" );
        return FAILURE;
      }

      /* init log file */
    }

    /* init pressure file */
    if( control->ensemble == NPT  ||
        control->ensemble == iNPT ||
        control->ensemble == sNPT ) {
      sprintf( temp, "%s.prs", control->sim_name );
      if( (out_control->prs = fopen( temp, "w" )) != NULL ) {
        fprintf(out_control->prs,"%8s%13s%13s%13s%13s%13s%13s%13s\n",
                "step", "Pint/norm[x]", "Pint/norm[y]", "Pint/norm[z]",
                "Pext/Ptot[x]", "Pext/Ptot[y]", "Pext/Ptot[z]", "Pkin/V" );
        fflush( out_control->prs );
      }
      else {
        strcpy(msg,"init_out_controls: .prs file couldn't be opened\n");
        return FAILURE;
      }
    }
  }

  return SUCCESS;
}


/************************ close output files ************************/
int Close_Output_Files( reax_system *system, control_params *control,
                        output_controls *out_control, mpi_datatypes *mpi_data )
{
  if( out_control->write_steps > 0 )
    End_Traj( system->my_rank, out_control );

  if( system->my_rank == MASTER_NODE ) {
    if( out_control->energy_update_freq > 0 ) {
      fclose( out_control->pot );
    }

    if( control->ensemble == NPT || control->ensemble == iNPT ||
        control->ensemble == sNPT )
      fclose( out_control->prs );
  }

  return SUCCESS;
}



void Print_Box( simulation_box* box, char *name, FILE *out )
{
  // int i, j;

  fprintf( out, "%s:\n", name );
  fprintf( out, "\tmin[%8.3f %8.3f %8.3f]\n",
           box->min[0], box->min[1], box->min[2] );
  fprintf( out, "\tmax[%8.3f %8.3f %8.3f]\n",
           box->max[0], box->max[1], box->max[2] );
  fprintf( out, "\tdims[%8.3f%8.3f%8.3f]\n",
           box->box_norms[0], box->box_norms[1], box->box_norms[2] );

}



void Print_Grid( grid* g, FILE *out )
{
  int x, y, z, gc_type;
  ivec gc_str;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  fprintf( out, "\tnumber of grid cells: %d %d %d\n",
           g->ncells[0], g->ncells[1], g->ncells[2] );
  fprintf( out, "\tgcell lengths: %8.3f %8.3f %8.3f\n",
           g->cell_len[0], g->cell_len[1], g->cell_len[2] );
  fprintf( out, "\tinverses of gcell lengths: %8.3f %8.3f %8.3f\n",
           g->inv_len[0], g->inv_len[1], g->inv_len[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tnumber of native gcells: %d %d %d\n",
           g->native_cells[0], g->native_cells[1], g->native_cells[2] );
  fprintf( out, "\tnative gcell span: %d-%d  %d-%d  %d-%d\n",
           g->native_str[0], g->native_end[0],
           g->native_str[1], g->native_end[1],
           g->native_str[2], g->native_end[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tvlist gcell stretch: %d %d %d\n",
           g->vlist_span[0], g->vlist_span[1], g->vlist_span[2] );
  fprintf( out, "\tnonbonded nbrs gcell stretch: %d %d %d\n",
           g->nonb_span[0], g->nonb_span[1], g->nonb_span[2] );
  fprintf( out, "\tbonded nbrs gcell stretch: %d %d %d\n",
           g->bond_span[0], g->bond_span[1], g->bond_span[2] );
  fprintf( out, "\t---------------------------------\n" );
  fprintf( out, "\tghost gcell span: %d %d %d\n",
           g->ghost_span[0], g->ghost_span[1], g->ghost_span[2] );
  fprintf( out, "\tnonbonded ghost gcell span: %d %d %d\n",
           g->ghost_nonb_span[0],g->ghost_nonb_span[1],g->ghost_nonb_span[2]);
  fprintf(out, "\thbonded ghost gcell span: %d %d %d\n",
          g->ghost_hbond_span[0],g->ghost_hbond_span[1],g->ghost_hbond_span[2]);
  fprintf( out, "\tbonded ghost gcell span: %d %d %d\n",
           g->ghost_bond_span[0],g->ghost_bond_span[1],g->ghost_bond_span[2]);
  fprintf( out, "\t---------------------------------\n" );

  fprintf( stderr, "GCELL MARKS:\n" );
  gc_type = g->cells[0][0][0].type;
  ivec_MakeZero( gc_str );

  x = y = z = 0;
  for( x = 0; x < g->ncells[0]; ++x )
    for( y = 0; y < g->ncells[1]; ++y )
      for( z = 0; z < g->ncells[2]; ++z )
        if( g->cells[x][y][z].type != gc_type ){
          fprintf( stderr,
                   "\tgcells from(%2d %2d %2d) to (%2d %2d %2d): %d - %s\n",
                   gc_str[0], gc_str[1], gc_str[2], x, y, z,
                   gc_type, gcell_type_text[gc_type] );
          gc_type = g->cells[x][y][z].type;
          gc_str[0] = x;
          gc_str[1] = y;
          gc_str[2] = z;
        }
  fprintf( stderr, "\tgcells from(%2d %2d %2d) to (%2d %2d %2d): %d - %s\n",
           gc_str[0], gc_str[1], gc_str[2], x, y, z,
           gc_type, gcell_type_text[gc_type] );
  fprintf( out, "-------------------------------------\n" );
}

void Print_Native_GCells( reax_system *system )
{
  int        i, j, k, l;
  char       fname[100];
  FILE      *f;
  grid      *g;
  grid_cell *gc;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  sprintf( fname, "native_gcells.%d", system->my_rank );
  f = fopen( fname, "w" );
  g = &(system->my_grid);

  for( i = g->native_str[0]; i < g->native_end[0]; i++ )
    for( j = g->native_str[1]; j < g->native_end[1]; j++ )
      for( k = g->native_str[2]; k < g->native_end[2]; k++ )
          {
            gc = &( g->cells[i][j][k] );

            fprintf( f, "p%d gcell(%2d %2d %2d) of type %d(%s)\n",
                     system->my_rank, i, j, k,
                   gc->type, gcell_type_text[gc->type] );

            fprintf( f, "\tatom list start: %d, end: %d\n\t", gc->str, gc->end );

          for( l = gc->str; l < gc->end; ++l )
              fprintf( f, TAGINT_FORMAT, system->my_atoms[l].orig_id );
            fprintf( f, "\n" );
          }

  fclose(f);
}



void Print_All_GCells( reax_system *system )
{
  int        i, j, k, l;
  char       fname[100];
  FILE      *f;
  grid      *g;
  grid_cell *gc;
  char gcell_type_text[10][12] =
    { "NO_NBRS", "NEAR_ONLY", "HBOND_ONLY", "FAR_ONLY",
      "NEAR_HBOND", "NEAR_FAR", "HBOND_FAR", "FULL_NBRS", "NATIVE" };

  sprintf( fname, "all_gcells.%d", system->my_rank );
  f = fopen( fname, "w" );
  g = &(system->my_grid);

  for( i = 0; i < g->ncells[0]; i++ )
    for( j = 0; j < g->ncells[1]; j++ )
      for( k = 0; k < g->ncells[2]; k++ )
          {
            gc = &( g->cells[i][j][k] );

            fprintf( f, "p%d gcell(%2d %2d %2d) of type %d(%s)\n",
                     system->my_rank, i, j, k,
                   gc->type, gcell_type_text[gc->type] );

            fprintf( f, "\tatom list start: %d, end: %d\n\t", gc->str, gc->end );

          for( l = gc->str; l < gc->end; ++l )
              fprintf( f, TAGINT_FORMAT, system->my_atoms[l].orig_id );
            fprintf( f, "\n" );
          }

  fclose(f);
}



void Print_My_Atoms( reax_system *system )
{
  int   i;
  char  fname[100];
  FILE *fh;

  sprintf( fname, "my_atoms.%d", system->my_rank );
  if( (fh = fopen( fname, "w" )) == NULL )
    {
      fprintf( stderr, "error in opening my_atoms file" );
      MPI_Abort( MPI_COMM_WORLD, FILE_NOT_FOUND );
    }

  for( i = 0; i < system->n; ++i )
    fprintf( fh, "p%-2d %-5d %2d %24.15e%24.15e%24.15e\n",
             system->my_rank,
             system->my_atoms[i].orig_id, system->my_atoms[i].type,
             system->my_atoms[i].x[0],
             system->my_atoms[i].x[1],
             system->my_atoms[i].x[2] );

  fclose( fh );
}


void Print_My_Ext_Atoms( reax_system *system )
{
  int   i;
  char  fname[100];
  FILE *fh;

  sprintf( fname, "my_ext_atoms.%d", system->my_rank );
  if( (fh = fopen( fname, "w" )) == NULL )
    {
      fprintf( stderr, "error in opening my_ext_atoms file" );
      MPI_Abort( MPI_COMM_WORLD, FILE_NOT_FOUND );
    }

  for( i = 0; i < system->N; ++i )
    fprintf( fh, "p%-2d %-5d imprt%-5d %2d %24.15e%24.15e%24.15e\n",
             system->my_rank, system->my_atoms[i].orig_id,
             system->my_atoms[i].imprt_id, system->my_atoms[i].type,
             system->my_atoms[i].x[0],
             system->my_atoms[i].x[1],
             system->my_atoms[i].x[2] );

  fclose( fh );
}


void Print_Far_Neighbors( reax_system *system, reax_list **lists,
                          control_params *control )
{
  char  fname[100];
  int   i, j, nbr, natoms;
  rc_tagint id_i, id_j;
  FILE *fout;
  reax_list *far_nbrs;

  sprintf( fname, "%s.far_nbrs.%d", control->sim_name, system->my_rank );
  fout      = fopen( fname, "w" );
  far_nbrs = (*lists) + FAR_NBRS;
  natoms = system->N;

  for( i = 0; i < natoms; ++i ) {
    id_i = system->my_atoms[i].orig_id;

    for( j = Start_Index(i,far_nbrs); j < End_Index(i,far_nbrs); ++j ) {
      nbr = far_nbrs->select.far_nbr_list[j].nbr;
      id_j = system->my_atoms[nbr].orig_id;

      fprintf( fout, "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
               id_i, id_j, far_nbrs->select.far_nbr_list[j].d,
               far_nbrs->select.far_nbr_list[j].dvec[0],
               far_nbrs->select.far_nbr_list[j].dvec[1],
               far_nbrs->select.far_nbr_list[j].dvec[2] );

      fprintf( fout, "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
               id_j, id_i, far_nbrs->select.far_nbr_list[j].d,
               -far_nbrs->select.far_nbr_list[j].dvec[0],
               -far_nbrs->select.far_nbr_list[j].dvec[1],
               -far_nbrs->select.far_nbr_list[j].dvec[2] );
    }
  }

  fclose( fout );
}


void Print_Sparse_Matrix( reax_system *system, sparse_matrix *A )
{
  int i, j;

  for( i = 0; i < A->n; ++i )
    for( j = A->start[i]; j < A->end[i]; ++j )
      fprintf( stderr, "%d %d %.15e\n",
               system->my_atoms[i].orig_id,
               system->my_atoms[A->entries[j].j].orig_id,
               A->entries[j].val );
}


void Print_Sparse_Matrix2( reax_system *system, sparse_matrix *A, char *fname )
{
  int i, j;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < A->n; ++i )
    for( j = A->start[i]; j < A->end[i]; ++j )
      fprintf( f, "%d %d %.15e\n",
               system->my_atoms[i].orig_id,
               system->my_atoms[A->entries[j].j].orig_id,
               A->entries[j].val );

  fclose(f);
}


void Print_Symmetric_Sparse(reax_system *system, sparse_matrix *A, char *fname)
{
  int i, j;
  reax_atom *ai, *aj;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < A->n; ++i ) {
    ai = &(system->my_atoms[i]);
    for( j = A->start[i]; j < A->end[i]; ++j ) {
      aj = &(system->my_atoms[A->entries[j].j]);
      fprintf( f, "%d %d %.15e\n",
               ai->renumber, aj->renumber, A->entries[j].val );
      if( A->entries[j].j < system->n && ai->renumber != aj->renumber )
        fprintf( f, "%d %d %.15e\n",
                 aj->renumber, ai->renumber, A->entries[j].val );
    }
  }

  fclose(f);
}


void Print_Linear_System( reax_system *system, control_params *control,
                          storage *workspace, int step )
{
  int   i;
  char  fname[100];
  reax_atom *ai;
  FILE *out;

  // print rhs and init guesses for QEq
  sprintf( fname, "%s.p%dstate%d", control->sim_name, system->my_rank, step );
  out = fopen( fname, "w" );
  for( i = 0; i < system->n; i++ ) {
    ai = &(system->my_atoms[i]);
    fprintf( out, "%6d%2d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
             ai->renumber, ai->type, ai->x[0], ai->x[1], ai->x[2],
             workspace->s[i], workspace->b_s[i],
             workspace->t[i], workspace->b_t[i] );
  }
  fclose( out );

  // print QEq coef matrix
  sprintf( fname, "%s.p%dH%d", control->sim_name, system->my_rank, step );
  Print_Symmetric_Sparse( system, workspace->H, fname );

}


void Print_LinSys_Soln( reax_system *system, real *x, real *b_prm, real *b )
{
  int    i;
  char   fname[100];
  FILE  *fout;

  sprintf( fname, "qeq.%d.out", system->my_rank );
  fout = fopen( fname, "w" );

  for( i = 0; i < system->n; ++i )
    fprintf( fout, "%6d%10.4f%10.4f%10.4f\n",
             system->my_atoms[i].orig_id, x[i], b_prm[i], b[i] );

  fclose( fout );
}


void Print_Charges( reax_system *system )
{
  int    i;
  char   fname[100];
  FILE  *fout;

  sprintf( fname, "q.%d.out", system->my_rank );
  fout = fopen( fname, "w" );

  for( i = 0; i < system->n; ++i )
    fprintf( fout, "%6d %10.7f %10.7f %10.7f\n",
             system->my_atoms[i].orig_id,
             system->my_atoms[i].s[0],
             system->my_atoms[i].t[0],
             system->my_atoms[i].q );

  fclose( fout );
}


void Print_Bonds( reax_system *system, reax_list *bonds, char *fname )
{
  int i, j, pj;
  bond_data *pbond;
  bond_order_data *bo_ij;
  FILE *f = fopen( fname, "w" );

  for( i = 0; i < system->N; ++i )
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      pbond = &(bonds->select.bond_list[pj]);
      bo_ij = &(pbond->bo_data);
      j = pbond->nbr;
      fprintf( f, "%8d%8d %24.15f %24.15f\n",
               i, j,//system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
               pbond->d, bo_ij->BO );
    }

  fclose(f);
}


int fn_qsort_intcmp( const void *a, const void *b )
{
  return( *(int *)a - *(int *)b );
}

void Print_Bond_List2( reax_system *system, reax_list *bonds, char *fname )
{
  int i,j, nbr, pj;
  rc_tagint id_i, id_j;
  FILE *f = fopen( fname, "w" );
  int temp[500];
  int num=0;

  for( i = 0; i < system->n; ++i ) {
    num=0;
    id_i = system->my_atoms[i].orig_id;
    fprintf( f, "%6d:", id_i);
    for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj ) {
      nbr = bonds->select.bond_list[pj].nbr;
      id_j = system->my_atoms[nbr].orig_id;
      if( id_i < id_j )
        temp[num++] = id_j;
    }

    qsort(&temp, num, sizeof(int), fn_qsort_intcmp);
    for(j=0; j < num; j++)
      fprintf(f, "%6d", temp[j] );
    fprintf(f, "\n");
  }
}


void Print_Total_Force( reax_system *system, simulation_data *data,
                        storage *workspace )
{
  int    i;

  fprintf( stderr, "step: %d\n", data->step );
  fprintf( stderr, "%6s\t%-38s\n", "atom", "atom.f[0,1,2]");

  for( i = 0; i < system->N; ++i )
    fprintf( stderr, "%6d %f %f %f\n",
             //"%6d%24.15e%24.15e%24.15e\n",
             system->my_atoms[i].orig_id,
             workspace->f[i][0], workspace->f[i][1], workspace->f[i][2] );
}

void Output_Results( reax_system *system, control_params *control,
                     simulation_data *data, reax_list **lists,
                     output_controls *out_control, mpi_datatypes *mpi_data )
{

  if((out_control->energy_update_freq > 0 &&
      data->step%out_control->energy_update_freq == 0) ||
     (out_control->write_steps > 0 &&
      data->step%out_control->write_steps == 0)){
    /* update system-wide energies */
    Compute_System_Energy( system, data, mpi_data->world );

    /* output energies */
    if( system->my_rank == MASTER_NODE &&
        out_control->energy_update_freq > 0 &&
        data->step % out_control->energy_update_freq == 0 ) {

      if( control->virial ){
        fprintf( out_control->prs,
                 "%8d%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n",
                 data->step,
                 data->int_press[0], data->int_press[1], data->int_press[2],
                 data->ext_press[0], data->ext_press[1], data->ext_press[2],
                 data->kin_press );

        fprintf( out_control->prs,
                 "%8s%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f%13.6f\n",
                 "",system->big_box.box_norms[0], system->big_box.box_norms[1],
                 system->big_box.box_norms[2],
                 data->tot_press[0], data->tot_press[1], data->tot_press[2],
                 system->big_box.V );

        fflush( out_control->prs);
      }
    }

    /* write current frame */
    if( out_control->write_steps > 0 &&
        (data->step-data->prev_steps) % out_control->write_steps == 0 ) {
      Append_Frame( system, control, data, lists, out_control, mpi_data );
    }
  }

}
