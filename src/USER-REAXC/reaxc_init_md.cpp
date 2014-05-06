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
#include "reaxc_init_md.h"
#include "reaxc_allocate.h"
#include "reaxc_forces.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_reset_tools.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

int Init_System( reax_system *system, control_params *control, char *msg )
{
  int i;
  reax_atom *atom;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  // determine the local and total capacity

  system->local_cap = MAX( (int)(system->n * safezone), mincap);
  system->total_cap = MAX( (int)(system->N * safezone), mincap);

  /* estimate numH and Hcap */
  system->numH = 0;
  if( control->hbond_cut > 0 )
    for( i = 0; i < system->n; ++i ) {
      atom = &(system->my_atoms[i]);
      if (atom->type < 0) continue;
      if( system->reax_param.sbp[ atom->type ].p_hbond == 1 )
        atom->Hindex = system->numH++;
      else atom->Hindex = -1;
    }
  system->Hcap = (int)(MAX( system->numH * saferzone, mincap ));

  return SUCCESS;
}


int Init_Simulation_Data( reax_system *system, control_params *control,
                          simulation_data *data, char *msg )
{
  Reset_Simulation_Data( data, control->virial );

  /* initialize the timer(s) */
  if( system->my_rank == MASTER_NODE ) {
    data->timing.start = Get_Time( );
  }

  data->step = data->prev_steps = 0;

  return SUCCESS;
}

void Init_Taper( control_params *control,  storage *workspace, MPI_Comm comm )
{
  real d1, d7;
  real swa, swa2, swa3;
  real swb, swb2, swb3;

  swa = control->nonb_low;
  swb = control->nonb_cut;

  if( fabs( swa ) > 0.01 )
    fprintf( stderr, "Warning: non-zero lower Taper-radius cutoff\n" );

  if( swb < 0 ) {
    fprintf( stderr, "Negative upper Taper-radius cutoff\n" );
    MPI_Abort( comm,  INVALID_INPUT );
  }
  else if( swb < 5 )
    fprintf( stderr, "Warning: very low Taper-radius cutoff: %f\n", swb );

  d1 = swb - swa;
  d7 = pow( d1, 7.0 );
  swa2 = SQR( swa );
  swa3 = CUBE( swa );
  swb2 = SQR( swb );
  swb3 = CUBE( swb );

  workspace->Tap[7] =  20.0 / d7;
  workspace->Tap[6] = -70.0 * (swa + swb) / d7;
  workspace->Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  workspace->Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / d7;
  workspace->Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / d7;
  workspace->Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  workspace->Tap[1] = 140.0 * swa3 * swb3 / d7;
  workspace->Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 +
                     7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;
}


int Init_Workspace( reax_system *system, control_params *control,
                    storage *workspace, MPI_Comm comm, char *msg )
{
  int ret;

  ret = Allocate_Workspace( system, control, workspace,
                            system->local_cap, system->total_cap, comm, msg );
  if( ret != SUCCESS )
    return ret;

  memset( &(workspace->realloc), 0, sizeof(reallocate_data) );
  Reset_Workspace( system, workspace );

  /* Initialize the Taper function */
  Init_Taper( control, workspace, comm );

  return SUCCESS;
}


/************** setup communication data structures  **************/
int Init_MPI_Datatypes( reax_system *system, storage *workspace,
                        mpi_datatypes *mpi_data, MPI_Comm comm, char *msg )
{

  /* setup the world */
  mpi_data->world = comm;
  MPI_Comm_size( comm, &(system->wsize) );

  return SUCCESS;
}

int  Init_Lists( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace, reax_list **lists,
                 mpi_datatypes *mpi_data, char *msg )
{
  int i, total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
  int *hb_top, *bond_top;
  MPI_Comm comm;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  comm = mpi_data->world;
  bond_top = (int*) calloc( system->total_cap, sizeof(int) );
  hb_top = (int*) calloc( system->local_cap, sizeof(int) );
  Estimate_Storages( system, control, lists,
                     &Htop, hb_top, bond_top, &num_3body, comm );

  if( control->hbond_cut > 0 ) {
    /* init H indexes */
    total_hbonds = 0;
    for( i = 0; i < system->n; ++i ) {
      system->my_atoms[i].num_hbonds = hb_top[i];
      total_hbonds += hb_top[i];
    }
    total_hbonds = (int)(MAX( total_hbonds*saferzone, mincap*MIN_HBONDS ));

    if( !Make_List( system->Hcap, total_hbonds, TYP_HBOND,
                    *lists+HBONDS, comm ) ) {
      fprintf( stderr, "not enough space for hbonds list. terminating!\n" );
      MPI_Abort( comm, INSUFFICIENT_MEMORY );
    }
  }

  total_bonds = 0;
  for( i = 0; i < system->N; ++i ) {
    system->my_atoms[i].num_bonds = bond_top[i];
    total_bonds += bond_top[i];
  }
  bond_cap = (int)(MAX( total_bonds*safezone, mincap*MIN_BONDS ));

  if( !Make_List( system->total_cap, bond_cap, TYP_BOND,
                  *lists+BONDS, comm ) ) {
    fprintf( stderr, "not enough space for bonds list. terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  /* 3bodies list */
  cap_3body = (int)(MAX( num_3body*safezone, MIN_3BODIES ));
  if( !Make_List( bond_cap, cap_3body, TYP_THREE_BODY,
                  *lists+THREE_BODIES, comm ) ){
    fprintf( stderr, "Problem in initializing angles list. Terminating!\n" );
    MPI_Abort( comm, INSUFFICIENT_MEMORY );
  }

  free( hb_top );
  free( bond_top );

  return SUCCESS;
}

void Initialize( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace,
                 reax_list **lists, output_controls *out_control,
                 mpi_datatypes *mpi_data, MPI_Comm comm )
{
  char msg[MAX_STR];


  if( Init_MPI_Datatypes(system, workspace, mpi_data, comm, msg) == FAILURE ) {
    fprintf( stderr, "p%d: init_mpi_datatypes: could not create datatypes\n",
             system->my_rank );
    fprintf( stderr, "p%d: mpi_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }

  if( Init_System(system, control, msg) == FAILURE ){
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }

  if( Init_Simulation_Data( system, control, data, msg ) == FAILURE ) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: sim_data couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }

  if( Init_Workspace( system, control, workspace, mpi_data->world, msg ) ==
      FAILURE ) {
    fprintf( stderr, "p%d:init_workspace: not enough memory\n",
             system->my_rank );
    fprintf( stderr, "p%d:workspace couldn't be initialized! terminating.\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }

  if( Init_Lists( system, control, data, workspace, lists, mpi_data, msg ) ==
      FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: system could not be initialized! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }

  if( Init_Output_Files(system,control,out_control,mpi_data,msg)== FAILURE) {
    fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
    fprintf( stderr, "p%d: could not open output files! terminating...\n",
             system->my_rank );
    MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
  }

  if( control->tabulate ) {
    if( Init_Lookup_Tables( system, control, workspace, mpi_data, msg ) == FAILURE ) {
      fprintf( stderr, "p%d: %s\n", system->my_rank, msg );
      fprintf( stderr, "p%d: couldn't create lookup table! terminating.\n",
               system->my_rank );
      MPI_Abort( mpi_data->world, CANNOT_INITIALIZE );
    }
  }


  Init_Force_Functions( control );
}
