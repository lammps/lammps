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

#include "error.h"

int Init_System( reax_system *system, control_params *control, char * /*msg*/ )
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
  if (control->hbond_cut > 0)
    for( i = 0; i < system->n; ++i ) {
      atom = &(system->my_atoms[i]);
      if (system->reax_param.sbp[ atom->type ].p_hbond == 1 && atom->type >= 0)
        atom->Hindex = system->numH++;
      else atom->Hindex = -1;
    }
  system->Hcap = (int)(MAX( system->numH * saferzone, mincap ));

  return SUCCESS;
}


int Init_Simulation_Data( reax_system *system, control_params *control,
                          simulation_data *data, char * /*msg*/ )
{
  Reset_Simulation_Data( data, control->virial );

  /* initialize the timer(s) */
  if (system->my_rank == MASTER_NODE) {
    data->timing.start = Get_Time( );
  }

  data->step = data->prev_steps = 0;

  return SUCCESS;
}

void Init_Taper( control_params *control,  storage *workspace )
{
  double d1, d7;
  double swa, swa2, swa3;
  double swb, swb2, swb3;

  swa = control->nonb_low;
  swb = control->nonb_cut;

  if (fabs( swa ) > 0.01 && control->me == 0)
    control->error_ptr->warning( FLERR, "Non-zero lower Taper-radius cutoff" );

  if (swb < 0) {
    control->error_ptr->all(FLERR,"Negative upper Taper-radius cutoff");
  }
  else if( swb < 5 && control->me == 0) {
    char errmsg[256];
    snprintf(errmsg, 256, "Very low Taper-radius cutoff: %f", swb );
    control->error_ptr->warning( FLERR, errmsg );
  }

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
  workspace->Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
                     7.0*swa*swb3*swb3 + swb3*swb3*swb ) / d7;
}


int Init_Workspace( reax_system *system, control_params *control,
                    storage *workspace, char *msg )
{
  int ret;

  ret = Allocate_Workspace( system, control, workspace,
                            system->local_cap, system->total_cap, msg );
  if (ret != SUCCESS)
    return ret;

  memset( &(workspace->realloc), 0, sizeof(reallocate_data) );
  Reset_Workspace( system, workspace );

  /* Initialize the Taper function */
  Init_Taper( control, workspace);

  return SUCCESS;
}


/************** setup communication data structures  **************/
int Init_MPI_Datatypes( reax_system *system, storage * /*workspace*/,
                        mpi_datatypes *mpi_data, MPI_Comm comm, char * /*msg*/ )
{

  /* setup the world */
  mpi_data->world = comm;
  MPI_Comm_size( comm, &(system->wsize) );

  return SUCCESS;
}

int  Init_Lists( reax_system *system, control_params *control,
                 simulation_data * /*data*/, storage * /*workspace*/, reax_list **lists,
                 mpi_datatypes * /*mpi_data*/, char * /*msg*/ )
{
  int i, total_hbonds, total_bonds, bond_cap, num_3body, cap_3body, Htop;
  int *hb_top, *bond_top;

  int mincap = system->mincap;
  double safezone = system->safezone;
  double saferzone = system->saferzone;

  bond_top = (int*) calloc( system->total_cap, sizeof(int) );
  hb_top = (int*) calloc( system->local_cap, sizeof(int) );
  Estimate_Storages( system, control, lists,
                     &Htop, hb_top, bond_top, &num_3body);

  if (control->hbond_cut > 0) {
    /* init H indexes */
    total_hbonds = 0;
    for( i = 0; i < system->n; ++i ) {
      system->my_atoms[i].num_hbonds = hb_top[i];
      total_hbonds += hb_top[i];
    }
    total_hbonds = (int)(MAX( total_hbonds*saferzone, mincap*MIN_HBONDS ));

    if( !Make_List( system->Hcap, total_hbonds, TYP_HBOND,
                    *lists+HBONDS ) ) {
      control->error_ptr->one(FLERR, "Not enough space for hbonds list.");
    }
    (*lists+HBONDS)->error_ptr = system->error_ptr;
  }

  total_bonds = 0;
  for( i = 0; i < system->N; ++i ) {
    system->my_atoms[i].num_bonds = bond_top[i];
    total_bonds += bond_top[i];
  }
  bond_cap = (int)(MAX( total_bonds*safezone, mincap*MIN_BONDS ));

  if( !Make_List( system->total_cap, bond_cap, TYP_BOND,
                  *lists+BONDS ) ) {
    control->error_ptr->one(FLERR, "Not enough space for bonds list.");
  }
  (*lists+BONDS)->error_ptr = system->error_ptr;

  /* 3bodies list */
  cap_3body = (int)(MAX( num_3body*safezone, MIN_3BODIES ));
  if( !Make_List( bond_cap, cap_3body, TYP_THREE_BODY,
                  *lists+THREE_BODIES ) ){
    control->error_ptr->one(FLERR,"Problem in initializing angles list.");
  }
  (*lists+THREE_BODIES)->error_ptr = system->error_ptr;

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

  if (Init_MPI_Datatypes(system, workspace, mpi_data, comm, msg) == FAILURE) {
    control->error_ptr->one(FLERR,"Could not create datatypes");
  }

  if (Init_System(system, control, msg) == FAILURE) {
    control->error_ptr->one(FLERR,"System could not be initialized");
  }

  if (Init_Simulation_Data( system, control, data, msg ) == FAILURE) {
    control->error_ptr->one(FLERR,"Sim_data could not be initialized");
  }

  if (Init_Workspace( system, control, workspace, msg ) ==
      FAILURE) {
    control->error_ptr->one(FLERR,"Workspace could not be initialized");
  }

  if (Init_Lists( system, control, data, workspace, lists, mpi_data, msg ) ==
      FAILURE) {
    control->error_ptr->one(FLERR,"Lists could not be initialized");
    }

  if (Init_Output_Files(system,control,out_control,mpi_data,msg)== FAILURE) {
    control->error_ptr->one(FLERR,"Could not open output files");
  }

  if (control->tabulate) {
    if (Init_Lookup_Tables( system, control, workspace, mpi_data, msg ) == FAILURE) {
    control->error_ptr->one(FLERR,"Lookup table could not be created");
    }
  }


  Init_Force_Functions( control );
}
