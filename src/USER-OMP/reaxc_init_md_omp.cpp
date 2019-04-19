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
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reaxc_omp.h"
#include "reaxc_init_md_omp.h"
#include "reaxc_allocate.h"
#include "reaxc_forces.h"
#include "reaxc_forces_omp.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_reset_tools.h"
#include "reaxc_system_props.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

// Functions defined in reaxc_init_md.cpp
extern int Init_MPI_Datatypes(reax_system*, storage*, mpi_datatypes*, MPI_Comm, char*);
extern int Init_System(reax_system*, control_params*, char*);
extern int Init_Simulation_Data(reax_system*, control_params*, simulation_data*, char*);
extern int Init_Workspace(reax_system*, control_params*, storage*, char*);

/* ---------------------------------------------------------------------- */

int Init_ListsOMP( reax_system *system, control_params *control,
                 simulation_data * /* data */, storage * /* workspace */,
                 reax_list **lists, mpi_datatypes *mpi_data, char * /* msg */)
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
                     &Htop, hb_top, bond_top, &num_3body );

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
      system->error_ptr->one( FLERR, "Not enough space for hbonds list. Terminating!" );
    }
  }

  total_bonds = 0;
  for( i = 0; i < system->N; ++i ) {
    system->my_atoms[i].num_bonds = bond_top[i];
    total_bonds += bond_top[i];
  }
  bond_cap = (int)(MAX( total_bonds*safezone, mincap*MIN_BONDS ));

  if( !Make_List( system->total_cap, bond_cap, TYP_BOND,
                  *lists+BONDS ) ) {
    system->error_ptr->one( FLERR, "Not enough space for bonds list. Terminating!\n" );
  }

  int nthreads = control->nthreads;
  reax_list *bonds = (*lists)+BONDS;

  for (i = 0; i < bonds->num_intrs; ++i)
    bonds->select.bond_list[i].bo_data.CdboReduction =
      (double*) smalloc(system->error_ptr, sizeof(double)*nthreads, "CdboReduction");

  /* 3bodies list */
  cap_3body = (int)(MAX( num_3body*safezone, MIN_3BODIES ));
  if( !Make_List( bond_cap, cap_3body, TYP_THREE_BODY,
                  *lists+THREE_BODIES ) ){

    system->error_ptr->one( FLERR, "Problem in initializing angles list. Terminating!" );
  }

  free( hb_top );
  free( bond_top );

  return SUCCESS;
}

/* ---------------------------------------------------------------------- */

// The only difference with the MPI-only function is calls to Init_ListsOMP and Init_Force_FunctionsOMP().
void InitializeOMP( reax_system *system, control_params *control,
                 simulation_data *data, storage *workspace,
                 reax_list **lists, output_controls *out_control,
                 mpi_datatypes *mpi_data, MPI_Comm comm )
{
  char msg[MAX_STR];
  char errmsg[512];


  if (Init_MPI_Datatypes(system, workspace, mpi_data, comm, msg) == FAILURE) {
    system->error_ptr->one( FLERR, "init_mpi_datatypes: could not create datatypes. "
    "Mpi_data couldn't be initialized! Terminating.");
  }

  if (Init_System(system, control, msg) == FAILURE) {
    snprintf( errmsg, 512, "Error on: %s. "
    "System could not be initialized! Terminating.", msg );
    system->error_ptr->one(FLERR, errmsg);
  }

  if (Init_Simulation_Data( system, control, data, msg ) == FAILURE) {
    snprintf( errmsg, 512, "Error on: %s. "
    "Sim_data couldn't be initialized! Terminating.", msg );
    system->error_ptr->one(FLERR, errmsg);
  }

  if (Init_Workspace( system, control, workspace, msg ) ==
      FAILURE) {
    system->error_ptr->one(FLERR, "init_workspace: not enough memory. "
    "Workspace couldn't be initialized! Terminating.");
  }

  if (Init_ListsOMP( system, control, data, workspace, lists, mpi_data, msg ) ==
      FAILURE) {
      snprintf( errmsg, 512, "Error on: %s. "
      "System could not be initialized! Terminating.", msg );
      system->error_ptr->one(FLERR, errmsg);
    }

  if (Init_Output_Files(system,control,out_control,mpi_data,msg)== FAILURE) {
    snprintf( errmsg, 512, "Error on: %s"
    "Could not open output files! Terminating.", msg );
    system->error_ptr->one(FLERR, errmsg);
  }

  if (control->tabulate) {
    if (Init_Lookup_Tables( system, control, workspace, mpi_data, msg ) == FAILURE) {
      snprintf( errmsg, 512, "Error on: %s."
      " Couldn't create lookup table! Terminating.", msg );
      system->error_ptr->one(FLERR, errmsg);
    }
  }


  Init_Force_FunctionsOMP( control );
}

