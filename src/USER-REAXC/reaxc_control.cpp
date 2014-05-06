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
#include "reaxc_control.h"
#include "reaxc_tool_box.h"

char Read_Control_File( char *control_file, control_params* control,
                        output_controls *out_control )
{
  FILE *fp;
  char *s, **tmp;
  int   i,ival;
  real  val;

  /* open control file */
  if ( (fp = fopen( control_file, "r" ) ) == NULL ) {
    fprintf( stderr, "error opening the control file! terminating...\n" );
    MPI_Abort( MPI_COMM_WORLD,  FILE_NOT_FOUND );
  }

  /* assign default values */
  strcpy( control->sim_name, "simulate" );
  control->ensemble        = NVE;
  control->nsteps          = 0;
  control->dt              = 0.25;
  control->nprocs          = 1;
  control->procs_by_dim[0] = 1;
  control->procs_by_dim[1] = 1;
  control->procs_by_dim[2] = 1;
  control->geo_format = 1;

  control->restart          = 0;
  out_control->restart_format = WRITE_BINARY;
  out_control->restart_freq = 0;
  control->reposition_atoms = 0;
  control->restrict_bonds   = 0;
  control->remove_CoM_vel   = 25;
  out_control->debug_level  = 0;
  out_control->energy_update_freq = 0;

  control->reneighbor = 1;
  control->vlist_cut = control->nonb_cut;
  control->bond_cut = 5.0;
  control->bg_cut = 0.3;
  control->thb_cut = 0.001;
  control->thb_cutsq = 0.00001;
  control->hbond_cut = 7.5;

  control->tabulate = 0;

  control->qeq_freq = 1;
  control->q_err = 1e-6;
  control->refactor = 100;
  control->droptol = 1e-2;;

  control->T_init = 0.;
  control->T_final = 300.;
  control->Tau_T = 500.0;
  control->T_mode = 0;
  control->T_rate = 1.;
  control->T_freq = 1.;

  control->P[0] = control->P[1] = control->P[2] = 0.000101325;
  control->Tau_P[0] = control->Tau_P[1] = control->Tau_P[2] = 500.0;
  control->Tau_PT[0] = control->Tau_PT[1] = control->Tau_PT[2] = 500.0;
  control->compressibility = 1.0;
  control->press_mode = 0;
  control->virial = 0;

  out_control->write_steps = 0;
  out_control->traj_compress = 0;
  out_control->traj_method = REG_TRAJ;
  strcpy( out_control->traj_title, "default_title" );
  out_control->atom_info = 0;
  out_control->bond_info = 0;
  out_control->angle_info = 0;

  control->molecular_analysis = 0;
  control->dipole_anal = 0;
  control->freq_dipole_anal = 0;
  control->diffusion_coef = 0;
  control->freq_diffusion_coef = 0;
  control->restrict_type = 0;

  /* memory allocations */
  s = (char*) malloc(sizeof(char)*MAX_LINE);
  tmp = (char**) malloc(sizeof(char*)*MAX_TOKENS);
  for (i=0; i < MAX_TOKENS; i++)
    tmp[i] = (char*) malloc(sizeof(char)*MAX_LINE);

  /* read control parameters file */
  while (!feof(fp)) {
    fgets( s, MAX_LINE, fp );
    Tokenize( s, &tmp );

    if( strcmp(tmp[0], "simulation_name") == 0 ) {
      strcpy( control->sim_name, tmp[1] );
    }
    else if( strcmp(tmp[0], "ensemble_type") == 0 ) {
      ival = atoi(tmp[1]);
      control->ensemble = ival;
      if( ival == iNPT || ival ==sNPT || ival == NPT )
        control->virial = 1;
    }
    else if( strcmp(tmp[0], "nsteps") == 0 ) {
      ival = atoi(tmp[1]);
      control->nsteps = ival;
    }
    else if( strcmp(tmp[0], "dt") == 0) {
      val = atof(tmp[1]);
      control->dt = val * 1.e-3;  // convert dt from fs to ps!
    }
    else if( strcmp(tmp[0], "proc_by_dim") == 0 ) {
      ival = atoi(tmp[1]);
      control->procs_by_dim[0] = ival;
      ival = atoi(tmp[2]);
      control->procs_by_dim[1] = ival;
      ival = atoi(tmp[3]);
      control->procs_by_dim[2] = ival;

      control->nprocs = control->procs_by_dim[0]*control->procs_by_dim[1]*
        control->procs_by_dim[2];
    }
    else if( strcmp(tmp[0], "random_vel") == 0 ) {
      ival = atoi(tmp[1]);
      control->random_vel = ival;
    }
    else if( strcmp(tmp[0], "restart_format") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->restart_format = ival;
    }
    else if( strcmp(tmp[0], "restart_freq") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->restart_freq = ival;
    }
    else if( strcmp(tmp[0], "reposition_atoms") == 0 ) {
      ival = atoi(tmp[1]);
      control->reposition_atoms = ival;
    }
    else if( strcmp(tmp[0], "restrict_bonds") == 0 ) {
      ival = atoi( tmp[1] );
      control->restrict_bonds = ival;
    }
    else if( strcmp(tmp[0], "remove_CoM_vel") == 0 ) {
      ival = atoi(tmp[1]);
      control->remove_CoM_vel = ival;
    }
    else if( strcmp(tmp[0], "debug_level") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->debug_level = ival;
    }
    else if( strcmp(tmp[0], "energy_update_freq") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->energy_update_freq = ival;
    }
    else if( strcmp(tmp[0], "reneighbor") == 0 ) {
      ival = atoi( tmp[1] );
      control->reneighbor = ival;
    }
    else if( strcmp(tmp[0], "vlist_buffer") == 0 ) {
      val = atof(tmp[1]);
      control->vlist_cut= val + control->nonb_cut;
    }
    else if( strcmp(tmp[0], "nbrhood_cutoff") == 0 ) {
      val = atof(tmp[1]);
      control->bond_cut = val;
    }
    else if( strcmp(tmp[0], "bond_graph_cutoff") == 0 ) {
      val = atof(tmp[1]);
      control->bg_cut = val;
    }
    else if( strcmp(tmp[0], "thb_cutoff") == 0 ) {
      val = atof(tmp[1]);
      control->thb_cut = val;
    }
    else if( strcmp(tmp[0], "thb_cutoff_sq") == 0 ) {
      val = atof(tmp[1]);
      control->thb_cutsq = val;
    }
    else if( strcmp(tmp[0], "hbond_cutoff") == 0 ) {
      val = atof( tmp[1] );
      control->hbond_cut = val;
    }
    else if( strcmp(tmp[0], "ghost_cutoff") == 0 ) {
      val = atof(tmp[1]);
      control->user_ghost_cut = val;
    }
    else if( strcmp(tmp[0], "tabulate_long_range") == 0 ) {
      ival = atoi( tmp[1] );
      control->tabulate = ival;
    }
    else if( strcmp(tmp[0], "qeq_freq") == 0 ) {
      ival = atoi( tmp[1] );
      control->qeq_freq = ival;
    }
    else if( strcmp(tmp[0], "q_err") == 0 ) {
      val = atof( tmp[1] );
      control->q_err = val;
    }
    else if( strcmp(tmp[0], "ilu_refactor") == 0 ) {
      ival = atoi( tmp[1] );
      control->refactor = ival;
    }
    else if( strcmp(tmp[0], "ilu_droptol") == 0 ) {
      val = atof( tmp[1] );
      control->droptol = val;
    }
    else if( strcmp(tmp[0], "temp_init") == 0 ) {
      val = atof(tmp[1]);
      control->T_init = val;

      if( control->T_init < 0.1 )
        control->T_init = 0.1;
    }
    else if( strcmp(tmp[0], "temp_final") == 0 ) {
      val = atof(tmp[1]);
      control->T_final = val;

      if( control->T_final < 0.1 )
        control->T_final = 0.1;
    }
    else if( strcmp(tmp[0], "t_mass") == 0 ) {
      val = atof(tmp[1]);
      control->Tau_T = val * 1.e-3;    // convert t_mass from fs to ps
    }
    else if( strcmp(tmp[0], "t_mode") == 0 ) {
      ival = atoi(tmp[1]);
      control->T_mode = ival;
    }
    else if( strcmp(tmp[0], "t_rate") == 0 ) {
      val = atof(tmp[1]);
      control->T_rate = val;
    }
    else if( strcmp(tmp[0], "t_freq") == 0 ) {
      val = atof(tmp[1]);
      control->T_freq = val;
    }
    else if( strcmp(tmp[0], "pressure") == 0 ) {
      if( control->ensemble == iNPT ) {
        control->P[0] = control->P[1] = control->P[2] = atof(tmp[1]);
      }
      else if( control->ensemble == sNPT ) {
        control->P[0] = atof(tmp[1]);
        control->P[1] = atof(tmp[2]);
        control->P[2] = atof(tmp[3]);
      }
    }
    else if( strcmp(tmp[0], "p_mass") == 0 ) {
      // convert p_mass from fs to ps
      if( control->ensemble == iNPT ) {
        control->Tau_P[0] = control->Tau_P[1] = control->Tau_P[2] =
          atof(tmp[1]) * 1.e-3;
      }
      else if( control->ensemble == sNPT ) {
        control->Tau_P[0] = atof(tmp[1]) * 1.e-3;
        control->Tau_P[1] = atof(tmp[2]) * 1.e-3;
        control->Tau_P[2] = atof(tmp[3]) * 1.e-3;
      }
    }
    else if( strcmp(tmp[0], "pt_mass") == 0 ) {
      val = atof(tmp[1]);
      control->Tau_PT[0] = control->Tau_PT[1] = control->Tau_PT[2] =
        val * 1.e-3;  // convert pt_mass from fs to ps
    }
    else if( strcmp(tmp[0], "compress") == 0 ) {
      val = atof(tmp[1]);
      control->compressibility = val;
    }
    else if( strcmp(tmp[0], "press_mode") == 0 ) {
      ival = atoi(tmp[1]);
      control->press_mode = ival;
    }
    else if( strcmp(tmp[0], "geo_format") == 0 ) {
      ival = atoi( tmp[1] );
      control->geo_format = ival;
    }
    else if( strcmp(tmp[0], "write_freq") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->write_steps = ival;
    }
    else if( strcmp(tmp[0], "traj_compress") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->traj_compress = ival;
    }
    else if( strcmp(tmp[0], "traj_method") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->traj_method = ival;
    }
    else if( strcmp(tmp[0], "traj_title") == 0 ) {
      strcpy( out_control->traj_title, tmp[1] );
    }
    else if( strcmp(tmp[0], "atom_info") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 4;
    }
    else if( strcmp(tmp[0], "atom_velocities") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 2;
    }
    else if( strcmp(tmp[0], "atom_forces") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 1;
    }
    else if( strcmp(tmp[0], "bond_info") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->bond_info = ival;
    }
    else if( strcmp(tmp[0], "angle_info") == 0 ) {
      ival = atoi(tmp[1]);
      out_control->angle_info = ival;
    }
    else if( strcmp(tmp[0], "molecular_analysis") == 0 ) {
      ival = atoi(tmp[1]);
      control->molecular_analysis = ival;
    }
    else if( strcmp(tmp[0], "ignore") == 0 ) {
      control->num_ignored = atoi(tmp[1]);
      for( i = 0; i < control->num_ignored; ++i )
        control->ignore[atoi(tmp[i+2])] = 1;
    }
    else if( strcmp(tmp[0], "dipole_anal") == 0 ) {
      ival = atoi(tmp[1]);
      control->dipole_anal = ival;
    }
    else if( strcmp(tmp[0], "freq_dipole_anal") == 0 ) {
      ival = atoi(tmp[1]);
      control->freq_dipole_anal = ival;
    }
    else if( strcmp(tmp[0], "diffusion_coef") == 0 ) {
      ival = atoi(tmp[1]);
      control->diffusion_coef = ival;
    }
    else if( strcmp(tmp[0], "freq_diffusion_coef") == 0 ) {
      ival = atoi(tmp[1]);
      control->freq_diffusion_coef = ival;
    }
    else if( strcmp(tmp[0], "restrict_type") == 0 ) {
      ival = atoi(tmp[1]);
      control->restrict_type = ival;
    }
    else {
      fprintf( stderr, "WARNING: unknown parameter %s\n", tmp[0] );
      MPI_Abort( MPI_COMM_WORLD, 15 );
    }
  }

  /* determine target T */
  if( control->T_mode == 0 )
    control->T = control->T_final;
  else control->T = control->T_init;

  /* free memory allocations at the top */
  for( i = 0; i < MAX_TOKENS; i++ )
    free( tmp[i] );
  free( tmp );
  free( s );

  fclose(fp);

  return SUCCESS;
}
