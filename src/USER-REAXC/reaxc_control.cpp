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

#include "reaxc_control.h"
#include <cstdlib>
#include <cstring>
#include "reaxc_defs.h"
#include "reaxc_tool_box.h"

#include "error.h"
#include "utils.h"

using LAMMPS_NS::utils::getsyserror;

void Read_Control_File(const char *control_file, control_params *control,
                        output_controls *out_control)
{
  FILE *fp;
  char *s, **tmp;
  int   i,ival;
  double  val;

  /* open control file */
  fp = fopen(control_file, "r");
  if (!fp)
    control->error_ptr->one(FLERR,fmt::format("The control file {} cannot be "
                                              "opened: {}",control_file,
                                              getsyserror()));
  /* assign default values */
  strcpy( control->sim_name, "simulate" );
  control->nthreads        = 1;

  out_control->energy_update_freq = 0;

  control->bond_cut = 5.0;
  control->bg_cut = 0.3;
  control->thb_cut = 0.001;
  control->thb_cutsq = 0.00001;
  control->hbond_cut = 7.5;

  control->tabulate = 0;

  control->virial = 0;

  out_control->write_steps = 0;
  strcpy( out_control->traj_title, "default_title" );
  out_control->atom_info = 0;
  out_control->bond_info = 0;
  out_control->angle_info = 0;

  /* memory allocations */
  s = (char*) malloc(sizeof(char)*MAX_LINE);
  tmp = (char**) malloc(sizeof(char*)*MAX_TOKENS);
  for (i=0; i < MAX_TOKENS; i++)
    tmp[i] = (char*) malloc(sizeof(char)*MAX_LINE);

  /* read control parameters file */
  while (!feof(fp)) {
    fgets( s, MAX_LINE, fp );
    Tokenize( s, &tmp );

    if (strcmp(tmp[0], "simulation_name") == 0) {
      strcpy( control->sim_name, tmp[1] );
    }
    else if (strcmp(tmp[0], "ensemble_type") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "nsteps") == 0) {
      ; // ignore
    }
    else if ( strcmp(tmp[0], "dt") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "proc_by_dim") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "random_vel") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "restart_format") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "restart_freq") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "reposition_atoms") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "restrict_bonds") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "remove_CoM_vel") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "debug_level") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "energy_update_freq") == 0) {
      ival = atoi(tmp[1]);
      out_control->energy_update_freq = ival;
    }
    else if (strcmp(tmp[0], "reneighbor") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "vlist_buffer") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "nbrhood_cutoff") == 0) {
      val = atof(tmp[1]);
      control->bond_cut = val;
    }
    else if (strcmp(tmp[0], "bond_graph_cutoff") == 0) {
      val = atof(tmp[1]);
      control->bg_cut = val;
    }
    else if (strcmp(tmp[0], "thb_cutoff") == 0) {
      val = atof(tmp[1]);
      control->thb_cut = val;
    }
    else if (strcmp(tmp[0], "thb_cutoff_sq") == 0) {
      val = atof(tmp[1]);
      control->thb_cutsq = val;
    }
    else if (strcmp(tmp[0], "hbond_cutoff") == 0) {
      val = atof( tmp[1] );
      control->hbond_cut = val;
    }
    else if (strcmp(tmp[0], "ghost_cutoff") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "tabulate_long_range") == 0) {
      ival = atoi( tmp[1] );
      control->tabulate = ival;
    }
    else if (strcmp(tmp[0], "qeq_freq") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "q_err") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "ilu_refactor") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "ilu_droptol") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "temp_init") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "temp_final") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "t_mass") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "t_mode") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "t_rate") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "t_freq") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "pressure") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "p_mass") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "pt_mass") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "compress") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "press_mode") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "geo_format") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "write_freq") == 0) {
      ival = atoi(tmp[1]);
      out_control->write_steps = ival;
    }
    else if (strcmp(tmp[0], "traj_compress") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "traj_method") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "traj_title") == 0) {
      strcpy( out_control->traj_title, tmp[1] );
    }
    else if (strcmp(tmp[0], "atom_info") == 0) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 4;
    }
    else if (strcmp(tmp[0], "atom_velocities") == 0) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 2;
    }
    else if (strcmp(tmp[0], "atom_forces") == 0) {
      ival = atoi(tmp[1]);
      out_control->atom_info += ival * 1;
    }
    else if (strcmp(tmp[0], "bond_info") == 0) {
      ival = atoi(tmp[1]);
      out_control->bond_info = ival;
    }
    else if (strcmp(tmp[0], "angle_info") == 0) {
      ival = atoi(tmp[1]);
      out_control->angle_info = ival;
    }
    else if (strcmp(tmp[0], "molecular_analysis") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "ignore") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "dipole_anal") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "freq_dipole_anal") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "diffusion_coef") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "freq_diffusion_coef") == 0) {
      ; // ignore
    }
    else if (strcmp(tmp[0], "restrict_type") == 0) {
      ; // ignore
    }
    else {
      char errmsg[128];
      snprintf(errmsg,128,"Unknown parameter %s in the control file", tmp[0]);
      control->error_ptr->all(FLERR, errmsg);
    }
  }

  /* free memory allocations at the top */
  for (i = 0; i < MAX_TOKENS; i++)
    free( tmp[i] );
  free( tmp );
  free( s );

  fclose(fp);
}
