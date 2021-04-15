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

#include "reaxc_io_tools.h"
#include <cstdio>
#include <cstring>
#include "reaxc_defs.h"
#include "reaxc_system_props.h"
#include "reaxc_traj.h"

int Init_Output_Files(reax_system *system, control_params *control,
                      output_controls *out_control, MPI_Comm world, char *msg)
{
  int ret;

  if (out_control->write_steps > 0) {
    ret = Init_Traj( system, control, out_control, world, msg );
    if (ret == FAILURE)
      return ret;
  }
  return SUCCESS;
}

/************************ close output files ************************/
int Close_Output_Files(reax_system *system, output_controls *out_control)
{
  if (out_control->write_steps > 0)
    End_Traj( system->my_rank, out_control );

  return SUCCESS;
}


void Output_Results(reax_system *system, control_params *control,
                    simulation_data *data, reax_list **lists,
                    output_controls *out_control, MPI_Comm world)
{

  if ((out_control->energy_update_freq > 0 &&
      data->step%out_control->energy_update_freq == 0) ||
     (out_control->write_steps > 0 &&
      data->step%out_control->write_steps == 0)) {
    /* update system-wide energies */
    Compute_System_Energy(system, data, world);

    /* write current frame */
    if ( out_control->write_steps > 0 && data->step % out_control->write_steps == 0) {
      Append_Frame( system, control, data, lists, out_control, world);
    }
  }

}
