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

#include "reaxc_io_tools.h"
#include <cstdio>
#include <cstring>
#include "reaxc_defs.h"
#include "reaxc_system_props.h"
#include "reaxc_traj.h"

int Init_Output_Files( reax_system *system, control_params *control,
                       output_controls *out_control, mpi_datatypes *mpi_data,
                       char *msg )
{
  char temp[MAX_STR+8];
  int ret;

  if (out_control->write_steps > 0) {
    ret = Init_Traj( system, control, out_control, mpi_data, msg );
    if (ret == FAILURE)
      return ret;
  }

  if (system->my_rank == MASTER_NODE) {
    /* These files are written only by the master node */
    if (out_control->energy_update_freq > 0) {

      /* init potentials file */
      sprintf( temp, "%s.pot", control->sim_name );
      if ((out_control->pot = fopen( temp, "w" )) != NULL) {
        fflush( out_control->pot );
      } else {
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
      if ((out_control->prs = fopen( temp, "w" )) != NULL) {
        fprintf(out_control->prs,"%8s%13s%13s%13s%13s%13s%13s%13s\n",
                "step", "Pint/norm[x]", "Pint/norm[y]", "Pint/norm[z]",
                "Pext/Ptot[x]", "Pext/Ptot[y]", "Pext/Ptot[z]", "Pkin/V" );
        fflush( out_control->prs );
      } else {
        strcpy(msg,"init_out_controls: .prs file couldn't be opened\n");
        return FAILURE;
      }
    }
  }

  return SUCCESS;
}


/************************ close output files ************************/
int Close_Output_Files( reax_system *system, control_params * /* control */,
                        output_controls *out_control, mpi_datatypes * /*mpi_data*/ )
{
  if (out_control->write_steps > 0)
    End_Traj( system->my_rank, out_control );

  if (system->my_rank == MASTER_NODE) {
    if (out_control->pot) {
      fclose( out_control->pot );
      out_control->pot = NULL;
    }

    if (out_control->prs) {
      fclose(out_control->prs);
      out_control->prs = NULL;
    }
  }

  return SUCCESS;
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

      if (control->virial && out_control->prs) {
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
