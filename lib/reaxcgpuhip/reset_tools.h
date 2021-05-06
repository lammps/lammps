/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

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

#ifndef __RESET_TOOLS_H_
#define __RESET_TOOLS_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void Reset_Pressures_Host( simulation_data * const );

void Reset_Simulation_Data_Host( simulation_data * const );

void Reset_Timing_Host( reax_timing * const );

void Reset_Workspace_Host( reax_system * const, storage * const );

void Reset_Lists_Host( reax_system * const, control_params * const,
        storage * const, reax_list** const );

void Reset_Grid_Host( grid * const );

void Reset_Out_Buffers_Host( mpi_out_data * const, int );

void Reset_Host( reax_system * const, control_params * const,
        simulation_data * const, storage * const, reax_list** const );

#ifdef TEST_FORCES
void Reset_Test_Forces( reax_system * const, storage * const );
#endif

#ifdef __cplusplus
}
#endif


#endif
