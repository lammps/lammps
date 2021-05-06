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

#ifndef __CUDA_INTEGRATE_H_
#define __CUDA_INTEGRATE_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void bNVP_scale_velocities( reax_system *, real, rvec );

int Cuda_Velocity_Verlet_NVE( reax_system*, control_params*,
        simulation_data*, storage*, reax_list**, output_controls*,
        mpi_datatypes* );

int Cuda_Velocity_Verlet_Nose_Hoover_NVT_Klein( reax_system*, control_params*,
        simulation_data*, storage*, reax_list**, output_controls*,
        mpi_datatypes* );

int Cuda_Velocity_Verlet_Berendsen_NVT( reax_system*, control_params*,
        simulation_data*, storage*, reax_list**, output_controls*,
        mpi_datatypes* );

int Cuda_Velocity_Verlet_Berendsen_NPT( reax_system*, control_params*,
        simulation_data*, storage*, reax_list**, output_controls*,
        mpi_datatypes* );

#ifdef __cplusplus
}
#endif


#endif
