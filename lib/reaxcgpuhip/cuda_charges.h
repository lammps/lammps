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

#ifndef __CUDA_CHARGES_H_
#define __CUDA_CHARGES_H_

#include "reaxc_types.h"

#ifdef __cplusplus
extern "C" {
#endif


void Cuda_Init_MatVec( reax_system *, storage * );

void cuda_charges_x( reax_system *, rvec2 );

void cuda_charges_st( reax_system *, storage *, real *, real );

void cuda_charges_updateq( reax_system *, real * );

void Cuda_Compute_Charges( reax_system*, control_params*, simulation_data*,
        storage*, output_controls*, mpi_datatypes* );


#ifdef __cplusplus
}
#endif

#endif
