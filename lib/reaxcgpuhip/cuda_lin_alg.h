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

#ifndef __CUDA_LIN_ALG_H_
#define __CUDA_LIN_ALG_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void Cuda_Vector_Sum( real *, real, real *, real, real *, int );

void Cuda_CG_Preconditioner( real *, real *, real *, int );

void Cuda_CG_Diagonal_Preconditioner( storage *, rvec2 *, int );

void Cuda_DualCG_Preconditioner( storage *, rvec2 *, rvec2, int, rvec2 );

void Cuda_Norm( rvec2 *, int, rvec2 );

void Cuda_Dot( rvec2 *, rvec2 *, rvec2, int );

void Cuda_Vector_Sum_Rvec2( rvec2 *, rvec2 *, rvec2, rvec2 *, int );

void Cuda_RvecCopy_From( real *, rvec2 *, int, int );

void Cuda_RvecCopy_To( rvec2 *, real *, int, int );

void Cuda_Dual_Matvec( sparse_matrix *, rvec2 *, rvec2 *, int , int );

void Cuda_Matvec( sparse_matrix *, real *, real *, int , int );

int Cuda_dual_CG( reax_system*, control_params*, storage*, sparse_matrix*,
        rvec2*, real, rvec2*, mpi_datatypes*, FILE* , simulation_data * );

int Cuda_CG( reax_system*, control_params*, storage*, sparse_matrix*,
        real*, real, real*, mpi_datatypes* );

#ifdef __cplusplus
}
#endif


#endif
