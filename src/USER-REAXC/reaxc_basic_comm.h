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

#ifndef __BASIC_COMM_H_
#define __BASIC_COMM_H_

#include "reaxc_types.h"

void real_packer( void*, mpi_out_data* );
void rvec_packer( void*, mpi_out_data* );
void rvec2_packer( void*, mpi_out_data* );
void Dist(reax_system*, mpi_datatypes*, void*, MPI_Datatype, int, dist_packer);

void real_unpacker( void*, void*, mpi_out_data* );
void rvec_unpacker( void*, void*, mpi_out_data* );
void rvec2_unpacker( void*, void*, mpi_out_data* );
void Coll( reax_system*, mpi_datatypes*, void*, MPI_Datatype,
           int, coll_unpacker );

real Parallel_Norm( real*, int, MPI_Comm );
real Parallel_Dot( real*, real*, int, MPI_Comm );
real Parallel_Vector_Acc( real*, int, MPI_Comm );

#if defined(TEST_FORCES)
void Coll_ids_at_Master( reax_system*, storage*, mpi_datatypes* );
void Coll_rvecs_at_Master( reax_system*, storage*, mpi_datatypes*, rvec* );
#endif

#endif
