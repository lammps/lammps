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
#include "reaxc_basic_comm.h"
#include "reaxc_vector.h"

real Parallel_Norm( real *v, int n, MPI_Comm comm )
{
  int  i;
  real my_sum, norm_sqr;

  my_sum = 0;
  for( i = 0; i < n; ++i )
    my_sum += SQR( v[i] );

  MPI_Allreduce( &my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, comm );

  return sqrt( norm_sqr );
}

real Parallel_Dot( real *v1, real *v2, int n, MPI_Comm comm )
{
  int  i;
  real my_dot, res;

  my_dot = 0;
  for( i = 0; i < n; ++i )
    my_dot += v1[i] * v2[i];

  MPI_Allreduce( &my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, comm );

  return res;
}

real Parallel_Vector_Acc( real *v, int n, MPI_Comm comm )
{
  int  i;
  real my_acc, res;

  my_acc = 0;
  for( i = 0; i < n; ++i )
    my_acc += v[i];

  MPI_Allreduce( &my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, comm );

  return res;
}
