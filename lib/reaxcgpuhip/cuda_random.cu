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

#include "cuda_random.h"


/* System random number generator used linear congruance method with
   large periodicity for generation of pseudo random number. function
   Random returns this random number appropriately scaled so that
   0 <= Random(range) < range */
CUDA_DEVICE double Cuda_Random( double range )
{
    //TODO: use cuRAND
//    return (random( ) * range) / 2147483647L;
    return 0.0;
}


/* This function seeds the system pseudo random number generator with
   current time. Use this function once in the begining to initialize
   the system */
void Cuda_Randomize( )
{
    //TODO: use cuRAND
//    hiprandState_t state;
//
//    hiprand_init( time(NULL), 0, 0, &state );
}


/* GRandom return random number with gaussian distribution with mean
   and standard deviation "sigma" */
CUDA_DEVICE double Cuda_GRandom( double mean, double sigma )
{
    double v1 = Cuda_Random(2.0) - 1.0;
    double v2 = Cuda_Random(2.0) - 1.0;
    double rsq = v1 * v1 + v2 * v2;

    while (rsq >= 1.0 || rsq == 0.0)
    {
        v1 = Cuda_Random(2.0) - 1.0;
        v2 = Cuda_Random(2.0) - 1.0;
        rsq = v1 * v1 + v2 * v2;
    }

    return mean + v1 * sigma * SQRT(-2.0 * LOG(rsq) / rsq);
}
