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

#ifndef __CUDA_VECTOR_H_
#define __CUDA_VECTOR_H_

#include "reaxc_types.h"

#include "cuda_random.h"


CUDA_DEVICE static inline void cuda_rvec_Random( rvec v )
{
//    v[0] = Cuda_Random( 2.0 ) - 1.0;
//    v[1] = Cuda_Random( 2.0 ) - 1.0;
//    v[2] = Cuda_Random( 2.0 ) - 1.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}


#endif
