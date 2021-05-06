/*----------------------------------------------------------------------
 *   PuReMD - Purdue ReaxFF Molecular Dynamics Program
 *
 *   Copyright (2010) Purdue University
 *   Hasan Metin Aktulga, haktulga@cs.purdue.edu
 *   Joseph Fogarty, jcfogart@mail.usf.edu
 *   Sagar Pandit, pandit@usf.edu
 *   Ananth Y Grama, ayg@cs.purdue.edu
 *
 *   This program is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU General Public License as
 *   published by the Free Software Foundation; either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details:
 *   <http://www.gnu.org/licenses/>.
 *----------------------------------------------------------------------*/

#ifndef __CUDA_SHUFFLE_H_
#define __CUDA_SHUFFLE_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

#if defined(__SM_35__)

/* Part of the code is taken from this site.
 * And the other is taken from the download in the PGPuReMD folder on CUPID
 * http://wenda.baba.io/questions/4481817/overloading-the-cuda-shuffle-function-makes-the-original-ones-invisible.html
 */
CUDA_DEVICE static inline real shfl(real x, int lane)
{
    // Split the double number into 2 32b registers.
    int lo, hi;
    asm volatile( "mov.b64 {%0,%1}, %2;" : "=r"(lo), "=r"(hi) : "d"(x) );

    // Shuffle the two 32b registers.
    lo = __shfl_xor( lo, lane );
    hi = __shfl_xor( hi, lane );

    // Recreate the 64b number.
    //asm volatile( "mov.b64 %0, {%1,%2};" : "=d(x)" : "r"(lo), "r"(hi) );
    //return x;
    return __hiloint2double( hi, lo );
}

#endif

#ifdef __cplusplus
}
#endif


#endif
