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

#ifndef __RANDOM_H_
#define __RANDOM_H_

#include "reaxc_types.h"


#ifdef _cplusplus
extern "C" {
#endif

/* System random number generator used linear congruance method with
   large periodicity for generation of pseudo random number. function
   Random returns this random number appropriately scaled so that
   0 <= Random(range) < range */
double Random( double );

/* This function seeds the system pseudo random number generator with
   current time. Use this function once in the begining to initialize
   the system */
void Randomize( );

/* GRandom return random number with gaussian distribution with mean
   and standard deviation "sigma" */
double GRandom( double, double );

#ifdef _cplusplus
}
#endif


#endif
