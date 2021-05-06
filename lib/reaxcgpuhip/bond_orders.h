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

#ifndef __BOND_ORDERS_H_
#define __BOND_ORDERS_H_

#include "reaxc_types.h"


typedef struct
{
    real C1dbo;
    real C2dbo;
    real C3dbo;
    real C1dbopi;
    real C2dbopi;
    real C3dbopi;
    real C4dbopi;
    real C1dbopi2;
    real C2dbopi2;
    real C3dbopi2;
    real C4dbopi2;
    real C1dDelta;
    real C2dDelta;
    real C3dDelta;
} dbond_coefficients;


#ifdef __cplusplus
extern "C" {
#endif



#ifdef __cplusplus
}
#endif


#endif
