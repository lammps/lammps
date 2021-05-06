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

#ifndef __CUDA_MULTI_BODY_H_
#define __CUDA_MULTI_BODY_H_

#include "reaxc_types.h"


CUDA_GLOBAL void Cuda_Atom_Energy( reax_atom *, global_parameters,
        single_body_parameters *, two_body_parameters *, storage,
        reax_list, int, int, real *, real *, real *);

CUDA_GLOBAL void Cuda_Atom_Energy_PostProcess( reax_list, storage, int );


#endif
