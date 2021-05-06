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

#ifndef __CUDA_TORSION_ANGLES_H_
#define __CUDA_TORSION_ANGLES_H_

#include "reaxc_types.h"


CUDA_GLOBAL void Cuda_Torsion_Angles( reax_atom *, global_parameters,
        four_body_header *, control_params *, reax_list, reax_list,
        storage, int, int, real *, real *, rvec * );

CUDA_GLOBAL void Cuda_Torsion_Angles_PostProcess( reax_atom *,
        storage, reax_list, int );


#endif
