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

#ifndef __BONDS_H_
#define __BONDS_H_

#include "reaxc_types.h"

void Bonds( reax_system*, control_params*, simulation_data*, 
	    storage*, reax_list**, output_controls* );
#endif
