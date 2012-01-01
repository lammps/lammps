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

#ifndef __BOND_ORDERS_H_
#define __BOND_ORDERS_H_

#include "reaxc_types.h"

typedef struct{
  real C1dbo, C2dbo, C3dbo;
  real C1dbopi, C2dbopi, C3dbopi, C4dbopi;
  real C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
  real C1dDelta, C2dDelta, C3dDelta;
}dbond_coefficients;

#ifdef TEST_FORCES
void Get_dBO( reax_system*, reax_list**, int, int, real, rvec* );
void Get_dBOpinpi2( reax_system*, reax_list**, 
		    int, int, real, real, rvec*, rvec* );

void Add_dBO( reax_system*, reax_list**, int, int, real, rvec* );
void Add_dBOpinpi2( reax_system*, reax_list**, 
		    int, int, real, real, rvec*, rvec* );

void Add_dBO_to_Forces( reax_system*, reax_list**, int, int, real );
void Add_dBOpinpi2_to_Forces( reax_system*, reax_list**, 
			      int, int, real, real );

void Add_dDelta( reax_system*, reax_list**, int, real, rvec* );
void Add_dDelta_to_Forces( reax_system *, reax_list**, int, real );
#endif

void Add_dBond_to_Forces( reax_system*, int, int, storage*, reax_list** );
void Add_dBond_to_Forces_NPT( int, int, simulation_data*, 
			      storage*, reax_list** );
int BOp(storage*, reax_list*, real, int, int, far_neighbor_data*,
	single_body_parameters*, single_body_parameters*, two_body_parameters*);
void BO( reax_system*, control_params*, simulation_data*,
	 storage*, reax_list**, output_controls* );
#endif
