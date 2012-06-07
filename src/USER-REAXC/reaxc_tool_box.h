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

#ifndef __TOOL_BOX_H_
#define __TOOL_BOX_H_

#include "reaxc_types.h"
#include "reaxc_defs.h"

/* from comm_tools.h */
int SumScan( int, int, int, MPI_Comm );
void SumScanB( int, int, int, int, MPI_Comm, int* );

/* from box.h */
void Transform_to_UnitBox( rvec, simulation_box*, char, rvec );
void Fit_to_Periodic_Box( simulation_box*, rvec* );
void Box_Touch_Point( simulation_box*, ivec, rvec );
int  is_Inside_Box( simulation_box*, rvec );
int  iown_midpoint( simulation_box*, rvec, rvec );

/* from grid.h */
void GridCell_Closest_Point( grid_cell*, grid_cell*, ivec, ivec, rvec );
void GridCell_to_Box_Points( grid_cell*, ivec, rvec, rvec );
real DistSqr_between_Special_Points( rvec, rvec );
real DistSqr_to_Special_Point( rvec, rvec );
int Relative_Coord_Encoding( ivec );

/* from geo_tools.h */
void Make_Point( real, real, real, rvec* );
int is_Valid_Serial( storage*, int );
int Check_Input_Range( int, int, int, char*, MPI_Comm );
void Trim_Spaces( char* );

/* from system_props.h */
real Get_Time( );
real Get_Timing_Info( real );
void Update_Timing_Info( real*, real* );

/* from io_tools.h */
int   Get_Atom_Type( reax_interaction*, char*, MPI_Comm );
char *Get_Element( reax_system*, int );
char *Get_Atom_Name( reax_system*, int );
int   Allocate_Tokenizer_Space( char**, char**, char*** );
int   Tokenize( char*, char*** );

/* from lammps */
void *smalloc( long, char*, MPI_Comm );
void *scalloc( int, int, char*, MPI_Comm );
void sfree( void*, char* );

#endif
