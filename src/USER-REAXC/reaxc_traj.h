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

#ifndef __TRAJ_H__
#define __TRAJ_H__

#include "reaxc_types.h"

#define MAX_TRJ_LINE_LEN     120
#define MAX_TRJ_BUFFER_SIZE  (MAX_TRJ_LINE_LEN * 100)

#define NUM_HEADER_LINES 37
#define HEADER_LINE_LEN 62
#define STR_LINE  "%-37s%-24s\n"
#define INT_LINE  "%-37s%-24d\n"
#define INT2_LINE  "%-36s%-12d,%-12d\n"
#define REAL_LINE "%-37s%-24.3f\n"
#define SCI_LINE  "%-37s%-24g\n"
#define REAL3_LINE "%-32s%9.3f,%9.3f,%9.3f\n"

#define INIT_DESC "%9d%3d%9s%10.3f\n" //AtomID - AtomType, AtomName, AtomMass
#define INIT_DESC_LEN 32

#define SIZE_INFO_LINE2 "%-10d%-10d\n"
#define SIZE_INFO_LEN2 21

#define SIZE_INFO_LINE3 "%-10d%-10d%-10d\n"
#define SIZE_INFO_LEN3 31

#define ATOM_BASIC "%9d%10.3f%10.3f%10.3f%10.3f\n" //AtomID AtomType (X Y Z) Charge
#define ATOM_BASIC_LEN 50
#define ATOM_wV    "%9d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" //AtomID (X Y Z) (Vx Vy Vz) Charge
#define ATOM_wV_LEN 80
#define ATOM_wF    "%9d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" //AtomID (X Y Z) (Fx Fy Fz) Charge
#define ATOM_wF_LEN 80
#define ATOM_FULL  "%9d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" //AtomID (X Y Z) (Vx Vy Vz) (Fx Fy Fz) Charge
#define ATOM_FULL_LEN 110

#define BOND_BASIC "%9d%9d%10.3f%10.3f\n" // Atom1 Atom2 Dist Total_BO
#define BOND_BASIC_LEN 39
#define BOND_FULL  "%9d%9d%10.3f%10.3f%10.3f%10.3f%10.3f\n" // Atom1 Atom2 Dist Total_BO BOs BOpi BOpi2
#define BOND_FULL_LEN 69

#define ANGLE_BASIC "%9d%9d%9d%10.3f\n" // Atom1 Atom2 Atom3 Theta
#define ANGLE_BASIC_LEN 38

enum ATOM_LINE_OPTS  { OPT_NOATOM = 0, OPT_ATOM_BASIC = 4, OPT_ATOM_wF = 5, OPT_ATOM_wV = 6, OPT_ATOM_FULL = 7, NR_OPT_ATOM = 8 };
enum BOND_LINE_OPTS  { OPT_NOBOND, OPT_BOND_BASIC, OPT_BOND_FULL, NR_OPT_BOND };
enum ANGLE_LINE_OPTS { OPT_NOANGLE, OPT_ANGLE_BASIC, NR_OPT_ANGLE };


int  Init_Traj( reax_system*, control_params*, output_controls*, 
		mpi_datatypes*, char* );
int  End_Traj( int, output_controls* );

int  Append_Frame( reax_system*, control_params*, simulation_data*, 
		   reax_list**, output_controls*, mpi_datatypes* );

#endif
