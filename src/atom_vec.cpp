/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "atom_vec.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVec::AtomVec(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  nmax = 0;
  molecular = 0;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 0;
  mass_type = dipole_type = 0;
  comm_x_only = comm_f_only = 1;
  size_comm = size_reverse = size_border = 0;
}

/* ----------------------------------------------------------------------
   unpack a single line from Velocity section of data file
   individual style may override this
------------------------------------------------------------------------- */

void AtomVec::data_vel(int m, char *line, int ihybrid)
{
  int tmp;
  double **v = atom->v;
  sscanf(line,"%d %lg %lg %lg",&tmp,&v[m][0],&v[m][1],&v[m][2]);
}
