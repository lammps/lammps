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

/* ----------------------------------------------------------------------
   Contributing author: Carsten Svaneborg (SDU)
------------------------------------------------------------------------- */

#include <cstring>
#include "dihedral_zero.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DihedralZero::DihedralZero(LAMMPS *lmp) : Dihedral(lmp), coeffflag(1) {}

/* ---------------------------------------------------------------------- */

DihedralZero::~DihedralZero()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralZero::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);
}

/* ---------------------------------------------------------------------- */

void DihedralZero::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1))
    error->all(FLERR,"Illegal dihedral_style command");

  if (narg == 1) {
    if (strcmp("nocoeff",arg[0]) == 0) coeffflag=0;
    else error->all(FLERR,"Illegal dihedral_style command");
  }
}

/* ---------------------------------------------------------------------- */

void DihedralZero::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void DihedralZero::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 1))
    error->all(FLERR,"Incorrect args for dihedral coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralZero::write_restart(FILE * /*fp*/) {}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralZero::read_restart(FILE * /*fp*/)
{
  allocate();
  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralZero::write_data(FILE *fp) {
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp,"%d\n",i);
}
