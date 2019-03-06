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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "improper_zero.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ImproperZero::ImproperZero(LAMMPS *lmp) : Improper(lmp), coeffflag(1) {}

/* ---------------------------------------------------------------------- */

ImproperZero::~ImproperZero()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperZero::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;
}

/* ---------------------------------------------------------------------- */

void ImproperZero::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1))
    error->all(FLERR,"Illegal improper_style command");

  if (narg == 1) {
    if (strcmp("nocoeff",arg[0]) == 0) coeffflag=0;
    else error->all(FLERR,"Illegal improper_style command");
  }
}

/* ---------------------------------------------------------------------- */

void ImproperZero::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void ImproperZero::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 1))
    error->all(FLERR,"Incorrect args for improper coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperZero::write_restart(FILE * /*fp*/) {}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperZero::read_restart(FILE * /*fp*/)
{
  allocate();
  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperZero::write_data(FILE *fp) {
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d\n",i);
}

