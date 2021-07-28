/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "improper_zero2.h"

#include "atom.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ImproperZero2::ImproperZero2(LAMMPS *lmp) : Improper(lmp), coeffflag(1)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperZero2::~ImproperZero2()
{
  if (allocated && !copymode) { memory->destroy(setflag); }
}

/* ---------------------------------------------------------------------- */

void ImproperZero2::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/* ---------------------------------------------------------------------- */

void ImproperZero2::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1)) error->all(FLERR, "Illegal improper_style command");

  if (narg == 1) {
    if (strcmp("nocoeff", arg[0]) == 0)
      coeffflag = 0;
    else
      error->all(FLERR, "Illegal improper_style command");
  }
}

/* ---------------------------------------------------------------------- */

void ImproperZero2::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(setflag, n + 1, "improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void ImproperZero2::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 1))
    error->all(FLERR, "Incorrect args for improper coefficients");

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nimpropertypes, ilo, ihi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperZero2::write_restart(FILE * /*fp*/) {}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperZero2::read_restart(FILE * /*fp*/)
{
  allocate();
  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperZero2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++) fprintf(fp, "%d\n", i);
}
