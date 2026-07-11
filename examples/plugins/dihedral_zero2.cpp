/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Carsten Svaneborg (SDU)
------------------------------------------------------------------------- */

#include "dihedral_zero2.h"

#include "atom.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DihedralZero2::DihedralZero2(LAMMPS *lmp) : Dihedral(lmp), coeffflag(1)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralZero2::~DihedralZero2()
{
  if (allocated && !copymode) { memory->destroy(setflag); }
}

/* ---------------------------------------------------------------------- */

void DihedralZero2::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/* ---------------------------------------------------------------------- */

void DihedralZero2::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1)) error->all(FLERR, "Illegal dihedral_style command");

  if (narg == 1) {
    if (strcmp("nocoeff", arg[0]) == 0)
      coeffflag = 0;
    else
      error->all(FLERR, "Illegal dihedral_style command");
  }
}

/* ---------------------------------------------------------------------- */

void DihedralZero2::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(setflag, n + 1, "dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void DihedralZero2::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 1))
    error->all(FLERR, "Incorrect args for dihedral coefficients");

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralZero2::write_restart(FILE * /*fp*/) {}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralZero2::read_restart(FILE * /*fp*/)
{
  allocate();
  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralZero2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++) fprintf(fp, "%d\n", i);
}
