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

#include "bond_zero2.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondZero2::BondZero2(LAMMPS *lmp) : Bond(lmp), coeffflag(1) {}

/* ---------------------------------------------------------------------- */

BondZero2::~BondZero2()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondZero2::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/* ---------------------------------------------------------------------- */

void BondZero2::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1)) error->all(FLERR, "Illegal bond_style command");

  if (narg == 1) {
    if (strcmp("nocoeff", arg[0]) == 0)
      coeffflag = 0;
    else
      error->all(FLERR, "Illegal bond_style command");
  }
}

/* ---------------------------------------------------------------------- */

void BondZero2::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(r0, n + 1, "bond:r0");
  memory->create(setflag, n + 1, "bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondZero2::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 2))
    error->all(FLERR, "Incorrect args for bond coefficients");

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double r0_one = 0.0;
  if (coeffflag && (narg == 2)) r0_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    r0[i] = r0_one;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondZero2::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondZero2::write_restart(FILE *fp)
{
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondZero2::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondZero2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) fprintf(fp, "%d %g\n", i, r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondZero2::single(int /*type*/, double /*rsq*/, int /*i*/, int /*j*/, double & /*fforce*/)
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void *BondZero2::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "r0") == 0) return (void *) r0;
  return nullptr;
}
