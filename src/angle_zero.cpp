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

#include "angle_zero.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

/* ---------------------------------------------------------------------- */

AngleZero::AngleZero(LAMMPS *_lmp) : Angle(_lmp), coeffflag(1) {}

/* ---------------------------------------------------------------------- */

AngleZero::~AngleZero()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(theta0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleZero::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/* ---------------------------------------------------------------------- */

void AngleZero::settings(int narg, char **arg)
{
  if ((narg != 0) && (narg != 1)) error->all(FLERR, "Illegal angle_style command");

  if (narg == 1) {
    if (strcmp("nocoeff", arg[0]) == 0)
      coeffflag = 0;
    else
      error->all(FLERR, "Illegal angle_style command");
  }
}

/* ---------------------------------------------------------------------- */

void AngleZero::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(theta0, np1, "angle:theta0");
  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleZero::coeff(int narg, char **arg)
{
  if ((narg < 1) || (coeffflag && narg > 2))
    error->all(FLERR, "Incorrect args for angle coefficients");

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double theta0_one = 0.0;
  if (coeffflag && (narg == 2)) theta0_one = utils::numeric(FLERR, arg[1], false, lmp);

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    setflag[i] = 1;
    theta0[i] = DEG2RAD * theta0_one;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleZero::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleZero::write_restart(FILE *fp)
{
  fwrite(&theta0[1], sizeof(double), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleZero::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &theta0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&theta0[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}
/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleZero::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++) fprintf(fp, "%d %g\n", i, RAD2DEG * theta0[i]);
}

/* ---------------------------------------------------------------------- */

double AngleZero::single(int /*type*/, int /*i1*/, int /*i2*/, int /*i3*/)
{
  return 0.0;
}
