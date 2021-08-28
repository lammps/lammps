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
   Contributing authors: Mario Orsi & Wei Ding (QMUL), m.orsi@qmul.ac.uk
------------------------------------------------------------------------- */

#include "angle_dipole.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AngleDipole::AngleDipole(LAMMPS *lmp) : Angle(lmp)
{
  k = nullptr;
  gamma0 = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleDipole::~AngleDipole()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(gamma0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleDipole::compute(int eflag, int vflag)
{
  int iRef, iDip, iDummy, n, type;
  double delx, dely, delz;
  double eangle, tangle, fi[3], fj[3];
  double r, cosGamma, deltaGamma, kdg, rmu;

  eangle = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;      // position vector
  double **mu = atom->mu;    // point-dipole components and moment magnitude
  double **torque = atom->torque;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **f = atom->f;
  double delTx, delTy, delTz;
  double fx, fy, fz, fmod, fmod_sqrtff;

  if (!newton_bond) error->all(FLERR, "'newton' flag for bonded interactions must be 'on'");

  for (n = 0; n < nanglelist; n++) {
    iDip = anglelist[n][0];      // dipole whose orientation is to be restrained
    iRef = anglelist[n][1];      // reference atom toward which dipole will point
    iDummy = anglelist[n][2];    // dummy atom - irrelevant to the interaction
    type = anglelist[n][3];

    delx = x[iRef][0] - x[iDip][0];
    dely = x[iRef][1] - x[iDip][1];
    delz = x[iRef][2] - x[iDip][2];

    r = sqrt(delx * delx + dely * dely + delz * delz);

    rmu = r * mu[iDip][3];
    cosGamma = (mu[iDip][0] * delx + mu[iDip][1] * dely + mu[iDip][2] * delz) / rmu;
    deltaGamma = cosGamma - cos(gamma0[type]);
    kdg = k[type] * deltaGamma;

    if (eflag) eangle = kdg * deltaGamma;    // energy

    tangle = 2.0 * kdg / rmu;

    delTx = tangle * (dely * mu[iDip][2] - delz * mu[iDip][1]);
    delTy = tangle * (delz * mu[iDip][0] - delx * mu[iDip][2]);
    delTz = tangle * (delx * mu[iDip][1] - dely * mu[iDip][0]);

    torque[iDip][0] += delTx;
    torque[iDip][1] += delTy;
    torque[iDip][2] += delTz;

    // Force couple that counterbalances dipolar torque
    fx = dely * delTz - delz * delTy;    // direction (fi): - r x (-T)
    fy = delz * delTx - delx * delTz;
    fz = delx * delTy - dely * delTx;

    fmod = sqrt(delTx * delTx + delTy * delTy + delTz * delTz) / r;    // magnitude
    fmod_sqrtff = fmod / sqrt(fx * fx + fy * fy + fz * fz);

    fi[0] = fx * fmod_sqrtff;
    fi[1] = fy * fmod_sqrtff;
    fi[2] = fz * fmod_sqrtff;

    fj[0] = -fi[0];
    fj[1] = -fi[1];
    fj[2] = -fi[2];

    f[iDip][0] += fj[0];
    f[iDip][1] += fj[1];
    f[iDip][2] += fj[2];

    f[iRef][0] += fi[0];
    f[iRef][1] += fi[1];
    f[iRef][2] += fi[2];

    if (evflag)    // virial = rij.fi = 0 (fj = -fi & fk = 0)
      ev_tally(iRef, iDip, iDummy, nlocal, newton_bond, eangle, fj, fi, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleDipole::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k, n + 1, "angle:k");
  memory->create(gamma0, n + 1, "angle:gamma0");

  memory->create(setflag, n + 1, "angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleDipole::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double gamma0_one = utils::numeric(FLERR, arg[2], false, lmp);

  // convert gamma0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    gamma0[i] = gamma0_one / 180.0 * MY_PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleDipole::init_style()
{
  if (!atom->mu_flag || !atom->torque_flag)
    error->all(FLERR, "Angle style dipole requires atom attributes mu and torque");
}

/* ----------------------------------------------------------------------
   used by SHAKE
------------------------------------------------------------------------- */

double AngleDipole::equilibrium_angle(int i)
{
  return gamma0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleDipole::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&gamma0[1], sizeof(double), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleDipole::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &gamma0[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&gamma0[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleDipole::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++) fprintf(fp, "%d %g %g\n", i, k[i], gamma0[i]);
}

/* ----------------------------------------------------------------------
   used by ComputeAngleLocal
------------------------------------------------------------------------- */

double AngleDipole::single(int type, int iRef, int iDip, int /*iDummy*/)
{
  double **x = atom->x;      // position vector
  double **mu = atom->mu;    // point-dipole components and moment magnitude

  double delx = x[iRef][0] - x[iDip][0];
  double dely = x[iRef][1] - x[iDip][1];
  double delz = x[iRef][2] - x[iDip][2];

  domain->minimum_image(delx, dely, delz);

  double r = sqrt(delx * delx + dely * dely + delz * delz);
  double rmu = r * mu[iDip][3];
  double cosGamma = (mu[iDip][0] * delx + mu[iDip][1] * dely + mu[iDip][2] * delz) / rmu;
  double deltaGamma = cosGamma - cos(gamma0[type]);
  double kdg = k[type] * deltaGamma;

  return kdg * deltaGamma;    // energy
}
