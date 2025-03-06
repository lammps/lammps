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
   Contributing author: James D. Farrell (IoP CAS) j.d.farrell at gmail.com
   [ based on angle_cosine.cpp ]
------------------------------------------------------------------------- */

#include "angle_mwlc.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

AngleMWLC::AngleMWLC(LAMMPS *_lmp) : Angle(_lmp)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

AngleMWLC::~AngleMWLC()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k1);
    memory->destroy(k2);
    memory->destroy(mu);
    memory->destroy(temp);
  }
}

/* ---------------------------------------------------------------------- */

void AngleMWLC::compute(int eflag, int vflag)
{
  int i1, i2, i3, n, type;
  double delx1, dely1, delz1, delx2, dely2, delz2;
  double eangle, f1[3], f3[3];
  double rsq1, rsq2, r1, r2, c, a, a11, a12, a22;
  double q, qm, Q;

  eangle = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double kbt, v_min;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    kbt = temp[type] * force->boltz;
    v_min = -kbt * log(1 + exp(-mu[type] / kbt));

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    r2 = sqrt(rsq2);

    // c = cosine of angle

    c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy

    q = exp(-k1[type] * (1.0 + c) / kbt);
    qm = exp((-k2[type] * (1.0 + c) - mu[type]) / kbt);
    Q = q + qm;

    if (eflag) eangle = -kbt * log(Q) - v_min;

    a = (k1[type] * q + k2[type] * qm) / Q;
    a11 = a * c / rsq1;
    a12 = -a / (r1 * r2);
    a22 = a * c / rsq2;

    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag)
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, f1, f3, delx1, dely1, delz1, delx2, dely2,
               delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleMWLC::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(k1, np1, "angle:k1");
  memory->create(k2, np1, "angle:k2");
  memory->create(mu, np1, "angle:mu");
  memory->create(temp, np1, "angle:temp");
  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void AngleMWLC::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR, "Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double k1_one = utils::numeric(FLERR, arg[1], false, lmp);
  double k2_one = utils::numeric(FLERR, arg[2], false, lmp);
  double mu_one = utils::numeric(FLERR, arg[3], false, lmp);
  double temp_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k1[i] = k1_one;
    k2[i] = k2_one;
    mu[i] = mu_one;
    temp[i] = temp_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleMWLC::equilibrium_angle(int /*i*/)
{
  return MY_PI;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleMWLC::write_restart(FILE *fp)
{
  fwrite(&k1[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&k2[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&mu[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&temp[1], sizeof(double), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleMWLC::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k1[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &k2[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &mu[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &temp[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&k1[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&k2[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mu[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&temp[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleMWLC::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %g %g %g %g\n", i, k1[i], k2[i], mu[i], temp[i]);
}

/* ---------------------------------------------------------------------- */

double AngleMWLC::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double kbt = temp[type] * force->boltz;
  double v_min = -kbt * log(1 + exp(-mu[type] / kbt));
  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1, dely1, delz1);
  double r1 = sqrt(delx1 * delx1 + dely1 * dely1 + delz1 * delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2, dely2, delz2);
  double r2 = sqrt(delx2 * delx2 + dely2 * dely2 + delz2 * delz2);

  double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
  c /= r1 * r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double q = exp(-k1[type] * (1.0 + c) / kbt);
  double qm = exp((-k2[type] * (1.0 + c) - mu[type]) / kbt);
  return -kbt * log(q + qm) - v_min;
}

/* ---------------------------------------------------------------------- */

void AngleMWLC::born_matrix(int type, int i1, int i2, int i3, double &du, double &du2)
{
  double **x = atom->x;
  double kbt = temp[type] * force->boltz;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1, dely1, delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2, dely2, delz2);

  double c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
  c /= sqrt((delx1 * delx1 + dely1 * dely1 + delz1 * delz1) *
            (delx2 * delx2 + dely2 * dely2 + delz2 * delz2));
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  const double q = exp(-k1[type] * (1.0 + c) / kbt);
  const double qm = exp((-k2[type] * (1.0 + c) - mu[type]) / kbt);
  const double Q = q + qm;

  du = (k1[type] * q + k2[type] * qm) / Q;
  du2 = (k1[type] - k2[type]) / Q;
  du2 *= -du2 * q * qm / kbt;
}

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *AngleMWLC::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "k1") == 0) return (void *) k1;
  if (strcmp(str, "k2") == 0) return (void *) k2;
  if (strcmp(str, "mu") == 0) return (void *) mu;
  if (strcmp(str, "temp") == 0) return (void *) temp;
  return nullptr;
}
