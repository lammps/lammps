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
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#include "angle_mesocnt.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_2PI;
using MathConst::MY_PI;
using MathConst::RAD2DEG;

static constexpr double SMALL = 0.001;
static constexpr double A_CC = 1.421;

/* ---------------------------------------------------------------------- */

AngleMesoCNT::AngleMesoCNT(LAMMPS *_lmp) : Angle(_lmp)
{
  buckling = nullptr;
  kh = nullptr;
  kb = nullptr;
  thetab = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleMesoCNT::~AngleMesoCNT()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(buckling);
    memory->destroy(kh);
    memory->destroy(kb);
    memory->destroy(thetab);
    memory->destroy(theta0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleMesoCNT::compute(int eflag, int vflag)
{
  int i1, i2, i3, n, type;
  double delx1, dely1, delz1, delx2, dely2, delz2;
  double eangle, f1[3], f3[3];
  double dtheta, tk;
  double rsq1, rsq2, r1, r2, c, s, a, a11, a12, a22;

  eangle = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int flag, cols;
  int index = atom->find_custom("buckled", flag, cols);
  int *buckled = atom->ivector[index];

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

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

    // angle (cos and sin)

    c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
    c /= r1 * r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    // force & energy

    dtheta = acos(c) - theta0[type];

    // harmonic bending
    if (!buckling[type] || fabs(dtheta) < thetab[type]) {
      tk = kh[type] * dtheta;
      if (eflag) eangle = tk * dtheta;
      a = -2.0 * tk * s;

      buckled[i2] = 0;
    }
    // bending buckling
    else {
      if (eflag)
        eangle = kb[type] * fabs(dtheta) + thetab[type] * (kh[type] * thetab[type] - kb[type]);
      a = kb[type] * s;

      buckled[i2] = 1;
    }
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

void AngleMesoCNT::allocate()
{
  allocated = 1;
  const int np1 = atom->nangletypes + 1;

  memory->create(buckling, np1, "angle:buckling");
  memory->create(kh, np1, "angle:kh");
  memory->create(kb, np1, "angle:kb");
  memory->create(thetab, np1, "angle:thetab");
  memory->create(theta0, np1, "angle:theta0");

  memory->create(setflag, np1, "angle:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleMesoCNT::coeff(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "angle_coeff", error);

  int buckling_one;
  if (strcmp(arg[1], "buckling") == 0)
    buckling_one = 1;
  else if (strcmp(arg[1], "harmonic") == 0)
    buckling_one = 0;
  else
    error->all(FLERR,
               "Unknown first argument {} for angle coefficients, must be 'buckling' or 'harmonic'",
               arg[1]);

  // units, eV to energy unit conversion

  double ang = force->angstrom;
  double eunit;
  if (strcmp(update->unit_style, "real") == 0)
    eunit = 23.06054966;
  else if (strcmp(update->unit_style, "metal") == 0)
    eunit = 1.0;
  else if (strcmp(update->unit_style, "si") == 0)
    eunit = 1.6021765e-19;
  else if (strcmp(update->unit_style, "cgs") == 0)
    eunit = 1.6021765e-12;
  else if (strcmp(update->unit_style, "electron") == 0)
    eunit = 3.674932248e-2;
  else if (strcmp(update->unit_style, "micro") == 0)
    eunit = 1.6021765e-4;
  else if (strcmp(update->unit_style, "nano") == 0)
    eunit = 1.6021765e2;
  else
    error->all(FLERR, "Angle style mesocnt does not support {} units", update->unit_style);

  // set parameters

  double kh_one, kb_one, thetab_one;
  if (strcmp(arg[2], "custom") == 0) {
    if (buckling_one) {
      if (narg != 6) error->all(FLERR, "Incorrect number of args for 'custom' angle coefficients");
      kb_one = utils::numeric(FLERR, arg[4], false, lmp);
      thetab_one = utils::numeric(FLERR, arg[5], false, lmp);
    } else if (narg != 4)
      error->all(FLERR, "Incorrect number of args for 'custom' angle coefficients");

    kh_one = utils::numeric(FLERR, arg[3], false, lmp);
  } else if (strcmp(arg[2], "C") == 0) {
    if (narg != 6)
      error->all(FLERR, "Incorrect number of args for 'C' preset in angle coefficients");
    int n = utils::inumeric(FLERR, arg[3], false, lmp);
    int m = utils::inumeric(FLERR, arg[4], false, lmp);
    double l = utils::numeric(FLERR, arg[5], false, lmp);

    double r_ang = sqrt(3.0 * (n * n + n * m + m * m)) * A_CC / MY_2PI;

    // empirical parameters

    double k = 63.80 * pow(r_ang, 2.93) * eunit * ang;
    kh_one = 0.5 * k / l;

    if (buckling_one) {
      kb_one = 0.7 * k / (275.0 * ang);
      thetab_one = 180.0 / MY_PI * atan(l / (275.0 * ang));
    }
  } else
    error->all(FLERR, "Unknown {} preset in angle coefficients", arg[2]);

  // set safe default values for buckling parameters if buckling is disabled

  if (!buckling_one) {
    kb_one = 0.0;
    thetab_one = 180.0;
  }

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  // convert thetab from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    buckling[i] = buckling_one;
    kh[i] = kh_one;
    kb[i] = kb_one;
    thetab[i] = DEG2RAD * thetab_one;
    theta0[i] = DEG2RAD * 180.0;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Invalid angle type {}", arg[0]);
}

/* ---------------------------------------------------------------------- */

void AngleMesoCNT::init_style()
{
  std::string id_fix = "angle_mesocnt_buckled";
  if (!modify->get_fix_by_id(id_fix))
    modify->add_fix(id_fix + " all property/atom i_buckled ghost yes");
}

/* ---------------------------------------------------------------------- */

double AngleMesoCNT::equilibrium_angle(int /*i*/)
{
  return 180.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleMesoCNT::write_restart(FILE *fp)
{
  fwrite(&buckling[1], sizeof(int), atom->nangletypes, fp);
  fwrite(&kh[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&kb[1], sizeof(double), atom->nangletypes, fp);
  fwrite(&thetab[1], sizeof(double), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleMesoCNT::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &buckling[1], sizeof(int), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &kh[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &kb[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &thetab[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
  MPI_Bcast(&buckling[1], atom->nangletypes, MPI_INT, 0, world);
  MPI_Bcast(&kh[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&kb[1], atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&thetab[1], atom->nangletypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) {
    theta0[i] = 180.0;
    setflag[i] = 1;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleMesoCNT::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %d %g %g %g\n", i, buckling[i], kh[i], kb[i], RAD2DEG * thetab[i]);
}

/* ---------------------------------------------------------------------- */

double AngleMesoCNT::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

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

  double dtheta = acos(c) - theta0[type];

  // harmonic bending
  if (!buckling[type] || dtheta < thetab[type]) {
    double tk = kh[type] * dtheta;
    return tk * dtheta;
  }
  // bending buckling
  else
    return kb[type] * dtheta + thetab[type] * (kh[type] * thetab[type] - kb[type]);
}
