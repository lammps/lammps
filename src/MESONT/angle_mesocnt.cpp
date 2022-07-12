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

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::RAD2DEG;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleMesoCNT::AngleMesoCNT(LAMMPS *_lmp) : Angle(_lmp)
{
  kh = nullptr;
  kb = nullptr;
  thetab = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleMesoCNT::~AngleMesoCNT()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
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
    if (fabs(dtheta) < thetab[type]) {
      tk = kh[type] * dtheta;
      if (eflag) eangle = tk * dtheta;
      a = -2.0 * tk * s;

      buckled[i2] = 0;
    }
    // bending buckling
    else {
      if (eflag) eangle = kb[type] * fabs(dtheta) + thetab[type] * (kh[type] * thetab[type] - kb[type]);
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
  if (narg != 4) error->all(FLERR, "Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  double kh_one = utils::numeric(FLERR, arg[1], false, lmp);
  double kb_one = utils::numeric(FLERR, arg[2], false, lmp);
  double thetab_one = utils::numeric(FLERR, arg[3], false, lmp);

  // convert thetab from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kh[i] = kh_one;
    kb[i] = kb_one;
    thetab[i] = DEG2RAD * thetab_one;
    theta0[i] = DEG2RAD * 180.0;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for angle coefficients");
}

void AngleMesoCNT::init_style()
{
  char *id_fix = utils::strdup("angle_mesocnt_buckled");
  if (modify->find_fix(id_fix) < 0)
    modify->add_fix(std::string(id_fix)+" all property/atom i_buckled ghost yes");
}

/* ---------------------------------------------------------------------- */

double AngleMesoCNT::equilibrium_angle(int i)
{
  return 180.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleMesoCNT::write_restart(FILE *fp)
{
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
    utils::sfread(FLERR, &kh[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &kb[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
    utils::sfread(FLERR, &thetab[1], sizeof(double), atom->nangletypes, fp, nullptr, error);
  }
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
    fprintf(fp, "%d %g %g %g\n", i, kh[i], kb[i], RAD2DEG * thetab[i]);
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
  if (dtheta < thetab[type]) {
    double tk = kh[type] * dtheta;
    return tk * dtheta;
  }
  // bending buckling
  else return kb[type] * dtheta + thetab[type] * (kh[type] * thetab[type] - kb[type]);
}
