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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "dihedral_harmonic.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double TOLERANCE = 0.05;

/* ---------------------------------------------------------------------- */

DihedralHarmonic::DihedralHarmonic(LAMMPS *_lmp) : Dihedral(_lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralHarmonic::~DihedralHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(sign);
    memory->destroy(multiplicity);
    memory->destroy(cos_shift);
    memory->destroy(sin_shift);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHarmonic::compute(int eflag, int vflag)
{
  int i1, i2, i3, i4, i, m, n, type;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z, vb2xm, vb2ym, vb2zm;
  double edihedral, f1[3], f2[3], f3[3], f4[3];
  double ax, ay, az, bx, by, bz, rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv;
  double df, df1, ddf1, fg, hg, fga, hgb, gaa, gbb;
  double dtfx, dtfy, dtfz, dtgx, dtgy, dtgz, dthx, dthy, dthz;
  double c, s, p, sx2, sy2, sz2;

  edihedral = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    // c,s calculation

    ax = vb1y * vb2zm - vb1z * vb2ym;
    ay = vb1z * vb2xm - vb1x * vb2zm;
    az = vb1x * vb2ym - vb1y * vb2xm;
    bx = vb3y * vb2zm - vb3z * vb2ym;
    by = vb3z * vb2xm - vb3x * vb2zm;
    bz = vb3x * vb2ym - vb3y * vb2xm;

    rasq = ax * ax + ay * ay + az * az;
    rbsq = bx * bx + by * by + bz * bz;
    rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
    rg = sqrt(rgsq);

    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0 / rg;
    if (rasq > 0) ra2inv = 1.0 / rasq;
    if (rbsq > 0) rb2inv = 1.0 / rbsq;
    rabinv = sqrt(ra2inv * rb2inv);

    c = (ax * bx + ay * by + az * bz) * rabinv;
    s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) problem(FLERR, i1, i2, i3, i4);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    m = multiplicity[type];
    p = 1.0;
    ddf1 = df1 = 0.0;

    for (i = 0; i < m; i++) {
      ddf1 = p * c - df1 * s;
      df1 = p * s + df1 * c;
      p = ddf1;
    }

    p = p * cos_shift[type] + df1 * sin_shift[type];
    df1 = df1 * cos_shift[type] - ddf1 * sin_shift[type];
    df1 *= -m;
    p += 1.0;

    if (m == 0) {
      p = 1.0 + cos_shift[type];
      df1 = 0.0;
    }

    if (eflag) edihedral = k[type] * p;

    fg = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
    hg = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
    fga = fg * ra2inv * rginv;
    hgb = hg * rb2inv * rginv;
    gaa = -ra2inv * rg;
    gbb = rb2inv * rg;

    dtfx = gaa * ax;
    dtfy = gaa * ay;
    dtfz = gaa * az;
    dtgx = fga * ax - hgb * bx;
    dtgy = fga * ay - hgb * by;
    dtgz = fga * az - hgb * bz;
    dthx = gbb * bx;
    dthy = gbb * by;
    dthz = gbb * bz;

    df = -k[type] * df1;

    sx2 = df * dtgx;
    sy2 = df * dtgy;
    sz2 = df * dtgz;

    f1[0] = df * dtfx;
    f1[1] = df * dtfy;
    f1[2] = df * dtfz;

    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];

    f4[0] = df * dthx;
    f4[1] = df * dthy;
    f4[2] = df * dthz;

    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, edihedral, f1, f3, f4, vb1x, vb1y, vb1z, vb2x,
               vb2y, vb2z, vb3x, vb3y, vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHarmonic::allocate()
{
  allocated = 1;
  const int np1 = atom->ndihedraltypes + 1;

  memory->create(k, np1, "dihedral:k");
  memory->create(sign, np1, "dihedral:sign");
  memory->create(multiplicity, np1, "dihedral:multiplicity");
  memory->create(cos_shift, np1, "dihedral:cos_shift");
  memory->create(sin_shift, np1, "dihedral:sin_shift");

  memory->create(setflag, np1, "dihedral:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralHarmonic::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR, "Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  int sign_one = utils::inumeric(FLERR, arg[2], false, lmp);
  int multiplicity_one = utils::inumeric(FLERR, arg[3], false, lmp);

  // require sign = +/- 1 for backwards compatibility
  // arbitrary phase angle shift could be allowed, but would break
  //   backwards compatibility and is probably not needed

  if (sign_one != -1 && sign_one != 1)
    error->all(FLERR, "Incorrect sign arg for dihedral coefficients");
  if (multiplicity_one < 0)
    error->all(FLERR, "Incorrect multiplicity arg for dihedral coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    sign[i] = sign_one;
    if (sign[i] == 1) {
      cos_shift[i] = 1;
      sin_shift[i] = 0;
    } else {
      cos_shift[i] = -1;
      sin_shift[i] = 0;
    }
    multiplicity[i] = multiplicity_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->ndihedraltypes, fp);
  fwrite(&sign[1], sizeof(int), atom->ndihedraltypes, fp);
  fwrite(&multiplicity[1], sizeof(int), atom->ndihedraltypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->ndihedraltypes, fp, nullptr, error);
    utils::sfread(FLERR, &sign[1], sizeof(int), atom->ndihedraltypes, fp, nullptr, error);
    utils::sfread(FLERR, &multiplicity[1], sizeof(int), atom->ndihedraltypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sign[1], atom->ndihedraltypes, MPI_INT, 0, world);
  MPI_Bcast(&multiplicity[1], atom->ndihedraltypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
    if (sign[i] == 1) {
      cos_shift[i] = 1;
      sin_shift[i] = 0;
    } else {
      cos_shift[i] = -1;
      sin_shift[i] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp, "%d %g %d %d\n", i, k[i], sign[i], multiplicity[i]);
}
