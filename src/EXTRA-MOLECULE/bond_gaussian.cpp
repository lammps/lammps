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

#include "bond_gaussian.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 2.0e-308;

/* ---------------------------------------------------------------------- */

BondGaussian::BondGaussian(LAMMPS *lmp) :
    Bond(lmp), nterms(nullptr), bond_temperature(nullptr), alpha(nullptr), width(nullptr),
    r0(nullptr)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

BondGaussian::~BondGaussian()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(nterms);
    memory->destroy(bond_temperature);
    for (int i = 1; i <= atom->nbondtypes; i++) {
      delete[] alpha[i];
      delete[] width[i];
      delete[] r0[i];
    }
    delete[] alpha;
    delete[] width;
    delete[] r0;
  }
}

/* ---------------------------------------------------------------------- */

void BondGaussian::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, dr;
  double prefactor, exponent, g_i, sum_g_i, sum_numerator;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);

    sum_g_i = 0.0;
    sum_numerator = 0.0;
    for (int i = 0; i < nterms[type]; i++) {
      dr = r - r0[type][i];
      prefactor = (alpha[type][i] / (width[type][i] * sqrt(MY_PI2)));
      exponent = -2.0 * dr * dr / (width[type][i] * width[type][i]);
      g_i = prefactor * exp(exponent);
      sum_g_i += g_i;
      sum_numerator += g_i * dr / (width[type][i] * width[type][i]);
    }

    // avoid overflow
    if (sum_g_i < sum_numerator * SMALL) sum_g_i = sum_numerator * SMALL;

    // force & energy
    if (r > 0.0)
      fbond = -4.0 * (force->boltz * bond_temperature[type]) * (sum_numerator / sum_g_i) / r;
    else
      fbond = 0.0;

    if (eflag) ebond = -(force->boltz * bond_temperature[type]) * log(sum_g_i);

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx * fbond;
      f[i1][1] += dely * fbond;
      f[i1][2] += delz * fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx * fbond;
      f[i2][1] -= dely * fbond;
      f[i2][2] -= delz * fbond;
    }

    if (evflag) ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondGaussian::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes + 1;

  memory->create(nterms, n, "bond:nterms");
  memory->create(bond_temperature, n, "bond:bond_temperature");

  alpha = new double *[n];
  width = new double *[n];
  r0 = new double *[n];
  memset(alpha, 0, sizeof(double *) * n);
  memset(width, 0, sizeof(double *) * n);
  memset(r0, 0, sizeof(double *) * n);

  memory->create(setflag, n, "bond:setflag");
  memset(setflag, 0, sizeof(int) * n);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondGaussian::coeff(int narg, char **arg)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "bond_coeff", error);

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double bond_temp_one = utils::numeric(FLERR, arg[1], false, lmp);
  int n = utils::inumeric(FLERR, arg[2], false, lmp);
  if (n < 1) error->all(FLERR, "Invalid bond style gaussian value for n: {}", n);
  if (narg != 3 * n + 3) utils::missing_cmd_args(FLERR, "bond_coeff", error);

  if (!allocated) allocate();

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    bond_temperature[i] = bond_temp_one;
    nterms[i] = n;
    delete[] alpha[i];
    alpha[i] = new double[n];
    delete[] width[i];
    width[i] = new double[n];
    delete[] r0[i];
    r0[i] = new double[n];
    for (int j = 0; j < n; j++) {
      alpha[i][j] = utils::numeric(FLERR, arg[3 + 3 * j], false, lmp);
      if (alpha[i][j] <= 0.0) error->all(FLERR, "Invalid value for A_{}: {}", j, alpha[i][j]);
      width[i][j] = utils::numeric(FLERR, arg[4 + 3 * j], false, lmp);
      if (width[i][j] <= 0.0) error->all(FLERR, "Invalid value for w_{}: {}", j, width[i][j]);
      r0[i][j] = utils::numeric(FLERR, arg[5 + 3 * j], false, lmp);
      if (r0[i][j] <= 0.0) error->all(FLERR, "Invalid value for r0_{}: {}", j, r0[i][j]);
      setflag[i] = 1;
    }
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondGaussian::equilibrium_distance(int i)
{
  return r0[i][0];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondGaussian::write_restart(FILE *fp)
{
  fwrite(&bond_temperature[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&nterms[1], sizeof(int), atom->nbondtypes, fp);
  for (int i = 1; i <= atom->nbondtypes; i++) {
    fwrite(alpha[i], sizeof(double), nterms[i], fp);
    fwrite(width[i], sizeof(double), nterms[i], fp);
    fwrite(r0[i], sizeof(double), nterms[i], fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGaussian::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &bond_temperature[1], sizeof(double), atom->nbondtypes, fp, nullptr,
                  error);
    utils::sfread(FLERR, &nterms[1], sizeof(int), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&bond_temperature[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&nterms[1], atom->nbondtypes, MPI_INT, 0, world);

  // allocate
  for (int i = 1; i <= atom->nbondtypes; i++) {
    alpha[i] = new double[nterms[i]];
    width[i] = new double[nterms[i]];
    r0[i] = new double[nterms[i]];
  }

  if (comm->me == 0) {
    for (int i = 1; i <= atom->nbondtypes; i++) {
      utils::sfread(FLERR, alpha[i], sizeof(double), nterms[i], fp, nullptr, error);
      utils::sfread(FLERR, width[i], sizeof(double), nterms[i], fp, nullptr, error);
      utils::sfread(FLERR, r0[i], sizeof(double), nterms[i], fp, nullptr, error);
    }
  }

  for (int i = 1; i <= atom->nbondtypes; i++) {
    MPI_Bcast(alpha[i], nterms[i], MPI_DOUBLE, 0, world);
    MPI_Bcast(width[i], nterms[i], MPI_DOUBLE, 0, world);
    MPI_Bcast(r0[i], nterms[i], MPI_DOUBLE, 0, world);
  }

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondGaussian::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) {
    fprintf(fp, "%d %g %d", i, bond_temperature[i], nterms[i]);
    for (int j = 0; j < nterms[i]; j++) {
      fprintf(fp, " %g %g %g", alpha[i][j], width[i][j], r0[i][j]);
    }
    fprintf(fp, "\n");
  }
}

/* ---------------------------------------------------------------------- */

double BondGaussian::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r = sqrt(rsq);
  fforce = 0;

  double sum_g_i = 0.0;
  double sum_numerator = 0.0;
  for (int i = 0; i < nterms[type]; i++) {
    double dr = r - r0[type][i];
    double prefactor = (alpha[type][i] / (width[type][i] * sqrt(MY_PI2)));
    double exponent = -2 * dr * dr / (width[type][i] * width[type][i]);
    double g_i = prefactor * exp(exponent);
    sum_g_i += g_i;
    sum_numerator += g_i * dr / (width[type][i] * width[type][i]);
  }

  if (sum_g_i < SMALL) sum_g_i = SMALL;
  if (r > 0.0)
    fforce = -4.0 * (force->boltz * bond_temperature[type]) * (sum_numerator / sum_g_i) / r;

  return -(force->boltz * bond_temperature[type]) * log(sum_g_i);
}

/* ---------------------------------------------------------------------- */

void BondGaussian::born_matrix(int type, double rsq, int /*i*/, int /*j*/, double &du, double &du2)
{
  double r = sqrt(rsq);

  // first derivative of energy with respect to distance
  double sum_g_i = 0.0;
  double sum_numerator = 0.0;
  for (int i = 0; i < nterms[type]; i++) {
    double dr = r - r0[type][i];
    double prefactor = (alpha[type][i] / (width[type][i] * sqrt(MY_PI2)));
    double exponent = -2 * dr * dr / (width[type][i] * width[type][i]);
    double g_i = prefactor * exp(exponent);
    sum_g_i += g_i;
    sum_numerator += g_i * dr / (width[type][i] * width[type][i]);
  }

  if (sum_g_i < SMALL) sum_g_i = SMALL;
  du = 4.0 * (force->boltz * bond_temperature[type]) * (sum_numerator / sum_g_i);

  // second derivative of energy with respect to distance
  sum_g_i = 0.0;
  double sum_dg_i = 0.0;
  double sum_d2g_i = 0.0;
  for (int i = 0; i < nterms[type]; i++) {
    double dr = r - r0[type][i];
    double prefactor = (alpha[type][i] / (width[type][i] * sqrt(MY_PI2)));
    double exponent = -2 * dr * dr / (width[type][i] * width[type][i]);
    double g_i = prefactor * exp(exponent);
    sum_g_i += g_i;
    sum_dg_i -= 4.0 * g_i * dr / pow(width[type][i], 2);
    sum_d2g_i += 4.0 *  g_i * (4.0 * pow(r0[type][i], 2) - 8.0 * r0[type][i] * r - pow(width[type][i], 2) + 4.0 * r * r) / pow(width[type][i], 4) ;
  }

  if (sum_g_i < SMALL) sum_g_i = SMALL;
  double numerator = sum_d2g_i*sum_g_i - sum_dg_i*sum_dg_i;
  double denominator = sum_g_i * sum_g_i;

  du2 = - (force->boltz * bond_temperature[type]) * numerator / denominator;
}
