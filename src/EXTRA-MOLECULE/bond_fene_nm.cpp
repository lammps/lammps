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

#include "bond_fene_nm.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondFENENM::BondFENENM(LAMMPS *lmp) : BondFENE(lmp), nn(nullptr), mm(nullptr)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

BondFENENM::~BondFENENM()
{
  if (allocated && !copymode) {
    memory->destroy(nn);
    memory->destroy(mm);
  }
}

/* ---------------------------------------------------------------------- */

void BondFENENM::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r0sq, rlogarg;
  double r;

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

    // force from log term

    rsq = delx * delx + dely * dely + delz * delz;
    r0sq = r0[type] * r0[type];
    rlogarg = 1.0 - rsq / r0sq;

    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort
    // change cutuff from .1 to .02 so only bond lengths > 1.485 give the warning
    // and crash the run if rlogarg < -.21 rather than < -3
    // Don't print out warnings, only errors
    if (rlogarg < .02) {
      error->warning(FLERR, "fene/nm/split bond too long: {} {} {} {}", update->ntimestep,
                     atom->tag[i1], atom->tag[i2], sqrt(rsq));
      if (rlogarg <= -.21) error->one(FLERR, "Bad FENE bond");
      rlogarg = 0.02;
    }

    fbond = -k[type] / rlogarg;
    // force from n-m term
    if (rsq < sigma[type] * sigma[type]) {
      r = sqrt(rsq);
      fbond += epsilon[type] * (nn[type] * mm[type] / (nn[type] - mm[type])) *
          (pow(sigma[type] / r, nn[type]) - pow(sigma[type] / r, mm[type])) / rsq;
    }

    // energy

    if (eflag) {
      ebond = -0.5 * k[type] * r0sq * log(rlogarg);
      if (rsq < sigma[type] * sigma[type])
        ebond += (epsilon[type] / (nn[type] - mm[type])) *
            (mm[type] * pow(sigma[type] / r, nn[type]) - nn[type] * pow(sigma[type] / r, mm[type]));
    }
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

void BondFENENM::allocate()
{
  BondFENE::allocate();
  int n = atom->nbondtypes + 1;
  memory->create(nn, n, "bond:nn");
  memory->create(mm, n, "bond:mm");
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondFENENM::coeff(int narg, char **arg)
{
  if (narg != 7) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double epsilon_one = utils::numeric(FLERR, arg[3], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[4], false, lmp);
  double nn_one = utils::numeric(FLERR, arg[5], false, lmp);
  double mm_one = utils::numeric(FLERR, arg[6], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    epsilon[i] = epsilon_one;
    sigma[i] = sigma_one;
    nn[i] = nn_one;
    mm[i] = mm_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondFENENM::init_style()
{
  // special bonds should be 0 1 1

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 || force->special_lj[3] != 1.0) {
    if (comm->me == 0) error->warning(FLERR, "Use special bonds = 0,1,1 with bond style fene");
  }
}

/* ---------------------------------------------------------------------- */

double BondFENENM::equilibrium_distance(int i)
{
  return 0.97 * sigma[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondFENENM::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&epsilon[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&sigma[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&nn[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&mm[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondFENENM::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &epsilon[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &sigma[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &nn[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &mm[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&epsilon[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&sigma[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&nn[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mm[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondFENENM::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %g %g %g %g %g %g\n", i, k[i], r0[i], epsilon[i], sigma[i], nn[i], mm[i]);
}

/* ---------------------------------------------------------------------- */

double BondFENENM::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r0sq = r0[type] * r0[type];
  double rlogarg = 1.0 - rsq / r0sq;
  double r;
  // if r -> r0, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*r0 something serious is wrong, abort

  // change cutuff from .1 to .02 so only bond lengths > 1.485 give the warning
  // and crash the run if rlogarg < -.21 rather than < -3
  // Don't print out warnings, only errors
  if (rlogarg < 0.02) {
    error->warning(FLERR, "FENE bond too long: {} {:.8}", update->ntimestep, sqrt(rsq));
    if (rlogarg <= -.21) error->one(FLERR, "Bad FENE bond");
    rlogarg = 0.02;
  }

  double eng = -0.5 * k[type] * r0sq * log(rlogarg);
  fforce = -k[type] / rlogarg;

  if (rsq < sigma[type] * sigma[type]) {
    r = sqrt(rsq);
    fforce += epsilon[type] * (nn[type] * mm[type] / (nn[type] - mm[type])) *
        (pow(sigma[type] / r, nn[type]) - pow(sigma[type] / r, mm[type])) / rsq;
    eng += (epsilon[type] / (nn[type] - mm[type])) *
        (mm[type] * pow(sigma[type] / r, nn[type]) - nn[type] * pow(sigma[type] / r, mm[type]));
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void BondFENENM::born_matrix(int type, double rsq, int /*i*/, int /*j*/, double &du, double &du2)
{
  double r = sqrt(rsq);
  double r0sq = r0[type] * r0[type];
  double rlogarg = 1.0 - rsq / r0sq;

  // Contribution from the attractive term
  du = k[type] * r / rlogarg;
  du2 = k[type] * (1.0 + rsq / r0sq) / (rlogarg * rlogarg);

  // Contribution from the repulsive Lennard-Jones term
  if (rsq < sigma[type] * sigma[type]) {
    double prefactor = epsilon[type] * nn[type] * mm[type] / (nn[type] - mm[type]);
    du += prefactor * (pow(sigma[type] / r, mm[type]) - pow(sigma[type] / r, nn[type])) / r;
    du2 += prefactor * ((nn[type] + 1.0) * pow(sigma[type] / r, nn[type]) -
                    (mm[type] + 1.0) * pow(sigma[type] / r, mm[type])) / rsq;
  }
}

/* ---------------------------------------------------------------------- */

void *BondFENENM::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "k") == 0) return (void *) k;
  if (strcmp(str, "r0") == 0) return (void *) r0;
  return nullptr;
}
