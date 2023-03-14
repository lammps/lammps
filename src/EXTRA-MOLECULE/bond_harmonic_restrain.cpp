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

#include "bond_harmonic_restrain.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHarmonicRestrain::BondHarmonicRestrain(LAMMPS *_lmp) : Bond(_lmp), initial(nullptr)
{
  writedata = 0;
  natoms = -1;
}

/* ---------------------------------------------------------------------- */

BondHarmonicRestrain::~BondHarmonicRestrain()
{
  if (initial) modify->delete_fix("BOND_RESTRAIN_X0");
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonicRestrain::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, r0, dr, rk;

  ebond = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **x0 = initial->astore;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x0[i1][0] - x0[i2][0];
    dely = x0[i1][1] - x0[i2][1];
    delz = x0[i1][2] - x0[i2][2];
    domain->minimum_image(delx, dely, delz);
    rsq = delx * delx + dely * dely + delz * delz;
    r0 = sqrt(rsq);

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx * delx + dely * dely + delz * delz;
    r = sqrt(rsq);
    dr = r - r0;
    rk = k[type] * dr;

    // force & energy

    if (r > 0.0)
      fbond = -2.0 * rk / r;
    else
      fbond = 0.0;

    if (eflag) ebond = rk * dr;

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

void BondHarmonicRestrain::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonicRestrain::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   initialize custom data storage
------------------------------------------------------------------------- */

void BondHarmonicRestrain::init_style()
{
  // store initial positions

  if (natoms < 0) {

    // create internal fix to store initial positions
    initial = dynamic_cast<FixStoreAtom *>(
      modify->add_fix("BOND_RESTRAIN_X0 all STORE/ATOM 3 0 1 1"));
    if (!initial) error->all(FLERR, "Failure to create internal per-atom storage");

    natoms = atom->natoms;
    double **x0 = initial->astore;
    const double *const *const x = atom->x;
    for (int i = 0; i < atom->nlocal; ++i)
      for (int j = 0; j < 3; ++j) x0[i][j] = x[i][j];

  } else {

    // after a restart, natoms is set but initial is a null pointer.
    // we add the fix, but do not initialize it.  It will pull the data from the restart.

    if (!initial) {
      initial = dynamic_cast<FixStoreAtom *>(
        modify->add_fix("BOND_RESTRAIN_X0 all STORE/ATOM 3 0 1 1"));
      if (!initial) error->all(FLERR, "Failure to create internal per-atom storage");
    }
  }

  // must not add  atoms
  if (natoms < atom->natoms)
    error->all(FLERR, "Bond style harmonic/restrain does not support adding atoms");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonicRestrain::write_restart(FILE *fp)
{
  fwrite(&natoms, sizeof(bigint), 1, fp);
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonicRestrain::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &natoms, sizeof(bigint), 1, fp, nullptr, error);
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&natoms, 1, MPI_LMP_BIGINT, 0, world);
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonicRestrain::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) fprintf(fp, "%d %g\n", i, k[i]);
}

/* ---------------------------------------------------------------------- */

double BondHarmonicRestrain::single(int type, double rsq, int i, int j, double &fforce)
{
  double **x0 = initial->astore;

  double delx = x0[i][0] - x0[j][0];
  double dely = x0[i][1] - x0[j][1];
  double delz = x0[i][2] - x0[j][2];
  domain->minimum_image(delx, dely, delz);
  double r0 = sqrt(delx * delx + dely * dely + delz * delz);

  double r = sqrt(rsq);
  double dr = r - r0;
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0 * rk / r;
  return rk * dr;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *BondHarmonicRestrain::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "k") == 0) return (void *) k;
  return nullptr;
}
