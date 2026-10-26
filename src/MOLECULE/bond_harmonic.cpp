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

#include "bond_harmonic.h"

#include "atom.h"
#include "bond_degradation.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHarmonic::BondHarmonic(LAMMPS *_lmp) : Bond(_lmp), k(nullptr), r0(nullptr)
{
  born_matrix_enable = 1;
}
/* ---------------------------------------------------------------------- */

BondHarmonic::~BondHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonic::compute(int eflag, int vflag)
{
  int i1, i2, n, type;
  double delx, dely, delz, ebond, fbond;
  double rsq, r, dr, rk;

  ebond = 0.0;
  ev_init(eflag, vflag);
  const auto degradation_sources = get_degradation_sources();

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
    dr = r - r0[type];
    rk = degradation_scale(i1, i2, type, degradation_sources) * k[type] * dr;

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

void BondHarmonic::allocate()
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(k, np1, "bond:k");
  memory->create(r0, np1, "bond:r0");

  memory->create(setflag, np1, "bond:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonic::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  double k_one = utils::numeric(FLERR, arg[1], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHarmonic::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1], sizeof(double), atom->nbondtypes, fp);
  fwrite(&r0[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR, &k[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
    utils::sfread(FLERR, &r0[1], sizeof(double), atom->nbondtypes, fp, nullptr, error);
  }
  MPI_Bcast(&k[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(&r0[1], atom->nbondtypes, MPI_DOUBLE, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) fprintf(fp, "%d %g %g\n", i, k[i], r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondHarmonic::single(int type, double rsq, int i, int j, double &fforce)
{
  const auto degradation_sources = get_degradation_sources();
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = degradation_scale(i, j, type, degradation_sources) * k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0 * rk / r;
  return rk * dr;
}

/* ---------------------------------------------------------------------- */

void BondHarmonic::born_matrix(int type, double rsq, int i, int j, double &du, double &du2)
{
  const auto degradation_sources = get_degradation_sources();
  double r = sqrt(rsq);
  double dr = r - r0[type];
  du2 = 0.0;
  du = 0.0;
  du2 = 2 * degradation_scale(i, j, type, degradation_sources) * k[type];
  if (r > 0.0) du = du2 * dr;
}

/* ----------------------------------------------------------------------
   discover providers at each entry point so adding/replacing/removing a
   fix cannot leave stale pointers, including between run commands
------------------------------------------------------------------------- */

std::vector<const BondDegradation *> BondHarmonic::get_degradation_sources() const
{
  std::vector<const BondDegradation *> sources;
  for (int n = 0; n < modify->nfix; n++) {
    int dim = 0;
    auto *source = static_cast<BondDegradation *>(
        modify->fix[n]->extract("bond_degradation", dim));
    if (source) sources.push_back(source);
  }
  return sources;
}

/* ---------------------------------------------------------------------- */

double BondHarmonic::degradation_scale(
    int i, int j, int type, const std::vector<const BondDegradation *> &sources) const
{
  // bond_write can pass atom types in place of atom indices. Its probes
  // describe the original type potential, independently of damaged bonds.
  if (single_type_only || sources.empty() || i == j || i < 0 || j < 0 ||
      i >= atom->nlocal + atom->nghost || j >= atom->nlocal + atom->nghost)
    return 1.0;

  double scale = 1.0;
  for (const auto *source : sources)
    scale = std::min(scale, source->bond_scale(atom->tag[i], atom->tag[j], type));
  return scale;
}

/* ----------------------------------------------------------------------
   return ptr to internal members upon request
------------------------------------------------------------------------ */

void *BondHarmonic::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str, "k") == 0) return (void *) k;
  if (strcmp(str, "r0") == 0) return (void *) r0;
  return nullptr;
}
