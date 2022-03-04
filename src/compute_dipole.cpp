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

#include "compute_dipole.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum { MASSCENTER, GEOMCENTER };

/* ---------------------------------------------------------------------- */

ComputeDipole::ComputeDipole(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if ((narg < 3) || (narg > 4)) error->all(FLERR, "Illegal compute dipole command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  extscalar = 0;
  extvector = 0;

  vector = new double[size_vector];
  vector[0] = vector[1] = vector[2] = 0.0;
  usecenter = MASSCENTER;

  if (narg == 4) {
    if (utils::strmatch(arg[3], "^geom"))
      usecenter = GEOMCENTER;
    else if (strcmp(arg[3], "mass") == 0)
      usecenter = MASSCENTER;
    else
      error->all(FLERR, "Illegal compute dipole command");
  }
}

/* ---------------------------------------------------------------------- */

ComputeDipole::~ComputeDipole()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeDipole::compute_vector()
{
  invoked_vector = update->ntimestep;

  const auto x = atom->x;
  const auto mask = atom->mask;
  const auto type = atom->type;
  const auto image = atom->image;
  const auto mass = atom->mass;
  const auto rmass = atom->rmass;
  const auto q = atom->q;
  const auto mu = atom->mu;
  const auto nlocal = atom->nlocal;

  double dipole[3] = {0.0, 0.0, 0.0};
  double comproc[3] = {0.0, 0.0, 0.0};
  double com[3] = {0.0, 0.0, 0.0};
  double masstotal = 0.0;
  double chrgtotal = 0.0;
  double massproc = 0.0;
  double chrgproc = 0.0;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      double unwrap[3];
      double massone = 1.0;    // for usecenter == GEOMCENTER
      if (usecenter == MASSCENTER) {
        if (rmass)
          massone = rmass[i];
        else
          massone = mass[type[i]];
      }
      massproc += massone;
      if (atom->q_flag) chrgproc += q[i];
      domain->unmap(x[i], image[i], unwrap);
      comproc[0] += unwrap[0] * massone;
      comproc[1] += unwrap[1] * massone;
      comproc[2] += unwrap[2] * massone;
    }
  }
  MPI_Allreduce(&massproc, &masstotal, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&chrgproc, &chrgtotal, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(comproc, com, 3, MPI_DOUBLE, MPI_SUM, world);

  if (masstotal > 0.0) {
    com[0] /= masstotal;
    com[1] /= masstotal;
    com[2] /= masstotal;
  }

  // compute dipole moment

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double unwrap[3];
      domain->unmap(x[i], image[i], unwrap);
      if (atom->q_flag) {
        dipole[0] += q[i] * unwrap[0];
        dipole[1] += q[i] * unwrap[1];
        dipole[2] += q[i] * unwrap[2];
      }
      if (atom->mu_flag) {
        dipole[0] += mu[i][0];
        dipole[1] += mu[i][1];
        dipole[2] += mu[i][2];
      }
    }
  }

  MPI_Allreduce(dipole, vector, 3, MPI_DOUBLE, MPI_SUM, world);

  // correct for position dependence with a net charged group
  vector[0] -= chrgtotal * com[0];
  vector[1] -= chrgtotal * com[1];
  vector[2] -= chrgtotal * com[2];
}

/* ---------------------------------------------------------------------- */

double ComputeDipole::compute_scalar()
{
  if (invoked_vector != update->ntimestep) compute_vector();

  invoked_scalar = update->ntimestep;
  scalar = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  return scalar;
}
