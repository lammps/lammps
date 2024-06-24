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

#include "compute_erotate_sphere.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr double INERTIA = 0.4;    // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

ComputeERotateSphere::ComputeERotateSphere(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal compute erotate/sphere command");

  scalar_flag = 1;
  extscalar = 1;

  // error check

  if (!atom->omega_flag) error->all(FLERR, "Compute erotate/sphere requires atom attribute omega");
}

/* ---------------------------------------------------------------------- */

void ComputeERotateSphere::init()
{
  pfactor = 0.5 * force->mvv2e * INERTIA;
}

/* ---------------------------------------------------------------------- */

double ComputeERotateSphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // point particles will not contribute, due to radius = 0.0

  double erotate = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      erotate +=
          (omega[i][0] * omega[i][0] + omega[i][1] * omega[i][1] + omega[i][2] * omega[i][2]) *
          radius[i] * radius[i] * rmass[i];

  MPI_Allreduce(&erotate, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= pfactor;
  return scalar;
}
