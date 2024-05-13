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

#include "compute_gyration.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGyration::ComputeGyration(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Illegal compute gyration command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeGyration::~ComputeGyration()
{
  delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeGyration::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

double ComputeGyration::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double xcm[3];
  if (group->dynamic[igroup]) masstotal = group->mass(igroup);
  group->xcm(igroup, masstotal, xcm);
  scalar = group->gyration(igroup, masstotal, xcm);
  return scalar;
}

/* ----------------------------------------------------------------------
   compute the radius-of-gyration tensor of group of atoms
   around center-of-mass cm
   must unwrap atoms to compute Rg tensor correctly
------------------------------------------------------------------------- */

void ComputeGyration::compute_vector()
{
  invoked_vector = update->ntimestep;

  double xcm[3];
  group->xcm(igroup, masstotal, xcm);

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double dx, dy, dz, massone;
  double unwrap[3];

  double rg[6];
  rg[0] = rg[1] = rg[2] = rg[3] = rg[4] = rg[5] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass)
        massone = rmass[i];
      else
        massone = mass[type[i]];

      domain->unmap(x[i], image[i], unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];

      rg[0] += dx * dx * massone;
      rg[1] += dy * dy * massone;
      rg[2] += dz * dz * massone;
      rg[3] += dx * dy * massone;
      rg[4] += dx * dz * massone;
      rg[5] += dy * dz * massone;
    }
  MPI_Allreduce(rg, vector, 6, MPI_DOUBLE, MPI_SUM, world);

  if (masstotal > 0.0)
    for (int i = 0; i < 6; i++) vector[i] /= masstotal;
}
