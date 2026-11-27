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

#include "compute_erotate_sphere_atom_omp.h"

#include "atom.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeErotateSphereAtomOMP::ComputeErotateSphereAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeErotateSphereAtom(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void ComputeErotateSphereAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow erot array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(erot);
    nmax = atom->nmax;
    memory->create(erot, nmax, "erotate/sphere/atom:erot");
    vector_atom = erot;
  }

  // compute rotational kinetic energy for each atom in group
  // point particles will have erot = 0.0, due to radius = 0.0

  const double *const *const omega = atom->omega;
  const double *const radius = atom->radius;
  const double *const rmass = atom->rmass;
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      erot[i] =
          (omega[i][0] * omega[i][0] + omega[i][1] * omega[i][1] + omega[i][2] * omega[i][2]) *
          radius[i] * radius[i] * rmass[i];
      erot[i] *= pfactor;
    } else
      erot[i] = 0.0;
  }
}
