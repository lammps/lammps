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

#include "compute_ke_atom_omp.h"

#include "atom.h"
#include "force.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKEAtomOMP::ComputeKEAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeKEAtom(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void ComputeKEAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ke);
    nmax = atom->nmax;
    memory->create(ke, nmax, "ke/atom:ke");
    vector_atom = ke;
  }

  // compute kinetic energy for each atom in group

  const double mvv2e = force->mvv2e;
  const double *const *const v = atom->v;
  const double *const mass = atom->mass;
  const double *const rmass = atom->rmass;
  const int *const mask = atom->mask;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;

  if (rmass) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        ke[i] =
            0.5 * mvv2e * rmass[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
      } else
        ke[i] = 0.0;
    }

  } else {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        ke[i] = 0.5 * mvv2e * mass[type[i]] *
            (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
      } else
        ke[i] = 0.0;
    }
  }
}
