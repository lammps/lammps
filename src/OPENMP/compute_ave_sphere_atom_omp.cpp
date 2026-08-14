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

#include "compute_ave_sphere_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeAveSphereAtomOMP::ComputeAveSphereAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeAveSphereAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeAveSphereAtom::compute_peratom().  All per-atom
   scratch (p/vcom/vnet/count) is stack-local per iteration; ghost velocities
   are forward-communicated (serial) before the threaded loop.  Each atom
   writes only its own output row result[i], so the result is bit-identical
   to the serial compute.
------------------------------------------------------------------------- */

void ComputeAveSphereAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(result);
    nmax = atom->nmax;
    memory->create(result, nmax, 2, "ave/sphere/atom:result");
    array_atom = result;
  }

  // need velocities of ghost atoms

  comm->forward_comm(this);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const double *const *const x = atom->x;
  const double *const *const v = atom->v;
  const double *const mass = atom->mass;
  const double *const rmass = atom->rmass;
  const int *const type = atom->type;
  const int *const mask = atom->mask;

  const double adof = domain->dimension;
  const double mvv2e = force->mvv2e;
  const double mv2d = force->mv2d;
  const double boltz = force->boltz;

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];

    if (mask[i] & groupbit) {
      const double massone_i = rmass ? rmass[i] : mass[type[i]];

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      // i atom contribution

      int count = 1;
      double totalmass = massone_i;
      double p[3];
      p[0] = v[i][0] * massone_i;
      p[1] = v[i][1] * massone_i;
      p[2] = v[i][2] * massone_i;

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const double massone_j = rmass ? rmass[j] : mass[type[j]];

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          count++;
          totalmass += massone_j;
          p[0] += v[j][0] * massone_j;
          p[1] += v[j][1] * massone_j;
          p[2] += v[j][2] * massone_j;
        }
      }

      double vcom[3];
      vcom[0] = p[0] / totalmass;
      vcom[1] = p[1] / totalmass;
      vcom[2] = p[2] / totalmass;

      // i atom contribution

      double vnet[3];
      vnet[0] = v[i][0] - vcom[0];
      vnet[1] = v[i][1] - vcom[1];
      vnet[2] = v[i][2] - vcom[2];
      double ke_sum = massone_i * (vnet[0] * vnet[0] + vnet[1] * vnet[1] + vnet[2] * vnet[2]);

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const double massone_j = rmass ? rmass[j] : mass[type[j]];

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        if (rsq < cutsq) {
          vnet[0] = v[j][0] - vcom[0];
          vnet[1] = v[j][1] - vcom[1];
          vnet[2] = v[j][2] - vcom[2];
          ke_sum += massone_j * (vnet[0] * vnet[0] + vnet[1] * vnet[1] + vnet[2] * vnet[2]);
        }
      }
      const double density = mv2d * totalmass / volume;
      const double temp = mvv2e * ke_sum / (adof * count * boltz);
      result[i][0] = density;
      result[i][1] = temp;
    }
  }
}
