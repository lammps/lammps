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

#include "compute_entropy_atom_omp.h"

#include "atom.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeEntropyAtomOMP::ComputeEntropyAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeEntropyAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeEntropyAtom::compute_peratom().  The per-atom
   g(r) scratch (gofr/integrand) is allocated thread-locally; the read-only
   rbin/rbinsq tables are shared.  When avg_flag is set the averaging pass
   reads neighbors' pair_entropy[j] (including ghost atoms, which is why the
   first loop spans inum+gnum), so it runs as a second omp-for region whose
   implicit barrier guarantees all pair_entropy[] are filled first.  Each
   atom writes only its own output, so the result is bit-identical to the
   serial compute.
------------------------------------------------------------------------- */

void ComputeEntropyAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // read-only distance-bin tables (shared across threads)

  auto *const rbin = new double[nbin];
  auto *const rbinsq = new double[nbin];
  for (int i = 0; i < nbin; i++) {
    rbin[i] = i * deltar;
    rbinsq[i] = rbin[i] * rbin[i];
  }

  // grow pair_entropy and pair_entropy_avg array if necessary

  if (atom->nmax > nmax) {
    if (!avg_flag) {
      memory->destroy(pair_entropy);
      nmax = atom->nmax;
      memory->create(pair_entropy, nmax, "entropy/atom:pair_entropy");
      vector_atom = pair_entropy;
    } else {
      memory->destroy(pair_entropy);
      memory->destroy(pair_entropy_avg);
      nmax = atom->nmax;
      memory->create(pair_entropy, nmax, "entropy/atom:pair_entropy");
      memory->create(pair_entropy_avg, nmax, "entropy/atom:pair_entropy_avg");
      vector_atom = pair_entropy_avg;
    }
  }

  // invoke occasional neighbor list build (if not perpetual)

  if (!avg_flag) neighbor->build_one(list);

  const int inum = list->inum + list->gnum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // Compute some constants

  const double sigmasq2 = 2 * sigma * sigma;
  const double volume = domain->xprd * domain->yprd * domain->zprd;
  const double global_density = atom->natoms / volume;

  // compute pair entropy for each atom in group
  // use full neighbor list

  const double *const *const x = atom->x;
  const int *const mask = atom->mask;

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    auto *const gofr = new double[nbin];
    auto *const integrand = new double[nbin];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        // If local density is used, calculate it

        double density = global_density;
        if (local_flag) {
          const double neigh_cutoff = force->pair->cutforce + neighbor->skin;
          const double vol = (4. / 3.) * MY_PI * neigh_cutoff * neigh_cutoff * neigh_cutoff;
          density = jnum / vol;
        }

        // calculate kernel normalization
        // Normalization of g(r)
        double normConstantBase = 4 * MY_PI * density;
        // Normalization of gaussian
        normConstantBase *= sqrt(2. * MY_PI) * sigma;
        const double invNormConstantBase = 1. / normConstantBase;

        // loop over list of all neighbors within force cutoff
        // initialize gofr

        for (int k = 0; k < nbin; ++k) gofr[k] = 0.;

        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            // contribute to gofr
            const double r = sqrt(rsq);
            const int bin = floor(r / deltar);    // NOLINT
            int minbin, maxbin;
            minbin = bin - deltabin;
            if (minbin < 0) minbin = 0;
            if (minbin > (nbin - 1)) minbin = nbin - 1;
            maxbin = bin + deltabin;
            if (maxbin > (nbin - 1)) maxbin = nbin - 1;
            for (int k = minbin; k < maxbin + 1; k++) {
              const double invNormKernel = invNormConstantBase / rbinsq[k];
              const double distance = r - rbin[k];
              gofr[k] += invNormKernel * exp(-distance * distance / sigmasq2);
            }
          }
        }

        // Calculate integrand

        for (int k = 0; k < nbin; ++k) {
          if (gofr[k] < 1.e-10) {
            integrand[k] = rbinsq[k];
          } else {
            integrand[k] = (gofr[k] * log(gofr[k]) - gofr[k] + 1) * rbinsq[k];
          }
        }

        // Integrate with trapezoid rule

        double value = 0.;
        for (int k = 1; k < nbin - 1; ++k) { value += integrand[k]; }
        value += 0.5 * integrand[0];
        value += 0.5 * integrand[nbin - 1];
        value *= deltar;

        pair_entropy[i] = -2 * MY_PI * density * value;
      }
    }

    delete[] gofr;
    delete[] integrand;
  }

  if (avg_flag) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        pair_entropy_avg[i] = pair_entropy[i];
        double counter = 1;

        // loop over list of all neighbors within force cutoff
        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq2) {
            pair_entropy_avg[i] += pair_entropy[j];
            counter += 1;
          }
        }
        pair_entropy_avg[i] /= counter;
      }
    }
  }

  delete[] rbin;
  delete[] rbinsq;
}
