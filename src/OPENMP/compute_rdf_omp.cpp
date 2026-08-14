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

#include "compute_rdf_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeRDFOMP::ComputeRDFOMP(LAMMPS *lmp, int narg, char **arg) : ComputeRDF(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   threaded variant of ComputeRDF::compute_array().  Only the histogram
   tally loop is threaded: each thread accumulates into a private partial
   histogram, then the partials are summed into hist[] before the
   MPI_Allreduce.  Every tally is an increment of 1.0, so the per-thread
   sums are exact integers and the reduced result is bit-identical to the
   serial compute regardless of thread count.  Setup, normalization and the
   collective MPI_Allreduce remain serial/collective.
------------------------------------------------------------------------- */

void ComputeRDFOMP::compute_array()
{
  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // tally the RDF
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // consider I,J as one interaction even if neighbor pair is stored on 2 procs
  // tally I,J pair each time I is central atom, and each time J is central

  const double *const *const x = atom->x;
  const int *const type = atom->type;
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;

  const double *const special_coul = force->special_coul;
  const double *const special_lj = force->special_lj;
  const int newton_pair = force->newton_pair;

  // per-thread partial histograms; tallies are exact integer increments,
  // so summing the partials reproduces the serial histogram bit-for-bit

  const int nthreads = comm->nthreads;
  const int ntot = npairs * nbin;
  auto *const hist_thr = new double[(std::size_t) nthreads * ntot];

#if defined(_OPENMP)
#pragma omp parallel
#endif
  {
    int tid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    double *const myhist = hist_thr + (std::size_t) tid * ntot;
    for (int k = 0; k < ntot; k++) myhist[k] = 0.0;

#if defined(_OPENMP)
#pragma omp for schedule(dynamic)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        const double factor_lj = special_lj[sbmask(j)];
        const double factor_coul = special_coul[sbmask(j)];
        j &= NEIGHMASK;

        // if both weighting factors are 0, skip this pair
        // could be 0 and still be in neigh list for long-range Coulombics
        // want consistency with non-charged pairs which wouldn't be in list

        if (factor_lj == 0.0 && factor_coul == 0.0) continue;

        if (!(mask[j] & groupbit)) continue;
        const int jtype = type[j];
        const int ipair = nrdfpair[itype][jtype];
        const int jpair = nrdfpair[jtype][itype];
        if (!ipair && !jpair) continue;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double r = sqrt(delx * delx + dely * dely + delz * delz);
        const int ibin = static_cast<int>(r * delrinv);
        if (ibin >= nbin) continue;

        for (int ihisto = 0; ihisto < ipair; ihisto++) {
          const int m = rdfpair[ihisto][itype][jtype];
          myhist[m * nbin + ibin] += 1.0;
        }
        if (newton_pair || j < nlocal) {
          for (int ihisto = 0; ihisto < jpair; ihisto++) {
            const int m = rdfpair[ihisto][jtype][itype];
            myhist[m * nbin + ibin] += 1.0;
          }
        }
      }
    }
  }

  // reduce per-thread partial histograms into hist[] (hist[0] is contiguous)

  double *const hist_flat = hist[0];
  for (int k = 0; k < ntot; k++) {
    double sum = 0.0;
    for (int t = 0; t < nthreads; t++) sum += hist_thr[(std::size_t) t * ntot + k];
    hist_flat[k] = sum;
  }
  delete[] hist_thr;

  // sum histograms across procs

  MPI_Allreduce(hist[0], histall[0], npairs * nbin, MPI_DOUBLE, MPI_SUM, world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m
  // npairs = number of pairs, corrected for duplicates
  // duplicates = pairs in which both atoms are the same

  double constant, vfrac, gr, ncoord, rlower, rupper, normfac;

  if (domain->dimension == 3) {
    constant = 4.0 * MY_PI / (3.0 * domain->xprd * domain->yprd * domain->zprd);

    for (int m = 0; m < npairs; m++) {
      if (rev[m]) {    // swap i and j because they were entered in different order
        normfac = (jcount[m] > 0)
            ? static_cast<double>(icount[m]) - static_cast<double>(duplicates[m]) / jcount[m]
            : 0.0;
      } else {
        normfac = (icount[m] > 0)
            ? static_cast<double>(jcount[m]) - static_cast<double>(duplicates[m]) / icount[m]
            : 0.0;
      }
      ncoord = 0.0;
      for (int ibin = 0; ibin < nbin; ibin++) {
        rlower = ibin * delr;
        rupper = (ibin + 1) * delr;
        vfrac = constant * (rupper * rupper * rupper - rlower * rlower * rlower);
        gr = 0.0;
        if (vfrac * normfac != 0.0) {
          if (rev[m]) {    // swap i and j because they were entered in different order
            gr = histall[m][ibin] / (vfrac * normfac * jcount[m]);
            if (jcount[m] != 0) ncoord += gr * vfrac * normfac;
          } else {
            gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
            if (icount[m] != 0) ncoord += gr * vfrac * normfac;
          }
        }
        array[ibin][1 + 2 * m] = gr;
        array[ibin][2 + 2 * m] = ncoord;
      }
    }

  } else {
    constant = MY_PI / (domain->xprd * domain->yprd);

    for (int m = 0; m < npairs; m++) {
      ncoord = 0.0;
      normfac = (icount[m] > 0)
          ? static_cast<double>(jcount[m]) - static_cast<double>(duplicates[m]) / icount[m]
          : 0.0;
      for (int ibin = 0; ibin < nbin; ibin++) {
        rlower = ibin * delr;
        rupper = (ibin + 1) * delr;
        vfrac = constant * (rupper * rupper - rlower * rlower);
        if (vfrac * normfac != 0.0)
          gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
        else
          gr = 0.0;
        if (icount[m] != 0) ncoord += gr * vfrac * normfac;
        array[ibin][1 + 2 * m] = gr;
        array[ibin][2 + 2 * m] = ncoord;
      }
    }
  }
}
