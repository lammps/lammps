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

#include "compute_coord_atom_omp.h"

#include "atom.h"
#include "comm.h"
#include "compute_orientorder_atom.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCoordAtomOMP::ComputeCoordAtomOMP(LAMMPS *lmp, int narg, char **arg) :
    ComputeCoordAtom(lmp, narg, arg)
{
}

/* ----------------------------------------------------------------------
   threaded variant of ComputeCoordAtom::compute_peratom().  The serial
   setup (array (re)allocation, dependent orientorder/atom invocation,
   forward communication and neighbor-list build) is unchanged; only the
   per-atom loops are threaded.  Each atom writes solely to its own output
   slot (cvec[i] / carray[i]), so the result is bit-identical to the serial
   compute regardless of thread count or schedule.
------------------------------------------------------------------------- */

void ComputeCoordAtomOMP::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nmax > nmax) {
    if (ncol == 1) {
      memory->destroy(cvec);
      nmax = atom->nmax;
      memory->create(cvec, nmax, "coord/atom:cvec");
      vector_atom = cvec;
    } else {
      memory->destroy(carray);
      nmax = atom->nmax;
      memory->create(carray, nmax, ncol, "coord/atom:carray");
      array_atom = carray;
    }
  }

  if (cstyle == ORIENT) {
    if (!(c_orientorder->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_orientorder->compute_peratom();
      c_orientorder->invoked_flag |= Compute::INVOKED_PERATOM;
    }
    nqlist = c_orientorder->nqlist;
    normv = c_orientorder->array_atom;
    comm->forward_comm(this);
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // compute coordination number(s) for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const double *const *const x = atom->x;
  const int *const type = atom->type;
  const int *const mask = atom->mask;

  if (cstyle == CUTOFF) {

    if (ncol == 1) {

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
      for (int ii = 0; ii < inum; ii++) {
        const int i = ilist[ii];
        if (mask[i] & groupbit) {
          const double xtmp = x[i][0];
          const double ytmp = x[i][1];
          const double ztmp = x[i][2];
          const int *const jlist = firstneigh[i];
          const int jnum = numneigh[i];

          int n = 0;
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj] & NEIGHMASK;

            if (mask[j] & jgroupbit) {
              const int jtype = type[j];
              const double delx = xtmp - x[j][0];
              const double dely = ytmp - x[j][1];
              const double delz = ztmp - x[j][2];
              const double rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < cutsq && jtype >= typelo[0] && jtype <= typehi[0]) n++;
            }
          }

          cvec[i] = n;
        } else
          cvec[i] = 0.0;
      }

    } else {

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
      for (int ii = 0; ii < inum; ii++) {
        const int i = ilist[ii];
        double *const count = carray[i];
        for (int m = 0; m < ncol; m++) count[m] = 0.0;

        if (mask[i] & groupbit) {
          const double xtmp = x[i][0];
          const double ytmp = x[i][1];
          const double ztmp = x[i][2];
          const int *const jlist = firstneigh[i];
          const int jnum = numneigh[i];

          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj] & NEIGHMASK;

            if (mask[j] & jgroupbit) {
              const int jtype = type[j];
              const double delx = xtmp - x[j][0];
              const double dely = ytmp - x[j][1];
              const double delz = ztmp - x[j][2];
              const double rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < cutsq) {
                for (int m = 0; m < ncol; m++)
                  if (jtype >= typelo[m] && jtype <= typehi[m]) count[m] += 1.0;
              }
            }
          }
        }
      }
    }

  } else if (cstyle == ORIENT) {

    const int ncomp = 2 * (2 * l + 1);

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) {
        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const int *const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        int n = 0;
        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq < cutsq) {
            double dot_product = 0.0;
            for (int m = 0; m < ncomp; m++) {
              dot_product += normv[i][nqlist + m] * normv[j][nqlist + m];
            }
            if (dot_product > threshold) n++;
          }
        }
        cvec[i] = n;
      } else
        cvec[i] = 0.0;
    }
  }
}
