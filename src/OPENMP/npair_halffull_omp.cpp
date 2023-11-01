// clang-format off
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

#include "npair_halffull_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int NEWTON, int TRI, int TRIM>
NPairHalffullOmp<NEWTON, TRI, TRIM>::NPairHalffullOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   works if full list is a skip list

   Newtoff:
     pair stored by me if j is ghost (also stored by proc owning j)
     works for owned (non-ghost) list, also for ghost list
     if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
   Newton:
     if j is ghost, only store if j coords are "above and to the right" of i
     use i < j < nlocal to eliminate half the local/local interactions
   Newton + Triclinic:
     must use delta to eliminate half the local/ghost interactions
     cannot use I/J exact coord comparision as for orthog
       b/c transforming orthog -> lambda -> orthog for ghost atoms
       with an added PBC offset can shift all 3 coords by epsilon
------------------------------------------------------------------------- */

template<int NEWTON, int TRI, int TRIM>
void NPairHalffullOmp<NEWTON, TRI, TRIM>::build(NeighList *list)
{
  const int inum_full = list->listfull->inum;
  const double delta = 0.01 * force->angstrom;

  NPAIR_OMP_INIT;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(inum_full);

  int i, j, ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  double cutsq_custom = cutoff_custom * cutoff_custom;

  // loop over atoms in full list

  for (ii = ifrom; ii < ito; ii++) {

    n = 0;
    neighptr = ipage.vget();

    // loop over parent full list

    i = ilist_full[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;

      if (NEWTON) {
        if (j < nlocal) {
          if (i > j) continue;
        } else if (TRI) {
          if (fabs(x[j][2] - ztmp) > delta) {
            if (x[j][2] < ztmp) continue;
          } else if (fabs(x[j][1] - ytmp) > delta) {
            if (x[j][1] < ytmp) continue;
          } else {
            if (x[j][0] < xtmp) continue;
          }
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }

        if (TRIM) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq > cutsq_custom) continue;
        }
        neighptr[n++] = joriginal;
      } else {
        if (j > i) {
          if (TRIM) {
            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx * delx + dely * dely + delz * delz;

            if (rsq > cutsq_custom) continue;
          }
          neighptr[n++] = joriginal;
        }
      }
    }

    ilist[ii] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  NPAIR_OMP_CLOSE;
  list->inum = inum_full;
}

namespace LAMMPS_NS {
template class NPairHalffullOmp<0,0,0>;
template class NPairHalffullOmp<1,0,0>;
template class NPairHalffullOmp<1,1,0>;
template class NPairHalffullOmp<0,0,1>;
template class NPairHalffullOmp<1,0,1>;
template class NPairHalffullOmp<1,1,1>;
}
