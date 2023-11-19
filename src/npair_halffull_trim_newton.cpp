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

#include "npair_halffull_trim_newton.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalffullTrimNewton::NPairHalffullTrimNewton(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

void NPairHalffullTrimNewton::build(NeighList *list)
{
  int i, j, ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq;

  const double delta = 0.01 * force->angstrom;
  const int triclinic = domain->triclinic;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int *ilist_full = list->listfull->ilist;
  int *numneigh_full = list->listfull->numneigh;
  int **firstneigh_full = list->listfull->firstneigh;
  int inum_full = list->listfull->inum;

  int inum = 0;
  ipage->reset();

  double cutsq_custom = cutoff_custom * cutoff_custom;

  // loop over parent full list

  for (ii = 0; ii < inum_full; ii++) {
    n = 0;
    neighptr = ipage->vget();

    i = ilist_full[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over full neighbor list
    // use i < j < nlocal to eliminate half the local/local interactions
    // for triclinic, must use delta to eliminate half the local/ghost interactions
    // cannot use I/J exact coord comparision as for orthog
    //   b/c transforming orthog -> lambda -> orthog for ghost atoms
    //   with an added PBC offset can shift all 3 coords by epsilon

    jlist = firstneigh_full[i];
    jnum = numneigh_full[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;

      if (j < nlocal) {
        if (i > j) continue;
      } else if (triclinic) {
        if (fabs(x[j][2]-ztmp) > delta) {
          if (x[j][2] < ztmp) continue;
        } else if (fabs(x[j][1]-ytmp) > delta) {
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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq > cutsq_custom) continue;

      neighptr[n++] = joriginal;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  list->inum = inum;
}
