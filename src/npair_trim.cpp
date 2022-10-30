/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Trimright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_trim.h"
#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairTrim::NPairTrim(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   create list which is a trimmed version of parent list
------------------------------------------------------------------------- */

void NPairTrim::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  double cutsq_custom = cutoff_custom * cutoff_custom;

  int ii, jj, n, jnum, joriginal;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq;

  double **x = atom->x;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  ipage->reset();

  int *ilist_copy = listcopy->ilist;
  int *numneigh_copy = listcopy->numneigh;
  int **firstneigh_copy = listcopy->firstneigh;
  int inum = listcopy->inum;

  list->inum = inum;
  list->gnum = listcopy->gnum;

  for (ii = 0; ii < inum; ii++) {
    n = 0;
    neighptr = ipage->vget();

    const int i = ilist_copy[ii];
    ilist[i] = i;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over copy list with larger cutoff and trim to shorter cutoff

    jlist = firstneigh_copy[i];
    jnum = numneigh_copy[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      const int j = joriginal & NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq > cutsq_custom) continue;

      neighptr[n++] = joriginal;
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
}
