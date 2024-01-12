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

#include "npair_skip_size_off2on.h"

#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int TRIM>
NPairSkipSizeOff2onTemp<TRIM>::NPairSkipSizeOff2onTemp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   parent non-skip list used newton off, this skip list is newton on
------------------------------------------------------------------------- */

template<int TRIM>
void NPairSkipSizeOff2onTemp<TRIM>::build(NeighList *list)
{
  int i, j, ii, jj, n, itype, jnum, joriginal;
  tagint itag, jtag;
  int *neighptr, *jlist;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int *ilist_skip = list->listskip->ilist;
  int *numneigh_skip = list->listskip->numneigh;
  int **firstneigh_skip = list->listskip->firstneigh;
  int inum_skip = list->listskip->inum;

  int *iskip = list->iskip;
  int **ijskip = list->ijskip;

  int inum = 0;
  ipage->reset();

  double **x = atom->x;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq;
  double cutsq_custom = cutoff_custom * cutoff_custom;

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;
    itag = tag[i];

    if (TRIM) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
    }

    n = 0;
    neighptr = ipage->vget();

    // loop over parent non-skip size list and optionally its history info

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;

      // only keep I,J when J = ghost if Itag < Jtag

      jtag = tag[j];
      if (j >= nlocal && jtag < itag) continue;

      if (TRIM) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
        if (rsq > cutsq_custom) continue;
      }

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

namespace LAMMPS_NS {
template class NPairSkipSizeOff2onTemp<0>;
template class NPairSkipSizeOff2onTemp<1>;
}
