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

#include "npair_skip_trim_size_off2on_oneside.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairSkipTrimSizeOff2onOneside::NPairSkipTrimSizeOff2onOneside(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   parent non-skip list used newton off and was not onesided,
     this skip list is newton on and onesided
------------------------------------------------------------------------- */

void NPairSkipTrimSizeOff2onOneside::build(NeighList *list)
{
  int i,j,ii,jj,itype,jnum,joriginal,flip,tmp;
  int *surf,*jlist;

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

  if (domain->dimension == 2) surf = atom->line;
  else surf = atom->tri;

  int inum = 0;
  ipage->reset();

  double **x = atom->x;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz, rsq;
  double cutsq_custom = cutoff_custom * cutoff_custom;

  // two loops over parent list required, one to count, one to store
  // because onesided constraint means pair I,J may be stored with I or J
  // so don't know in advance how much space to alloc for each atom's neighs

  // first loop over atoms in other list to count neighbors
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over parent non-skip size list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq_custom) continue;

      // flip I,J if necessary to satisfy onesided constraint
      // do not keep if I is now ghost

      if (surf[i] >= 0) {
        if (j >= nlocal) continue;
        tmp = i;
        i = j;
        j = tmp;
        flip = 1;
      } else flip = 0;

      numneigh[i]++;
      if (flip) i = j;
    }
  }

  // allocate all per-atom neigh list chunks

  for (i = 0; i < nlocal; i++) {
    if (numneigh[i] == 0) continue;
    firstneigh[i] = ipage->get(numneigh[i]);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  // second loop over atoms in other list to store neighbors
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over parent non-skip size list and optionally its history info

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq > cutsq_custom) continue;

      // flip I,J if necessary to satisfy onesided constraint
      // do not keep if I is now ghost

      if (surf[i] >= 0) {
        if (j >= nlocal) continue;
        tmp = i;
        i = j;
        j = tmp;
        flip = 1;
      } else flip = 0;

      // store j in neigh list, not joriginal, like other neigh methods
      // OK, b/c there is no special list flagging for surfs

      firstneigh[i][numneigh[i]] = j;
      numneigh[i]++;
      if (flip) i = j;
    }

    // only add atom I to ilist if it has neighbors

    if (numneigh[i]) ilist[inum++] = i;
  }

  list->inum = inum;
}
