/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include "npair_skip_size_off2on_oneside.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "fix_shear_history.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairSkipSizeOff2onOneside::NPairSkipSizeOff2onOneside(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   parent non-skip list used newton off and was not onesided,
     this skip list is newton on and onesided
   if list requests it, preserve shear history via fix shear/history
------------------------------------------------------------------------- */

void NPairSkipSizeOff2onOneside::build(NeighList *list)
{
  int i,j,ii,jj,m,n,nn,itype,jnum,joriginal,flip,dnum,dnumbytes,tmp;
  tagint jtag;
  int *surf,*jlist;

  int *npartner;
  tagint **partner;
  double **shearpartner;
  int **firsttouch;
  double **firstshear;
  MyPage<int> *ipage_touch;
  MyPage<double> *dpage_shear;
  NeighList *listhistory;

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

  if (domain->dimension == 2) surf = atom->line;
  else surf = atom->tri;

  FixShearHistory *fix_history = (FixShearHistory *) list->fix_history;
  if (fix_history) {
    fix_history->nlocal_neigh = nlocal;
    fix_history->nall_neigh = nlocal + atom->nghost;
    npartner = fix_history->npartner;
    partner = fix_history->partner;
    shearpartner = fix_history->shearpartner;
    listhistory = list->listhistory;
    firsttouch = listhistory->firstneigh;
    firstshear = listhistory->firstdouble;
    ipage_touch = listhistory->ipage;
    dpage_shear = listhistory->dpage;
    dnum = listhistory->dnum;
    dnumbytes = dnum * sizeof(double);
  }

  int inum = 0;
  ipage->reset();
  if (fix_history) {
    ipage_touch->reset();
    dpage_shear->reset();
  }

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

    n = 0;

    // loop over parent non-skip size list

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;

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

  // allocate all per-atom neigh list chunks, including history

  for (i = 0; i < nlocal; i++) {
    if (numneigh[i] == 0) continue;
    n = numneigh[i];
    firstneigh[i] = ipage->get(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    if (fix_history) {
      firsttouch[i] = ipage_touch->get(n);
      firstshear[i] = dpage_shear->get(dnum*n);
    }
  }

  // second loop over atoms in other list to store neighbors
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (i = 0; i < nlocal; i++) numneigh[i] = 0;

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;

    // loop over parent non-skip size list and optionally its history info

    jlist = firstneigh_skip[i];
    jnum = numneigh_skip[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      if (ijskip[itype][type[j]]) continue;

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

      // no numeric test for current touch
      // just use FSH partner list to infer it
      // would require complex calculation for surfs

      if (fix_history) {
        jtag = tag[j];
        n = numneigh[i];
        nn = dnum*n;
        for (m = 0; m < npartner[i]; m++)
          if (partner[i][m] == jtag) break;
        if (m < npartner[i]) {
          firsttouch[i][n] = 1;
          memcpy(&firstshear[i][nn],&shearpartner[i][dnum*m],dnumbytes);
        } else {
          firsttouch[i][n] = 0;
          memcpy(&firstshear[i][nn],zeroes,dnumbytes);
        }
      }

      numneigh[i]++;
      if (flip) i = j;
    }

    // only add atom I to ilist if it has neighbors
    // fix shear/history allows for this in pre_exchange_onesided()

    if (numneigh[i]) ilist[inum++] = i;
  }

  list->inum = inum;
}
