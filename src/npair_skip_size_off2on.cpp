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
#include "npair_skip_size_off2on.h"
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

NPairSkipSizeOff2on::NPairSkipSizeOff2on(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   iskip and ijskip flag which atom types and type pairs to skip
   parent non-skip list used newton off, this skip list is newton on
   if list requests it, preserve shear history via fix shear/history
------------------------------------------------------------------------- */

void NPairSkipSizeOff2on::build(NeighList *list)
{
  int i,j,ii,jj,m,n,nn,itype,jnum,joriginal,dnum,dnumbytes;
  tagint itag,jtag;
  int *neighptr,*jlist,*touchptr;
  double *shearptr;

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

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  for (ii = 0; ii < inum_skip; ii++) {
    i = ilist_skip[ii];
    itype = type[i];
    if (iskip[itype]) continue;
    itag = tag[i];

    n = 0;
    neighptr = ipage->vget();
    if (fix_history) {
      nn = 0;
      touchptr = ipage_touch->vget();
      shearptr = dpage_shear->vget();
    }

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

      neighptr[n] = joriginal;

      // no numeric test for current touch
      // just use FSH partner list to infer it
      // would require distance calculation for spheres
      // more complex calculation for surfs

      if (fix_history) {
        for (m = 0; m < npartner[i]; m++)
          if (partner[i][m] == jtag) break;
        if (m < npartner[i]) {
          touchptr[n] = 1;
          memcpy(&shearptr[nn],&shearpartner[i][dnum*m],dnumbytes);
          nn += dnum;
        } else {
          touchptr[n] = 0;
          memcpy(&shearptr[nn],zeroes,dnumbytes);
          nn += dnum;
        }
      }

      n++;
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (fix_history) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
      ipage_touch->vgot(n);
      dpage_shear->vgot(nn);
    }
  }

  list->inum = inum;
}
