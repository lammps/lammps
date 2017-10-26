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

#include "npair_half_respa_nsq_newtoff_omp.h"
#include "npair_omp.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "group.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfRespaNsqNewtoffOmp::NPairHalfRespaNsqNewtoffOmp(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   multiple respa lists
   N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
   pair added to list if atoms i and j are both owned and i < j
   pair added if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void NPairHalfRespaNsqNewtoffOmp::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == 2) ? 1 : 0;

  NPAIR_OMP_INIT;

  const int respamiddle = list->respamiddle;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i,j,n,itype,jtype,n_inner,n_middle,imol,iatom;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*neighptr_inner,*neighptr_middle;

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nall = atom->nlocal + atom->nghost;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int *ilist_inner = list->ilist_inner;
  int *numneigh_inner = list->numneigh_inner;
  int **firstneigh_inner = list->firstneigh_inner;

  int *ilist_middle,*numneigh_middle,**firstneigh_middle;
  if (respamiddle) {
    ilist_middle = list->ilist_middle;
    numneigh_middle = list->numneigh_middle;
    firstneigh_middle = list->firstneigh_middle;
  }

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  MyPage<int> &ipage_inner = list->ipage_inner[tid];
  ipage.reset();
  ipage_inner.reset();

  MyPage<int> *ipage_middle;
  if (respamiddle) {
    ipage_middle = list->ipage_middle + tid;
    ipage_middle->reset();
  }

  int which = 0;
  int minchange = 0;

  for (i = ifrom; i < ito; i++) {

    n = n_inner = 0;
    neighptr = ipage.vget();
    neighptr_inner = ipage_inner.vget();
    if (respamiddle) {
      n_middle = 0;
      neighptr_middle = ipage_middle->vget();
    }

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over remaining atoms, owned and ghost

    for (j = i+1; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;
      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
        if (molecular) {
          if (!moltemplate)
            which = find_special(special[i],nspecial[i],tag[j]);
          else if (imol >=0)
            which = find_special(onemols[imol]->special[iatom],
                                 onemols[imol]->nspecial[iatom],
                                 tag[j]-tagprev);
          else which = 0;
          if (which == 0) neighptr[n++] = j;
          else if ((minchange = domain->minimum_image_check(delx,dely,delz)))
            neighptr[n++] = j;
          else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
        } else neighptr[n++] = j;

        if (rsq < cut_inner_sq) {
          if (which == 0) neighptr_inner[n_inner++] = j;
          else if (minchange) neighptr_inner[n_inner++] = j;
          else if (which > 0) neighptr_inner[n_inner++] = j ^ (which << SBBITS);
        }

        if (respamiddle && rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
          if (which == 0) neighptr_middle[n_middle++] = j;
          else if (minchange) neighptr_middle[n_middle++] = j;
          else if (which > 0)
            neighptr_middle[n_middle++] = j ^ (which << SBBITS);
        }
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    ilist_inner[i] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    ipage.vgot(n_inner);
    if (ipage_inner.status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[i] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      ipage_middle->vgot(n_middle);
      if (ipage_middle->status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }
  }
  NPAIR_OMP_CLOSE;
  list->inum = nlocal;
  list->inum_inner = nlocal;
  if (respamiddle) list->inum_middle = nlocal;
}
