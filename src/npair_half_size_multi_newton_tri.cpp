// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_half_size_multi_newton_tri.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "molecule.h"
#include "my_page.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfSizeMultiNewtonTri::NPairHalfSizeMultiNewtonTri(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   size particles
   binned neighbor list construction with Newton's 3rd law for triclinic
   multi stencil is icollection-jcollection dependent
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfSizeMultiNewtonTri::build(NeighList *list)
{
  int i,j,jh,k,n,itype,jtype,icollection,jcollection,ibin,jbin,ns,js;
  int which,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutdistsq;
  int *neighptr,*s;

  int *collection = neighbor->collection;
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;

  int history = list->history;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int mask_history = 1 << HISTBITS;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    itype = type[i];
    icollection = collection[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    ibin = atom2bin[i];

    // loop through stencils for all collections
    for (jcollection = 0; jcollection < ncollections; jcollection++) {

      // if same collection use own bin
      if (icollection == jcollection) jbin = ibin;
          else jbin = coord2bin(x[i], jcollection);

      // loop over all atoms in bins in stencil
      // stencil is empty if i larger than j
      // stencil is half if i same size as j
      // stencil is full if i smaller than j
      // if half: pairs for atoms j "below" i are excluded
      // below = lower z or (equal z and lower y) or (equal zy and lower x)
      //         (equal zyx and j <= i)
      // latter excludes self-self interaction but allows superposed atoms

          s = stencil_multi[icollection][jcollection];
          ns = nstencil_multi[icollection][jcollection];

          for (k = 0; k < ns; k++) {
            js = binhead_multi[jcollection][jbin + s[k]];
            for (j = js; j >= 0; j = bins[j]) {

          // if same size (same collection), use half stencil
          if (cutcollectionsq[icollection][icollection] == cutcollectionsq[jcollection][jcollection]){
            if (x[j][2] < ztmp) continue;
            if (x[j][2] == ztmp) {
              if (x[j][1] < ytmp) continue;
              if (x[j][1] == ytmp) {
                if (x[j][0] < xtmp) continue;
                if (x[j][0] == xtmp && j <= i) continue;
              }
            }
          }

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              radsum = radi + radius[j];
              cutdistsq = (radsum+skin) * (radsum+skin);

              if (rsq <= cutdistsq) {
                jh = j;
                if (history && rsq < radsum*radsum)
                  jh = jh ^ mask_history;

                if (molecular != Atom::ATOMIC) {
                  if (!moltemplate)
                    which = find_special(special[i],nspecial[i],tag[j]);
                  else if (imol >= 0)
                    which = find_special(onemols[imol]->special[iatom],
                                         onemols[imol]->nspecial[iatom],
                                         tag[j]-tagprev);
                  else which = 0;
                  if (which == 0) neighptr[n++] = jh;
                  else if (domain->minimum_image_check(delx,dely,delz))
                    neighptr[n++] = jh;
                  else if (which > 0) neighptr[n++] = jh ^ (which << SBBITS);
                } else neighptr[n++] = jh;
              }
            }
          }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}
