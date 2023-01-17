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

#include "npair_half_bin_newton_tri.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonTri::NPairHalfBinNewtonTri(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfBinNewtonTri::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  double **x = atom->x;
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

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over all atoms in bins in stencil
    // pairs for atoms j "below" i are excluded
    // below = lower z or (equal z and lower y) or (equal zy and lower x)
    //         (equal zyx and j <= i)
    // latter excludes self-self interaction but allows superposed atoms

    ibin = atom2bin[i];
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp) {
          if (x[j][1] < ytmp) continue;
          if (x[j][1] == ytmp) {
            if (x[j][0] < xtmp) continue;
            if (x[j][0] == xtmp && j <= i) continue;
          }
        }

        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular != Atom::ATOMIC) {
            if (!moltemplate)
              which = find_special(special[i],nspecial[i],tag[j]);
            else if (imol >= 0)
              which = find_special(onemols[imol]->special[iatom],
                                   onemols[imol]->nspecial[iatom],
                                   tag[j]-tagprev);
            else which = 0;
            if (which == 0) neighptr[n++] = j;
            else if (domain->minimum_image_check(delx,dely,delz))
              neighptr[n++] = j;
            else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
          } else neighptr[n++] = j;
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
