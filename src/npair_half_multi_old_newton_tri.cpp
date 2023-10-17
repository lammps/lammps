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

#include "npair_half_multi_old_newton_tri.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfMultiOldNewtonTri::NPairHalfMultiOldNewtonTri(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   multi-type stencil is itype dependent and is distance checked
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfMultiOldNewtonTri::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which,ns,imol,iatom,moltemplate;
  tagint itag,jtag,tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  double *cutsq,*distsq;

  const double delta = 0.01 * force->angstrom;

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

    itag = tag[i];
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
    // for triclinic, bin stencil is full in all 3 dims
    // must use itag/jtag to eliminate half the I/J interactions
    // cannot use I/J exact coord comparision
    //   b/c transforming orthog -> lambda -> orthog for ghost atoms
    //   with an added PBC offset can shift all 3 coords by epsilon

    ibin = atom2bin[i];
    s = stencil_multi_old[itype];
    distsq = distsq_multi_old[itype];
    cutsq = cutneighsq[itype];
    ns = nstencil_multi_old[itype];
    for (k = 0; k < ns; k++) {
      for (j = binhead[ibin+s[k]]; j >= 0; j = bins[j]) {
        jtype = type[j];
        if (cutsq[jtype] < distsq[k]) continue;

        if (j <= i) continue;
        if (j >= nlocal) {
          jtag = tag[j];
          if (itag > jtag) {
            if ((itag+jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag+jtag) % 2 == 1) continue;
          } else {
            if (fabs(x[j][2]-ztmp) > delta) {
              if (x[j][2] < ztmp) continue;
            } else if (fabs(x[j][1]-ytmp) > delta) {
              if (x[j][1] < ytmp) continue;
            } else {
              if (x[j][0] < xtmp) continue;
            }
          }
        }

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
