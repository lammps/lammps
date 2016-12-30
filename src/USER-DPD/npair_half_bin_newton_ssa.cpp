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

/* ----------------------------------------------------------------------
   Contributing authors:
   James Larentzos and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "npair_half_bin_newton_ssa.h"
#include "neighbor.h"
#include "nstencil_ssa.h"
#include "nbin_ssa.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonSSA::NPairHalfBinNewtonSSA(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   for use by Shardlow Spliting Algorithm
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfBinNewtonSSA::build(NeighList *list)
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
  int molecular = atom->molecular;
  if (molecular == 2) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  NStencilSSA *ns_ssa = dynamic_cast<NStencilSSA*>(ns);
  if (!ns_ssa) error->one(FLERR, "NStencil wasn't a NStencilSSA object");
  int nstencil_half = ns_ssa->nstencil_half;
  int nstencil_full = ns_ssa->nstencil;

  NBinSSA *nb_ssa = dynamic_cast<NBinSSA*>(nb);
  if (!nb_ssa) error->one(FLERR, "NBin wasn't a NBinSSA object");
  int *bins_ssa = nb_ssa->bins_ssa;
  int *binhead_ssa = nb_ssa->binhead_ssa;
  int *gairhead_ssa = &(nb_ssa->gairhead_ssa[0]);

  int inum = 0;
  int gnum = 0;
  int xbin,ybin,zbin,xbin2,ybin2,zbin2;
  int **stencilxyz = ns_ssa->stencilxyz;

  ipage->reset();

  // loop over owned atoms, storing half of the neighbors

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

    // loop over rest of local atoms in i's bin
    // just store them, since j is beyond i in linked list

    for (j = bins_ssa[i]; j >= 0; j = bins_ssa[j]) {

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

    ibin = coord2bin(x[i]);

    // loop over all local atoms in other bins in "half" stencil

    for (k = 0; k < nstencil_half; k++) {
      for (j = binhead_ssa[ibin+stencil[k]]; j >= 0;
           j = bins_ssa[j]) {

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

  // loop over AIR ghost atoms, storing their local neighbors
  // since these are ghosts, must check if stencil bin is out of bounds
  for (int airnum = 2; airnum <= 8; airnum++) {
    for (i = gairhead_ssa[airnum]; i >= 0; i = bins_ssa[i]) {
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

      ibin = coord2bin(x[i],xbin,ybin,zbin);

      // loop over AIR ghost atoms in all bins in "full" stencil
      // Note: the non-AIR ghost atoms have already been filtered out
      for (k = 0; k < nstencil_full; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        // since we only care about ghost to local neighbors, these "bounds" could be inset
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) continue;
        for (j = binhead_ssa[ibin+stencil[k]]; j >= 0; j = bins_ssa[j]) {

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

      if (n > 0) ilist[inum + (gnum++)] = i;
      firstneigh[i] = neighptr;
      numneigh[i] = n;
      ipage->vgot(n);
      if (ipage->status())
        error->one(FLERR,"Neighbor (ghost) list overflow, boost neigh_modify one");
    }
  }
  list->gnum = gnum;
}
