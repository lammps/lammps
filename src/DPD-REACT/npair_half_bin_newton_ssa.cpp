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

/* ----------------------------------------------------------------------
   Contributing authors:
   James Larentzos and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "npair_half_bin_newton_ssa.h"
#include "nstencil_ssa.h"
#include "nbin_ssa.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "memory.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonSSA::NPairHalfBinNewtonSSA(LAMMPS *lmp) : NPair(lmp)
{
  ssa_maxPhaseCt = 0;
  ssa_maxPhaseLen = 0;
  ssa_phaseCt = 0;
  ssa_phaseLen = nullptr;
  ssa_itemLoc = nullptr;
  ssa_itemLen = nullptr;
  ssa_gphaseCt = 7;
  memory->create(ssa_gphaseLen,ssa_gphaseCt,"NPairHalfBinNewtonSSA:ssa_gphaseLen");
  memory->create(ssa_gitemLoc,ssa_gphaseCt,1,"NPairHalfBinNewtonSSA:ssa_gitemLoc");
  memory->create(ssa_gitemLen,ssa_gphaseCt,1,"NPairHalfBinNewtonSSA:ssa_gitemLen");
}

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonSSA::~NPairHalfBinNewtonSSA()
{
  ssa_maxPhaseCt = 0;
  ssa_maxPhaseLen = 0;
  ssa_phaseCt = 0;
  memory->destroy(ssa_phaseLen);
  memory->destroy(ssa_itemLoc);
  memory->destroy(ssa_itemLen);
  ssa_gphaseCt = 0;
  memory->destroy(ssa_gphaseLen);
  memory->destroy(ssa_gitemLoc);
  memory->destroy(ssa_gitemLen);
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   for use by Shardlow Splitting Algorithm
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

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  int molecular = atom->molecular;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  auto ns_ssa = dynamic_cast<NStencilSSA*>(ns);
  if (!ns_ssa) error->one(FLERR, "NStencil wasn't a NStencilSSA object");
  int *nstencil_ssa = &(ns_ssa->nstencil_ssa[0]);
  int nstencil_full = ns_ssa->nstencil;

  auto nb_ssa = dynamic_cast<NBinSSA*>(nb);
  if (!nb_ssa) error->one(FLERR, "NBin wasn't a NBinSSA object");
  int *bins = nb_ssa->bins;
  int *binhead = nb_ssa->binhead;
  int *gairhead_ssa = &(nb_ssa->gairhead_ssa[0]);

  int inum = 0;
  int gnum = 0;
  int xbin,ybin,zbin,xbin2,ybin2,zbin2;
  int **stencilxyz = ns_ssa->stencilxyz;
  int lbinxlo = nb_ssa->lbinxlo;
  int lbinxhi = nb_ssa->lbinxhi;
  int lbinylo = nb_ssa->lbinylo;
  int lbinyhi = nb_ssa->lbinyhi;
  int lbinzlo = nb_ssa->lbinzlo;
  int lbinzhi = nb_ssa->lbinzhi;

  int sx1 = ns_ssa->sx + 1;
  int sy1 = ns_ssa->sy + 1;
  int sz1 = ns_ssa->sz + 1;

  ssa_phaseCt = sz1*sy1*sx1;

  xbin = (lbinxhi - lbinxlo + sx1 - 1) / sx1 + 1;
  ybin = (lbinyhi - lbinylo + sy1 - 1) / sy1 + 1;
  zbin = (lbinzhi - lbinzlo + sz1 - 1) / sz1 + 1;

  int phaseLenEstimate = xbin*ybin*zbin;

  if (ssa_phaseCt > ssa_maxPhaseCt) {
    ssa_maxPhaseCt = ssa_phaseCt;
    ssa_maxPhaseLen = 0;
    memory->destroy(ssa_phaseLen);
    memory->destroy(ssa_itemLoc);
    memory->destroy(ssa_itemLen);
    memory->create(ssa_phaseLen,ssa_maxPhaseCt,"NPairHalfBinNewtonSSA:ssa_phaseLen");
  }

  if (phaseLenEstimate > ssa_maxPhaseLen) {
    ssa_maxPhaseLen = phaseLenEstimate;
    memory->destroy(ssa_itemLoc);
    memory->destroy(ssa_itemLen);
    memory->create(ssa_itemLoc,ssa_maxPhaseCt,ssa_maxPhaseLen,"NPairHalfBinNewtonSSA:ssa_itemLoc");
    memory->create(ssa_itemLen,ssa_maxPhaseCt,ssa_maxPhaseLen,"NPairHalfBinNewtonSSA:ssa_itemLen");
  }

  ipage->reset();

  int workPhase = 0;
  // loop over bins with local atoms, storing half of the neighbors
  for (int zoff = ns_ssa->sz; zoff >= 0; --zoff) {
  for (int yoff = ns_ssa->sy; yoff >= 0; --yoff) {
  for (int xoff = ns_ssa->sx; xoff >= 0; --xoff) {
    int workItem = 0;
  for (zbin = lbinzlo + zoff; zbin < lbinzhi; zbin += sz1) {
  for (ybin = lbinylo + yoff - ns_ssa->sy; ybin < lbinyhi; ybin += sy1) {
  for (xbin = lbinxlo + xoff - ns_ssa->sx; xbin < lbinxhi; xbin += sx1) {
    if (workItem >= phaseLenEstimate) error->one(FLERR,"phaseLenEstimate was too small");
    ssa_itemLoc[workPhase][workItem] = inum; // record where workItem starts in ilist

    for (int subphase = 0; subphase < 4; subphase++) {
      int s_ybin = ybin + ((subphase & 0x2) ? ns_ssa->sy : 0);
      int s_xbin = xbin + ((subphase & 0x1) ? ns_ssa->sx : 0);
      int ibin;

      if ((s_ybin < lbinylo) || (s_ybin >= lbinyhi)) continue;
      if ((s_xbin < lbinxlo) || (s_xbin >= lbinxhi)) continue;
      ibin = zbin*nb_ssa->mbiny*nb_ssa->mbinx
           + s_ybin*nb_ssa->mbinx
           + s_xbin;

      for (i = binhead[ibin]; i >= 0; i = bins[i]) {
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

        // loop over all local atoms in the current stencil "subphase"
        for (k = nstencil_ssa[subphase]; k < nstencil_ssa[subphase+1]; k++) {
          const int jbin = ibin+stencil[k];
          if (jbin != ibin) j = binhead[jbin];
          else j = bins[i]; // same bin as i, so start just past i in the bin
          for (; j >= 0; j = bins[j]) {
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

        if (n > 0) {
          firstneigh[inum] = neighptr;
          numneigh[inum] = n;
          ilist[inum++] = i;
        }
        ipage->vgot(n);
        if (ipage->status())
          error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
      }
    }
    // record where workItem ends in ilist
    ssa_itemLen[workPhase][workItem] = inum - ssa_itemLoc[workPhase][workItem];
    if (ssa_itemLen[workPhase][workItem] > 0) workItem++;
  }
  }
  }

    // record where workPhase ends
    ssa_phaseLen[workPhase++] = workItem;
  }
  }
  }

  if (ssa_phaseCt != workPhase) error->one(FLERR,"ssa_phaseCt was wrong");

  list->inum = inum;

  // loop over AIR ghost atoms, storing their local neighbors
  // since these are ghosts, must check if stencil bin is out of bounds
  for (workPhase = 0; workPhase < ssa_gphaseCt; workPhase++) {
    int locAIRct = 0;
    ssa_gitemLoc[workPhase][0] = inum + gnum; // record where workItem starts in ilist
    for (i = gairhead_ssa[workPhase+1]; i >= 0; i = bins[i]) {
      n = 0;
      neighptr = ipage->vget();

      itype = type[i];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      ibin = coord2bin(x[i],xbin,ybin,zbin);

      // loop over AIR ghost atoms in all bins in "full" stencil
      // Note: the non-AIR ghost atoms have already been filtered out
      for (k = 0; k < nstencil_full; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        // Skip it if this bin is outside the extent of local bins
        if (xbin2 < lbinxlo || xbin2 >= lbinxhi ||
            ybin2 < lbinylo || ybin2 >= lbinyhi ||
            zbin2 < lbinzlo || zbin2 >= lbinzhi) continue;
        for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {

          jtype = type[j];
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq <= cutneighsq[itype][jtype]) {
            if (molecular != Atom::ATOMIC) {
              if (!moltemplate)
                which = find_special(special[j],nspecial[j],tag[i]);
              else {
                int jmol = molindex[j];
                if (jmol >= 0) {
                  int jatom = molatom[j];
                  which = find_special(onemols[jmol]->special[jatom],
                                     onemols[jmol]->nspecial[jatom],
                                     tag[i] - (tag[j] - jatom - 1));
                } else which = 0;
              }
              if (which == 0) neighptr[n++] = j;
              else if (domain->minimum_image_check(delx,dely,delz))
                neighptr[n++] = j;
              else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
            } else neighptr[n++] = j;
          }
        }
      }

      if (n > 0) {
        firstneigh[inum + gnum] = neighptr;
        numneigh[inum + gnum] = n;
        ilist[inum + (gnum++)] = i;
        ++locAIRct;
      }
      ipage->vgot(n);
      if (ipage->status())
        error->one(FLERR,"Neighbor (ghost) list overflow, boost neigh_modify one");
    }
    ssa_gitemLen[workPhase][0] = locAIRct;
    ssa_gphaseLen[workPhase] = 1;
  }
  list->gnum = gnum;
}
