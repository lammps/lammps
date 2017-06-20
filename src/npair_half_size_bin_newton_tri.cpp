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
#include "npair_half_size_bin_newton_tri.h"
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

NPairHalfSizeBinNewtonTri::NPairHalfSizeBinNewtonTri(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   size particles
   binned neighbor list construction with Newton's 3rd law for triclinic
   shear history must be accounted for when a neighbor pair is added
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfSizeBinNewtonTri::build(NeighList *list)
{
  int i,j,k,m,n,nn,ibin,dnum,dnumbytes;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;

  int *npartner;
  tagint **partner;
  double **shearpartner;
  int **firsttouch;
  double **firstshear;
  MyPage<int> *ipage_touch;
  MyPage<double> *dpage_shear;
  NeighList *listhistory;

  double **x = atom->x;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

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

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (fix_history) {
      nn = 0;
      touchptr = ipage_touch->vget();
      shearptr = dpage_shear->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

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

        if (exclude && exclusion(i,j,type[i],type[j],mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radi + radius[j];
        cutsq = (radsum+skin) * (radsum+skin);

        if (rsq <= cutsq) {
          neighptr[n++] = j;

          if (fix_history) {
            if (rsq < radsum*radsum) {
              for (m = 0; m < npartner[i]; m++)
                if (partner[i][m] == tag[j]) break;
              if (m < npartner[i]) {
                touchptr[n] = 1;
                memcpy(&shearptr[nn],&shearpartner[i][dnum*m],dnumbytes);
                nn += dnum;
              } else {
                touchptr[n] = 0;
                memcpy(&shearptr[nn],zeroes,dnumbytes);
                nn += dnum;
              }
            } else {
              touchptr[n] = 0;
              memcpy(&shearptr[nn],zeroes,dnumbytes);
              nn += dnum;
            }
          }

          n++;
        }
      }
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
