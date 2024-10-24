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

#include "npair_multi_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE, int CUSTOMCHECK, int ATOMONLY>
NPairMultiOmp<HALF, NEWTON, TRI, SIZE, CUSTOMCHECK, ATOMONLY>::NPairMultiOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   multi stencil is icollection-jcollection dependent
   Full:
     binned neighbor list construction for all neighbors
     every neighbor pair appears in list of both atoms i and j
   Half + Newtoff:
     binned neighbor list construction with partial Newton's 3rd law
     each owned atom i checks own bin and other bins in stencil
     pair stored once if i,j are both owned and i < j
     pair stored by me if j is ghost (also stored by proc owning j)
   Half + Newton:
     binned neighbor list construction with full Newton's 3rd law
     each owned atom i checks its own bin and other bins in Newton stencil
     every pair stored exactly once by some processor
------------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE, int CUSTOMCHECK, int ATOMONLY>
void NPairMultiOmp<HALF, NEWTON, TRI, SIZE, CUSTOMCHECK, ATOMONLY>::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  const double delta = 0.01 * force->angstrom;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i, j, jh, js, k, n, itype, jtype, ibin, jbin, icollection, jcollection, which, ns, imol, iatom;
  tagint itag, jtag, tagprev, neigh_check;
  double xtmp, ytmp, ztmp, rtmp, delx, dely, delz, rsq, radsum, cut, cutsq;
  int *neighptr, *s;

  int *collection = neighbor->collection;
  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int history = list->history;
  int mask_history = 1 << HISTBITS;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itag = tag[i];
    itype = type[i];
    icollection = collection[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    if (!ATOMONLY) {
      if (moltemplate) {
        imol = molindex[i];
        iatom = molatom[i];
        tagprev = tag[i] - iatom - 1;
      }
    }

    if (SIZE)
      rtmp = radius[i];

    ibin = atom2bin[i];

    // loop through stencils for all collections

    for (jcollection = 0; jcollection < ncollections; jcollection++) {

      // Use own bin for same collection
      if (icollection == jcollection) jbin = ibin;
      else jbin = coord2bin(x[i], jcollection);

      s = stencil_multi[icollection][jcollection];
      ns = nstencil_multi[icollection][jcollection];

      for (k = 0; k < ns; k++) {
        js = binhead_multi[jcollection][jbin + s[k]];

        // For half-newton-ortho, first check self bin (k == 0, always half)
        // if checking its own binlist, skip all before i in linked list
        if (HALF && NEWTON && !TRI)
          if ((k == 0) && (icollection == jcollection)) js = bins[i];

        for (j = js; j >= 0; j = bins[j]) {
          if (!HALF) {
            // Full neighbor list, only uses full stencils
            // only skip i = j
            if (i == j) continue;
          } else if (!NEWTON) {
            // Half neighbor list, newton off, only uses full stencils
            // only store pair if i < j
            // stores own/own pairs only once
            // stores own/ghost pairs on both procs
            if (j <= i) continue;
          } else if (TRI) {
            // Half neighbor list, newton on, triclinic, only uses full stencils
            // If different sizes -> full stencil (accept all, one-way search)
            // If same size -> half stencil, exclude half of interactions
            //     stencil is empty if i larger than j
            //     stencil is full if i smaller than j
            //     stencil is full if i same size as j
            //   for i smaller than j:
            //     must use itag/jtag to eliminate half the I/J interactions
            //     cannot use I/J exact coord comparision
            //       b/c transforming orthog -> lambda -> orthog for ghost atoms
            //   with an added PBC offset can shift all 3 coords by epsilon

            if (flag_same_multi[icollection][jcollection]) {
              if (j <= i) continue;
              if (j >= nlocal) {
                jtag = tag[j];
                if (itag > jtag) {
                  if ((itag + jtag) % 2 == 0) continue;
                } else if (itag < jtag) {
                  if ((itag + jtag) % 2 == 1) continue;
                } else {
                  if (fabs(x[j][2] - ztmp) > delta) {
                    if (x[j][2] < ztmp) continue;
                  } else if (fabs(x[j][1] - ytmp) > delta) {
                    if (x[j][1] < ytmp) continue;
                  } else {
                    if (x[j][0] < xtmp) continue;
                  }
                }
              }
            }
          } else {
            // Half neighbor list, newton on, orthonormal, uses a mix of stencils
            // If different sizes -> full stencil (accept all, one-way search)
            // If same size -> half stencil (first includes a self bin search)
            if (k == 0 && flag_same_multi[icollection][jcollection]) {
              // if same collection,
              //   if j is owned atom, store it, since j is beyond i in linked list
              //   if j is ghost, only store if j coords are "above and to the right" of i

              // if different collections,
              //   if j is owned atom, store it if j > i
              //   if j is ghost, only store if j coords are "above and to the right" of i

              if ((icollection != jcollection) && (j < i)) continue;

              if (j >= nlocal) {
                if (x[j][2] < ztmp) continue;
                if (x[j][2] == ztmp) {
                  if (x[j][1] < ytmp) continue;
                  if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
                }
              }
            }
          }

          jtype = type[j];
          if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (CUSTOMCHECK) {
            neigh_check = pair->neigh_check(i, j, skin, rsq);
          } else if (SIZE) {
            radsum = rtmp + radius[j];
            cut = radsum + skin;
            cutsq = cut * cut;
            neigh_check = rsq <= cutsq;
          } else {
            neigh_check = rsq <= cutneighsq[itype][jtype];
          }

          if (!neigh_check) continue;

          if (ATOMONLY) {
            if (SIZE && history && (rsq < (radsum * radsum))) {
              neighptr[n++] = j ^ mask_history;
            } else {
              neighptr[n++] = j;
            }
          } else {
            if (molecular != Atom::ATOMIC) {
              if (!moltemplate)
                which = find_special(special[i], nspecial[i], tag[j]);
              else if (imol >= 0)
                which = find_special(onemols[imol]->special[iatom], onemols[imol]  ->nspecial[iatom],
                                     tag[j] - tagprev);
              else
                which = 0;

              if (SIZE && history && (rsq < (radsum * radsum))) {
                if (which == 0)
                  neighptr[n++] = j ^ mask_history;
                else if (domain->minimum_image_check(delx, dely, delz))
                  neighptr[n++] = j ^ mask_history;
                else if (which > 0)
                  neighptr[n++] = (j ^ mask_history) ^ (which <<   SBBITS);
              } else {
                if (which == 0)
                  neighptr[n++] = j;
                else if (domain->minimum_image_check(delx, dely, delz))
                  neighptr[n++] = j;
                else if (which > 0)
                  neighptr[n++] = j ^ (which << SBBITS);
              }
            } else {
              if (SIZE && history && (rsq < (radsum * radsum))) {
                neighptr[n++] = j ^ mask_history;
              } else {
                neighptr[n++] = j;
              }
            }
          }
        }
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
  NPAIR_OMP_CLOSE;
  list->inum = nlocal;
  list->gnum = 0;
}

namespace LAMMPS_NS {
template class NPairMultiOmp<0,1,0,0,0,0>;
template class NPairMultiOmp<1,0,0,0,0,0>;
template class NPairMultiOmp<1,1,0,0,0,0>;
template class NPairMultiOmp<1,1,1,0,0,0>;
template class NPairMultiOmp<0,1,0,1,0,0>;
template class NPairMultiOmp<1,0,0,1,0,0>;
template class NPairMultiOmp<1,1,0,1,0,0>;
template class NPairMultiOmp<1,1,1,1,0,0>;
template class NPairMultiOmp<0,1,0,0,1,0>;
template class NPairMultiOmp<1,0,0,0,1,0>;
template class NPairMultiOmp<1,1,0,0,1,0>;
template class NPairMultiOmp<1,1,1,0,1,0>;
template class NPairMultiOmp<0,1,0,0,0,1>;
template class NPairMultiOmp<1,0,0,0,0,1>;
template class NPairMultiOmp<1,1,0,0,0,1>;
template class NPairMultiOmp<1,1,1,0,0,1>;
template class NPairMultiOmp<0,1,0,1,0,1>;
template class NPairMultiOmp<1,0,0,1,0,1>;
template class NPairMultiOmp<1,1,0,1,0,1>;
template class NPairMultiOmp<1,1,1,1,0,1>;
template class NPairMultiOmp<0,1,0,0,1,1>;
template class NPairMultiOmp<1,0,0,0,1,1>;
template class NPairMultiOmp<1,1,0,0,1,1>;
template class NPairMultiOmp<1,1,1,0,1,1>;
}
