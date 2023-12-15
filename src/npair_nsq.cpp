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

#include "npair_nsq.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace NeighConst;

/* ---------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE>
NPairNsq<HALF, NEWTON, TRI, SIZE>::NPairNsq(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   Full:
     N^2 search for all neighbors
     every neighbor pair appears in list of both atoms i and j
   Half + Newtoff:
     N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
     pair stored once if i,j are both owned and i < j
     pair stored by me if j is ghost (also stored by proc owning j)
   Half + Newton:
     N^2 / 2 search for neighbor pairs with full Newton's 3rd law
     every pair stored exactly once by some processor
     decision on ghost atoms based on itag, jtag tests
   Half + Newton + Tri:
     use itag/jtap comparision to eliminate half the interactions
     for triclinic, must use delta to eliminate half the I/J interactions
     cannot use I/J exact coord comparision as for orthog
     b/c transforming orthog -> lambda -> orthog for ghost atoms
     with an added PBC offset can shift all 3 coords by epsilon
------------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE>
void NPairNsq<HALF, NEWTON, TRI, SIZE>::build(NeighList *list)
{
  int i, j, jh, jstart, n, itype, jtype, which, bitmask, imol, iatom, moltemplate;
  tagint itag, jtag, tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, radsum, cut, cutsq;
  int *neighptr;

  const double delta = 0.01 * force->angstrom;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) {
    nlocal = atom->nfirst;
    bitmask = group->bitmask[includegroup];
  }

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  if (molecular == Atom::TEMPLATE)
    moltemplate = 1;
  else
    moltemplate = 0;

  int history = list->history;
  int mask_history = 1 << HISTBITS;

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

    // Full: loop over all atoms, owned and ghost, skip i = j
    // Half: loop over remaining atoms, owned and ghost
    //   Newtoff: only store pair if i < j
    //   Newton: itag = jtag is possible for long cutoffs that include images of self

    if (!HALF) jstart = 0;
    else jstart = i + 1;

    for (j = jstart; j < nall; j++) {
      if (includegroup && !(mask[j] & bitmask)) continue;

      if (!HALF) {
        // Full neighbor list
        if (i == j) continue;
      } else if (NEWTON) {
        // Half neighbor list, newton on
        if (j >= nlocal) {
          jtag = tag[j];
          if (itag > jtag) {
            if ((itag + jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag + jtag) % 2 == 1) continue;
          } else if (TRI) {
            if (fabs(x[j][2] - ztmp) > delta) {
              if (x[j][2] < ztmp) continue;
            } else if (fabs(x[j][1] - ytmp) > delta) {
              if (x[j][1] < ytmp) continue;
            } else {
              if (x[j][0] < xtmp) continue;
            }
          } else {
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

      if (SIZE) {
        radsum = radius[i] + radius[j];
        cut = radsum + skin;
        cutsq = cut * cut;

        if (rsq <= cutsq) {
          jh = j;
          if (history && rsq < radsum * radsum) jh = jh ^ mask_history;

          if (molecular != Atom::ATOMIC) {
            if (!moltemplate)
              which = find_special(special[i], nspecial[i], tag[j]);
            else if (imol >= 0)
              which = find_special(onemols[imol]->special[iatom], onemols[imol]->nspecial[iatom],
                                   tag[j] - tagprev);
            else
              which = 0;
            if (which == 0)
              neighptr[n++] = jh;
            else if (domain->minimum_image_check(delx, dely, delz))
              neighptr[n++] = jh;
            else if (which > 0)
              neighptr[n++] = jh ^ (which << SBBITS);
          } else
            neighptr[n++] = jh;
        }
      } else {
        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular != Atom::ATOMIC) {
            if (!moltemplate)
              which = find_special(special[i], nspecial[i], tag[j]);
            else if (imol >= 0)
              which = find_special(onemols[imol]->special[iatom], onemols[imol]->nspecial[iatom],
                                   tag[j] - tagprev);
            else
              which = 0;
            if (which == 0)
              neighptr[n++] = j;
            else if (domain->minimum_image_check(delx, dely, delz))
              neighptr[n++] = j;
            else if (which > 0)
              neighptr[n++] = j ^ (which << SBBITS);
          } else
            neighptr[n++] = j;
        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
  if (!HALF) list->gnum = 0;
}

namespace LAMMPS_NS {
template class NPairNsq<0,1,0,0>;
template class NPairNsq<1,0,0,0>;
template class NPairNsq<1,1,0,0>;
template class NPairNsq<1,1,1,0>;
template class NPairNsq<0,1,0,1>;
template class NPairNsq<1,0,0,1>;
template class NPairNsq<1,1,0,1>;
template class NPairNsq<1,1,1,1>;
}
