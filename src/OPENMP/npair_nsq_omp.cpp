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

#include "npair_nsq_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

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

#include <cmath>

using namespace LAMMPS_NS;
using namespace NeighConst;

/* ---------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE>
NPairNsqOmp<HALF, NEWTON, TRI, SIZE>::NPairNsqOmp(LAMMPS *lmp) : NPair(lmp) {}

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
     decision on ghost atoms based on itag,jtag tests
   Half + Newton + Tri:
     use itag/jtap comparision to eliminate half the interactions
     for triclinic, must use delta to eliminate half the I/J interactions
     cannot use I/J exact coord comparision as for orthog
     b/c transforming orthog -> lambda -> orthog for ghost atoms
     with an added PBC offset can shift all 3 coords by epsilon
------------------------------------------------------------------------- */

template<int HALF, int NEWTON, int TRI, int SIZE>
void NPairNsqOmp<HALF, NEWTON, TRI, SIZE>::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int bitmask = (includegroup) ? group->bitmask[includegroup] : 0;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  const double delta = 0.01 * force->angstrom;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i, j, jh, jstart, n, itype, jtype, which, imol, iatom;
  tagint itag, jtag, tagprev, neigh_check;
  double xtmp, ytmp, ztmp, rtmp, delx, dely, delz, rsq, radsum, cut, cutsq;
  int *neighptr;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;

  int nall = atom->nlocal + atom->nghost;
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

  // loop over owned atoms, storing neighbors

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

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

    if (SIZE)
      rtmp = radius[i];

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
        radsum = rtmp + radius[j];
        cut = radsum + skin;
        cutsq = cut * cut;
        neigh_check = rsq <= cutsq;
      } else {
        neigh_check = rsq <= cutneighsq[itype][jtype];
      }

      if (!neigh_check) continue;

      if (molecular != Atom::ATOMIC) {
        if (!moltemplate)
          which = find_special(special[i], nspecial[i], tag[j]);
        else if (imol >= 0)
          which = find_special(onemols[imol]->special[iatom], onemols[imol]->nspecial[iatom],
                               tag[j] - tagprev);
        else
          which = 0;

        if (SIZE && history && (rsq < (radsum * radsum))) {
          if (which == 0)
            neighptr[n++] = j ^ mask_history;
          else if (domain->minimum_image_check(delx, dely, delz))
            neighptr[n++] = j ^ mask_history;
          else if (which > 0)
            neighptr[n++] = (j ^ mask_history) ^ (which << SBBITS);
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
template class NPairNsqOmp<0,1,0,0>;
template class NPairNsqOmp<1,0,0,0>;
template class NPairNsqOmp<1,1,0,0>;
template class NPairNsqOmp<1,1,1,0>;
template class NPairNsqOmp<0,1,0,1>;
template class NPairNsqOmp<1,0,0,1>;
template class NPairNsqOmp<1,1,0,1>;
template class NPairNsqOmp<1,1,1,1>;
}
