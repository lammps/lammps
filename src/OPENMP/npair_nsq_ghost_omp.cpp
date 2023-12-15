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


#include "npair_nsq_ghost_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "atom_vec.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF>
NPairNsqGhostOmp<HALF>::NPairNsqGhostOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   Full:
     N^2 search for all neighbors
     include neighbors of ghost atoms, but no "special neighbors" for ghosts
     every neighbor pair appears in list of both atoms i and j
   Half + Newtoff:
     N^2 / 2 search for neighbor pairs with partial Newton's 3rd law
     include neighbors of ghost atoms, but no "special neighbors" for ghosts
     pair stored once if i,j are both owned and i < j
     pair stored by me if i owned and j ghost (also stored by proc owning j)
     pair stored once if i,j are both ghost and i < j
------------------------------------------------------------------------- */

template<int HALF>
void NPairNsqGhostOmp<HALF>::build(NeighList *list)
{
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nall);

  int i, j, jstart, n, itype, jtype, which, imol, iatom;
  tagint tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
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

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  // loop over owned & ghost atoms, storing neighbors

  for (i = ifrom; i < ito; i++) {

    n = 0;
    neighptr = ipage.vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over all atoms, owned and ghost
    // Full:
    //   skip i = j
    // Half:
    //   only store pair if i < j
    //   stores own/own pairs only once
    //   stores own/ghost pairs with owned atom only, on both procs
    //   stores ghost/ghost pairs only once
    // no molecular test when i = ghost atom

    if (HALF) jstart = i + 1;
    else jstart = 0;

    if (i < nlocal) {
      for (j = jstart; j < nall; j++) {
        if (!HALF) {
          if (i == j) continue;
        }

        jtype = type[j];
        if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
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
    } else {
      for (j = jstart; j < nall; j++) {
        if (!HALF) {
          if (i == j) continue;
        }

        jtype = type[j];
        if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (HALF) {
          if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
        } else {
          if (rsq <= cutneighghostsq[itype][jtype]) neighptr[n++] = j;
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
  list->gnum = nall - nlocal;
}

namespace LAMMPS_NS {
template class NPairNsqGhostOmp<0>;
template class NPairNsqGhostOmp<1>;
}
