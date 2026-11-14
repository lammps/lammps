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

#include "npair_bin_ghost_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int HALF>
NPairBinGhostOmp<HALF>::NPairBinGhostOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   Full:
     binned neighbor list construction for all neighbors
     include neighbors of ghost atoms, but no "special neighbors" for ghosts
     every neighbor pair appears in list of both atoms i and j
   Half + Newtoff:
     binned neighbor list construction with partial Newton's 3rd law
     include neighbors of ghost atoms, but no "special neighbors" for ghosts
     owned and ghost atoms check own bin and other bins in stencil
     pair stored once if i,j are both owned and i < j
     pair stored by me if i owned and j ghost (also stored by proc owning j)
     pair stored once if i,j are both ghost and i < j
------------------------------------------------------------------------- */

template<int HALF>
void NPairBinGhostOmp<HALF>::build(NeighList *list)
{
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  int overflow = 0;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list,fix_package_omp) reduction(+:overflow)
#endif
  NPAIR_OMP_SETUP(nall);

  int i, j, k, n, itype, jtype, ibin, which, imol, iatom;
  tagint tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int xbin, ybin, zbin, xbin2, ybin2, zbin2;
  int *neighptr;

  const double (*const x)[3] = (const double (*)[3]) atom->x[0];
  const int *const type = atom->type;
  int *mask = atom->mask;
  const tagint *const tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  const bool atomic = (molecular == Atom::ATOMIC);

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  const int oneatom = list->oneatom;

  const int *const bins = this->bins;
  const int *const binhead = this->binhead;
  const int *const stencil = this->stencil;
  const int nstencil = this->nstencil;

  // whether type t has the same neighbor cutoff to all types, see copy_neighbor_info()
  const int *const ucut_type = this->ucut_type;
  const int *const ucut_ghost_type = this->ucut_ghost_type;

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

    // loop over all atoms in surrounding bins in stencil including self
    // when i is a ghost atom, must check if stencil bin is out of bounds
    // no molecular test when i = ghost atom

    if (i < nlocal) {
      ibin = atom2bin[i];
      const double *const cutneighsq_i = cutneighsq[itype];
      const double cutsq_u = cutneighsq_i[1];
      const bool uniform_cut = !exclude && (ucut_type[itype] != 0);

      for (k = 0; k < nstencil; k++) {
        for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
          if (HALF) {
            // Half neighbor list, newton off
            // only store pair if i < j
            // stores own/own pairs only once
            // stores own/ghost pairs on both procs
            // stores ghost/ghost pairs only once
            if (j <= i) continue;
          } else {
            // Full neighbor list
            // only skip i = j
            if (i == j) continue;
          }

          // skip type[j] load when !exclude and all cutoffs from itype match
          if (!uniform_cut) {
            jtype = type[j];
            if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;
          }

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          const double cutsq_ij = uniform_cut ? cutsq_u : cutneighsq_i[jtype];
          if (atomic) {
            if (n < oneatom) {
              neighptr[n] = j;
              n += (rsq <= cutsq_ij);
            } else if (rsq <= cutsq_ij)
              n++;
          } else if (rsq <= cutsq_ij) {
            if (!moltemplate)
              which = find_special(special[i], nspecial[i], tag[j]);
            else if ((imol >= 0) && onemols[imol]->special)
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
          }
        }
      }

    } else {
      ibin = coord2bin(atom->x[i], xbin, ybin, zbin);
      const double *const cutneighghostsq_i = cutneighghostsq[itype];
      const double cutsq_u = cutneighghostsq_i[1];
      const bool uniform_cut = !exclude && (ucut_ghost_type[itype] != 0);

      for (k = 0; k < nstencil; k++) {
        xbin2 = xbin + stencilxyz[k][0];
        ybin2 = ybin + stencilxyz[k][1];
        zbin2 = zbin + stencilxyz[k][2];
        if (xbin2 < 0 || xbin2 >= mbinx ||
            ybin2 < 0 || ybin2 >= mbiny ||
            zbin2 < 0 || zbin2 >= mbinz) continue;
        for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
          if (HALF) {
            if (j <= i) continue;
          } else {
            if (i == j) continue;
          }

          if (!uniform_cut) {
            jtype = type[j];
            if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;
          }

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          // the ghost cutoff is often much shorter than the stencil, so few candidates
          // are accepted and the branch is cheaper than an unconditional store
          const double cutsq_ij = uniform_cut ? cutsq_u : cutneighghostsq_i[jtype];
          if (rsq <= cutsq_ij) {
            if (n < oneatom) neighptr[n] = j;
            n++;
          }
        }
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status()) {
      overflow = 1;
      break;
    }
  }
  NPAIR_OMP_CLOSE;

  if (overflow > 0)
    error->one(FLERR, Error::NOLASTLINE,
               "Neighbor list overflow, boost neigh_modify one" + utils::errorurl(36));

  list->inum = nlocal;
  list->gnum = nall - nlocal;
}

namespace LAMMPS_NS {
template class NPairBinGhostOmp<0>;
template class NPairBinGhostOmp<1>;
}
