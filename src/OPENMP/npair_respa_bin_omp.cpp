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

#include "npair_respa_bin_omp.h"
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

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int NEWTON, int TRI>
NPairRespaBinOmp<NEWTON, TRI>::NPairRespaBinOmp(LAMMPS *lmp) :
  NPair(lmp) {}

/* ----------------------------------------------------------------------
   multiple respa lists
   Newtoff
     binned neighbor list construction with partial Newton's 3rd law
     each owned atom i checks own bin and surrounding bins in non-Newton stencil
     pair stored once if i,j are both owned and i < j
     pair stored by me if j is ghost (also stored by proc owning j)
  Newton
     binned neighbor list construction with full Newton's 3rd law
     each owned atom i checks its own bin and other bins in Newton stencil
     every pair stored exactly once by some processor
------------------------------------------------------------------------- */

template<int NEWTON, int TRI>
void NPairRespaBinOmp<NEWTON, TRI>::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  const double delta = 0.01 * force->angstrom;

  NPAIR_OMP_INIT;

  const int respamiddle = list->respamiddle;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i, j, k, n, itype, jtype, ibin, bin_start, n_inner, n_middle, imol, iatom;
  tagint itag, jtag, tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *neighptr, *neighptr_inner, *neighptr_middle;

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

  int *ilist_inner = list->ilist_inner;
  int *numneigh_inner = list->numneigh_inner;
  int **firstneigh_inner = list->firstneigh_inner;

  int *ilist_middle, *numneigh_middle, **firstneigh_middle;
  if (respamiddle) {
    ilist_middle = list->ilist_middle;
    numneigh_middle = list->numneigh_middle;
    firstneigh_middle = list->firstneigh_middle;
  }

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  MyPage<int> &ipage_inner = list->ipage_inner[tid];
  ipage.reset();
  ipage_inner.reset();

  MyPage<int> *ipage_middle;
  if (respamiddle) {
    ipage_middle = list->ipage_middle + tid;
    ipage_middle->reset();
  }

  int which = 0;
  int minchange = 0;

  for (i = ifrom; i < ito; i++) {

    n = n_inner = 0;
    neighptr = ipage.vget();
    neighptr_inner = ipage_inner.vget();
    if (respamiddle) {
      n_middle = 0;
      neighptr_middle = ipage_middle->vget();
    }

    itag = tag[i];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    ibin = atom2bin[i];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    for (k = 0; k < nstencil; k++) {
      bin_start = binhead[ibin + stencil[k]];
      if (NEWTON && (!TRI)) {
        if (k == 0) {
          // Half neighbor list, newton on, orthonormal
          // loop over rest of atoms in i's bin, ghosts are at end of linked list
          bin_start = bins[i];
        }
      }

      for (j = bin_start; j >= 0; j = bins[j]) {
        if (!NEWTON) {
          // Half neighbor list, newton off
          // only store pair if i < j
          // stores own/own pairs only once
          // stores own/ghost pairs on both procs
          if (j <= i) continue;
        } else if (TRI) {
          // Half neighbor list, newton on, triclinic
          // for triclinic, bin stencil is full in all 3 dims
          // must use itag/jtag to eliminate half the I/J interactions
          // cannot use I/J exact coord comparision
          //   b/c transforming orthog -> lambda -> orthog for ghost atoms
          //   with an added PBC offset can shift all 3 coords by epsilon
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
        } else {
          // Half neighbor list, newton on, orthonormal
          // store every pair for every bin in stencil,except for i's bin

          if (k == 0) {
            // if j is owned atom, store it, since j is beyond i in linked list
            // if j is ghost, only store if j coords are "above and to the "right" of i
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
            else if ((minchange = domain->minimum_image_check(delx, dely, delz)))
              neighptr[n++] = j;
            else if (which > 0)
              neighptr[n++] = j ^ (which << SBBITS);
          } else
            neighptr[n++] = j;

          if (rsq < cut_inner_sq) {
            if (which == 0)
              neighptr_inner[n_inner++] = j;
            else if (minchange)
              neighptr_inner[n_inner++] = j;
            else if (which > 0)
              neighptr_inner[n_inner++] = j ^ (which << SBBITS);
          }

          if (respamiddle &&
              rsq < cut_middle_sq && rsq > cut_middle_inside_sq) {
            if (which == 0)
              neighptr_middle[n_middle++] = j;
            else if (minchange)
              neighptr_middle[n_middle++] = j;
            else if (which > 0)
              neighptr_middle[n_middle++] = j ^ (which << SBBITS);
          }
        }
      }
    }

    ilist[i] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage.vgot(n);
    if (ipage.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

    ilist_inner[i] = i;
    firstneigh_inner[i] = neighptr_inner;
    numneigh_inner[i] = n_inner;
    ipage.vgot(n_inner);
    if (ipage_inner.status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");

    if (respamiddle) {
      ilist_middle[i] = i;
      firstneigh_middle[i] = neighptr_middle;
      numneigh_middle[i] = n_middle;
      ipage_middle->vgot(n_middle);
      if (ipage_middle->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
    }
  }
  NPAIR_OMP_CLOSE;
  list->inum = nlocal;
  list->inum_inner = nlocal;
  if (respamiddle) list->inum_middle = nlocal;
}

namespace LAMMPS_NS {
template class NPairRespaBinOmp<0,0>;
template class NPairRespaBinOmp<1,0>;
template class NPairRespaBinOmp<1,1>;
}
