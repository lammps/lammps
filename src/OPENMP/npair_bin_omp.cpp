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

#include "npair_bin_omp.h"
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

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
NPairBinOmp<HALF, NEWTON, TRI, SIZE, ATOMONLY>::NPairBinOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
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

template<int HALF, int NEWTON, int TRI, int SIZE, int ATOMONLY>
void NPairBinOmp<HALF, NEWTON, TRI, SIZE, ATOMONLY>::build(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  const int molecular = atom->molecular;
  const int moltemplate = (molecular == Atom::TEMPLATE) ? 1 : 0;
  const double delta = 0.01 * force->angstrom;
  int overflow = 0;

  NPAIR_OMP_INIT;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list,fix_package_omp) reduction(+:overflow)
#endif
  NPAIR_OMP_SETUP(nlocal);

  int i, j, jh, k, n, itype, jtype, ibin, bin_start, which, imol, iatom;
  tagint itag, jtag, tagprev;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, radsum, cut, cutsq;
  int * _noalias neighptr;

  const double (* _noalias const x)[3] = (const double (*)[3]) atom->x[0];
  double *radius = atom->radius;
  const int * _noalias const type = atom->type;
  int *mask = atom->mask;
  const tagint * _noalias const tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int history = list->history;
  int mask_history = 1 << HISTBITS;
  const int oneatom = list->oneatom;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  const int * _noalias const bins = this->bins;
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const stencil = this->stencil;
  const int * _noalias const atom2bin = this->atom2bin;
  const int nstencil = this->nstencil;

  // cache whether type t has identical neighbor cutoffs to all types
  // covers models up to 63 atom types on stack; larger counts fall back safely
  const int ntypes = atom->ntypes;
  bool ucut_type[64];
  if (!SIZE && !exclude && ntypes < 64) {
    for (int t = 1; t <= ntypes; t++) {
      bool u = true;
      for (int t2 = 2; u && t2 <= ntypes; t2++)
        if (cutneighsq[t][t2] != cutneighsq[t][1]) u = false;
      ucut_type[t] = u;
    }
  }

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
    if (!ATOMONLY) {
      if (moltemplate) {
        imol = molindex[i];
        iatom = molatom[i];
        tagprev = tag[i] - iatom - 1;
      }
    }

    // loop over all atoms in surrounding bins in stencil including self
    // skip i = j

    ibin = atom2bin[i];
    const double * _noalias const cutneighsq_i = cutneighsq[itype];
    const double cutsq_u = cutneighsq_i[1];
    const bool uniform_cut = (!SIZE && !exclude && ntypes < 64) ? ucut_type[itype] : false;

    for (k = 0; k < nstencil; k++) {
      bin_start = binhead[ibin + stencil[k]];
      if (HALF && NEWTON && (!TRI)) {
        if (k == 0) {
          // Half neighbor list, newton on, orthonormal
          // loop over rest of atoms in i's bin, ghosts are at end of linked list
          bin_start = bins[i];
        }
      }

      for (j = bin_start; j >= 0; j = bins[j]) {
        if (!HALF) {
          // Full neighbor list
          // only skip i = j
          if (i == j) continue;
        } else if (!NEWTON) {
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
          // store every pair for every bin in stencil, except for i's bin

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

        // skip type[j] load when !SIZE && !exclude and all cutoffs from itype match
        if (!uniform_cut) {
          jtype = type[j];
          if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;
        }

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        if (SIZE) {
          radsum = radius[i] + radius[j];
          cut = radsum + skin;
          cutsq = cut * cut;

          if (ATOMONLY) {
            if (rsq <= cutsq) {
              jh = j;
              if (history && rsq < radsum * radsum)
                jh = jh ^ mask_history;
              neighptr[n++] = jh;
            }
          } else {
            if (rsq <= cutsq) {
              jh = j;
              if (history && rsq < radsum * radsum)
                jh = jh ^ mask_history;

              if (molecular != Atom::ATOMIC) {
                if (!moltemplate)
                  which = find_special(special[i], nspecial[i], tag[j]);
                else if ((imol >= 0) && onemols[imol]->special)
                  which = find_special(onemols[imol]->special[iatom],
                                       onemols[imol]->nspecial[iatom], tag[j] - tagprev);
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
          }
        } else {
          const double cutsq_ij = uniform_cut ? cutsq_u : cutneighsq_i[jtype];
          if (ATOMONLY) {
            if (n < oneatom) {
              neighptr[n] = j;
              n += (rsq <= cutsq_ij);
            } else if (rsq <= cutsq_ij) n++;
          } else {
            if (molecular != Atom::ATOMIC) {
              if (rsq <= cutsq_ij) {
                if (!moltemplate)
                  which = find_special(special[i], nspecial[i], tag[j]);
                else if ((imol >= 0) && onemols[imol]->special)
                  which = find_special(onemols[imol]->special[iatom],
                                       onemols[imol]->nspecial[iatom], tag[j] - tagprev);
                else which = 0;
                if (which == 0)
                  neighptr[n++] = j;
                else if (domain->minimum_image_check(delx, dely, delz))
                  neighptr[n++] = j;
                else if (which > 0)
                  neighptr[n++] = j ^ (which << SBBITS);
              }
            } else {
              if (n < oneatom) {
                neighptr[n] = j;
                n += (rsq <= cutsq_ij);
              } else if (rsq <= cutsq_ij) n++;
            }
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
  if (!HALF) list->gnum = 0;
}

namespace LAMMPS_NS {
template class NPairBinOmp<0,1,0,0,0>;
template class NPairBinOmp<1,0,0,0,0>;
template class NPairBinOmp<1,1,0,0,0>;
template class NPairBinOmp<1,1,1,0,0>;
template class NPairBinOmp<0,1,0,1,0>;
template class NPairBinOmp<1,0,0,1,0>;
template class NPairBinOmp<1,1,0,1,0>;
template class NPairBinOmp<1,1,1,1,0>;
template class NPairBinOmp<0,1,0,0,1>;
template class NPairBinOmp<1,0,0,0,1>;
template class NPairBinOmp<1,1,0,0,1>;
template class NPairBinOmp<1,1,1,0,1>;
template class NPairBinOmp<0,1,0,1,1>;
template class NPairBinOmp<1,0,0,1,1>;
template class NPairBinOmp<1,1,0,1,1>;
template class NPairBinOmp<1,1,1,1,1>;
}
