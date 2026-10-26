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

#include "npair_trim_omp.h"
#include "npair_omp.h"
#include "omp_compat.h"

#include "atom.h"
#include "error.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairTrimOmp::NPairTrimOmp(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   trim from copy list to shorter cutoff
   if ghost, also trim the neighbors of ghost atoms and set inum,gnum
------------------------------------------------------------------------- */

void NPairTrimOmp::build(NeighList *list)
{
  const int inum_copy = list->listcopy->inum;
  const int gnum_copy = list->listcopy->gnum;
  const int inum_trim = list->ghost ? inum_copy + gnum_copy : inum_copy;
  int overflow = 0;

  NPAIR_OMP_INIT;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(list,fix_package_omp) reduction(+:overflow)
#endif
  NPAIR_OMP_SETUP(inum_trim);

  int i, j, ii, jj, n, jnum, joriginal, itype, jtype;
  int *neighptr, *jlist;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

  double **x = atom->x;
  int *type = atom->type;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *ilist_copy = list->listcopy->ilist;
  int *numneigh_copy = list->listcopy->numneigh;
  int **firstneigh_copy = list->listcopy->firstneigh;

  // each thread has its own page allocator
  MyPage<int> &ipage = list->ipage[tid];
  ipage.reset();

  double cutsq_custom = cutoff_custom * cutoff_custom;

  // loop over atoms in copy list

  for (ii = ifrom; ii < ito; ii++) {

    n = 0;
    neighptr = ipage.vget();

    // loop over parent copy list

    i = ilist_copy[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh_copy[i];
    jnum = numneigh_copy[i];

    for (jj = 0; jj < jnum; jj++) {
      joriginal = jlist[jj];
      j = joriginal & NEIGHMASK;
      jtype = type[j];

      // trim to shorter cutoff

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      double cutsq_trim = (cutsq_custom > 0.0) ? cutsq_custom : cutneighsq[itype][jtype];
      if (rsq > cutsq_trim) continue;

      neighptr[n++] = joriginal;
    }

    ilist[ii] = i;
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

  list->inum = inum_copy;
  list->gnum = gnum_copy;
}
