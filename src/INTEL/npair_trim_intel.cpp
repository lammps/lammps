// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "npair_trim_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairTrimIntel::NPairTrimIntel(LAMMPS *lmp) : NPair(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
}

/* ----------------------------------------------------------------------
   trim from copy list to shorter cutoff
------------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void NPairTrimIntel::build_t(NeighList *list,
                                       IntelBuffers<flt_t,acc_t> *buffers)
{
  const int inum_copy = list->listcopy->inum;
  const int nlocal = atom->nlocal;
  const int e_nall = nlocal + atom->nghost;
  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const ilist = list->ilist;
  int * _noalias const numneigh = list->numneigh;
  int ** _noalias const firstneigh = list->firstneigh;
  const int * _noalias const ilist_copy = list->listcopy->ilist;
  const int * _noalias const numneigh_copy = list->listcopy->numneigh;
  const int ** _noalias const firstneigh_copy = (const int ** const)list->listcopy->firstneigh;  // NOLINT

  const flt_t cutsq_custom = cutoff_custom * cutoff_custom;

  #if defined(_OPENMP)
  #pragma omp parallel
  #endif
  {
    int tid, ifrom, ito;
    IP_PRE_omp_range_id(ifrom, ito, tid, inum_copy, comm->nthreads);

    // each thread has its own page allocator
    MyPage<int> &ipage = list->ipage[tid];
    ipage.reset();

    // loop over parent copy list
    for (int ii = ifrom; ii < ito; ii++) {
      int n = 0;
      int *neighptr = ipage.vget();

      const int i = ilist_copy[ii];
      const flt_t xtmp = x[i].x;
      const flt_t ytmp = x[i].y;
      const flt_t ztmp = x[i].z;

      // loop over copy neighbor list

      const int * _noalias const jlist = firstneigh_copy[i];
      const int jnum = numneigh_copy[i];

      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma ivdep
      #endif
      for (int jj = 0; jj < jnum; jj++) {
        const int joriginal = jlist[jj];
        const int j = joriginal & NEIGHMASK;
        int addme = 1;

        // trim to shorter cutoff

        const flt_t delx = xtmp - x[j].x;
        const flt_t dely = ytmp - x[j].y;
        const flt_t delz = ztmp - x[j].z;
        const flt_t rsq = delx * delx + dely * dely + delz * delz;

        if (rsq > cutsq_custom) addme = 0;

        if (addme)
          neighptr[n++] = joriginal;
      }

      ilist[ii] = i;
      firstneigh[i] = neighptr;
      numneigh[i] = n;

      int pad_end = n;
      IP_PRE_neighbor_pad(pad_end, 0);
      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma loop_count min=1, max=INTEL_COMPILE_WIDTH-1, \
              avg=INTEL_COMPILE_WIDTH/2
      #endif
      for ( ; n < pad_end; n++)
        neighptr[n] = e_nall;

      ipage.vgot(n);
      if (ipage.status())
        error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
    }
  }
  list->inum = inum_copy;
}

/* ---------------------------------------------------------------------- */

void NPairTrimIntel::build(NeighList *list)
{
  if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
    build_t(list, _fix->get_mixed_buffers());
  else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    build_t(list, _fix->get_double_buffers());
  else
    build_t(list, _fix->get_single_buffers());
}
