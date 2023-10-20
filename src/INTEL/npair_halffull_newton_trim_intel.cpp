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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "npair_halffull_newton_trim_intel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "my_page.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalffullNewtonTrimIntel::NPairHalffullNewtonTrimIntel(LAMMPS *lmp) : NPair(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
}

/* ----------------------------------------------------------------------
   build half list from full list and trim to shorter cutoff
   pair stored once if i,j are both owned and i < j
   if j is ghost, only store if j coords are "above and to the right" of i
   works if full list is a skip list
------------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void NPairHalffullNewtonTrimIntel::build_t(NeighList *list,
                                       IntelBuffers<flt_t,acc_t> *buffers)
{
  const int inum_full = list->listfull->inum;
  const int nlocal = atom->nlocal;
  const int e_nall = nlocal + atom->nghost;
  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const ilist = list->ilist;
  int * _noalias const numneigh = list->numneigh;
  int ** _noalias const firstneigh = list->firstneigh;
  const int * _noalias const ilist_full = list->listfull->ilist;
  const int * _noalias const numneigh_full = list->listfull->numneigh;
  const int ** _noalias const firstneigh_full = (const int ** const)list->listfull->firstneigh;  // NOLINT

  const flt_t cutsq_custom = cutoff_custom * cutoff_custom;
  const double delta = 0.01 * force->angstrom;
  const int triclinic = domain->triclinic;

  #if defined(_OPENMP)
  #pragma omp parallel
  #endif
  {
    int tid, ifrom, ito;
    IP_PRE_omp_range_id(ifrom, ito, tid, inum_full, comm->nthreads);

    // each thread has its own page allocator
    MyPage<int> &ipage = list->ipage[tid];
    ipage.reset();

    // loop over parent full list
    for (int ii = ifrom; ii < ito; ii++) {
      int n = 0;
      int *neighptr = ipage.vget();

      const int i = ilist_full[ii];
      const flt_t xtmp = x[i].x;
      const flt_t ytmp = x[i].y;
      const flt_t ztmp = x[i].z;

      // loop over full neighbor list

      const int * _noalias const jlist = firstneigh_full[i];
      const int jnum = numneigh_full[i];

      if (!triclinic) {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;
          int addme = 1;
          if (j < nlocal) {
            if (i > j) addme = 0;
          } else {
            if (x[j].z < ztmp) addme = 0;
            if (x[j].z == ztmp) {
              if (x[j].y < ytmp) addme = 0;
              if (x[j].y == ytmp && x[j].x < xtmp) addme = 0;
            }
          }

          // trim to shorter cutoff

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          if (rsq > cutsq_custom) addme = 0;

          if (addme)
            neighptr[n++] = joriginal;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;
          int addme = 1;
          if (j < nlocal) {
            if (i > j) addme = 0;
          } else {
            if (fabs(x[j].z-ztmp) > delta) {
              if (x[j].z < ztmp) addme = 0;
            } else if (fabs(x[j].y-ytmp) > delta) {
              if (x[j].y < ytmp) addme = 0;
            } else {
              if (x[j].x < xtmp) addme = 0;
            }
          }

          // trim to shorter cutoff

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          if (rsq > cutsq_custom) addme = 0;

          if (addme)
            neighptr[n++] = joriginal;
        }
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
  list->inum = inum_full;
}

/* ----------------------------------------------------------------------
   build half list from full 3-body list and trim to shorter cutoff
   half list is already stored as first part of 3-body list
------------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void NPairHalffullNewtonTrimIntel::build_t3(NeighList *list, int *numhalf,
                                            IntelBuffers<flt_t,acc_t> *buffers)
{
  const int inum_full = list->listfull->inum;
  const int e_nall = atom->nlocal + atom->nghost;
  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const ilist = list->ilist;
  int * _noalias const numneigh = list->numneigh;
  int ** _noalias const firstneigh = list->firstneigh;
  const int * _noalias const ilist_full = list->listfull->ilist;
  const int * _noalias const numneigh_full = numhalf;
  const int ** _noalias const firstneigh_full = (const int ** const)list->listfull->firstneigh;  // NOLINT

  const flt_t cutsq_custom = cutoff_custom * cutoff_custom;

  int packthreads = 1;
  if (comm->nthreads > INTEL_HTHREADS) packthreads = comm->nthreads;

  #if defined(_OPENMP)
  #pragma omp parallel if (packthreads > 1)
  #endif
  {
    int tid, ifrom, ito;
    IP_PRE_omp_range_id(ifrom, ito, tid, inum_full, packthreads);

    // each thread has its own page allocator
    MyPage<int> &ipage = list->ipage[tid];
    ipage.reset();

    // loop over parent full list
    for (int ii = ifrom; ii < ito; ii++) {
      int n = 0;
      int *neighptr = ipage.vget();

      const int i = ilist_full[ii];
      const flt_t xtmp = x[i].x;
      const flt_t ytmp = x[i].y;
      const flt_t ztmp = x[i].z;

      // loop over full neighbor list

      const int * _noalias const jlist = firstneigh_full[i];
      const int jnum = numneigh_full[ii];

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
  list->inum = inum_full;
}

/* ---------------------------------------------------------------------- */

void NPairHalffullNewtonTrimIntel::build(NeighList *list)
{
  if (_fix->three_body_neighbor() == 0 || domain->triclinic) {
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
      build_t(list, _fix->get_mixed_buffers());
    else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      build_t(list, _fix->get_double_buffers());
    else
      build_t(list, _fix->get_single_buffers());
  } else {
    int *nhalf, *cnum;
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED) {
      _fix->get_mixed_buffers()->get_list_data3(list->listfull, nhalf, cnum);
      build_t3<float>(list, nhalf, _fix->get_mixed_buffers());
    } else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      _fix->get_double_buffers()->get_list_data3(list->listfull, nhalf, cnum);
      build_t3<double>(list, nhalf, _fix->get_double_buffers());
    } else {
      _fix->get_single_buffers()->get_list_data3(list->listfull, nhalf, cnum);
      build_t3<float>(list, nhalf, _fix->get_single_buffers());
    }
  }
}
