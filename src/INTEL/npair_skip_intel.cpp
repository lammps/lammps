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
   Contributing author: W. Michael Brown (Intel), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "npair_skip_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairSkipIntel::NPairSkipIntel(LAMMPS *lmp) : NPair(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
  _inum_starts = new int[comm->nthreads];
  _inum_counts = new int[comm->nthreads];
  _full_props = nullptr;
}

/* ---------------------------------------------------------------------- */

NPairSkipIntel::~NPairSkipIntel() {
  delete []_inum_starts;
  delete []_inum_counts;
  delete[] _full_props;
}

/* ---------------------------------------------------------------------- */

void NPairSkipIntel::copy_neighbor_info()
{
  NPair::copy_neighbor_info();
  // Only need to set _full_props once; npair object deleted for changes
  if (_full_props) return;
  _full_props = new int[neighbor->nrequest];
  for (int i = 0; i < neighbor->nrequest; i++)
    _full_props[i] = neighbor->requests[i]->full;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   works for half and full lists
   works for owned (non-ghost) list, also for ghost list
   iskip and ijskip flag which atom types and type pairs to skip
   if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

template<class flt_t, int THREE>
void NPairSkipIntel::build_t(NeighList *list, int *numhalf, int *cnumneigh,
                             int *numhalf_skip)
{
  const int nlocal = atom->nlocal;
  const int e_nall = nlocal + atom->nghost;
  const int * _noalias const type = atom->type;
  int * _noalias const ilist = list->ilist;
  int * _noalias const numneigh = list->numneigh;
  int ** _noalias const firstneigh = (int ** const)list->firstneigh;  // NOLINT
  const int * _noalias const ilist_skip = list->listskip->ilist;
  const int * _noalias const numneigh_skip = list->listskip->numneigh;
  const int ** _noalias const firstneigh_skip = (const int ** const)list->listskip->firstneigh;  // NOLINT
  const int * _noalias const iskip = list->iskip;
  const int **  _noalias const ijskip = (const int ** const)list->ijskip;  // NOLINT

  int num_skip = list->listskip->inum;
  if (list->ghost) num_skip += list->listskip->gnum;

  int packthreads;
  if (comm->nthreads > INTEL_HTHREADS && THREE==0)
    packthreads = comm->nthreads;
  else
    packthreads = 1;

  #if defined(_OPENMP)
  #pragma omp parallel if (packthreads > 1)
  #endif
  {
    int tid, ifrom, ito;
    IP_PRE_omp_range_id(ifrom, ito, tid, num_skip, packthreads);

    // each thread has its own page allocator
    MyPage<int> &ipage = list->ipage[tid];
    ipage.reset();

    int my_inum = ifrom;
    _inum_starts[tid] = ifrom;

    // loop over parent full list
    for (int ii = ifrom; ii < ito; ii++) {
      const int i = ilist_skip[ii];
      const int itype = type[i];
      if (iskip[itype]) continue;

      int n = 0;
      int *neighptr = ipage.vget();

      // loop over parent non-skip list

      const int * _noalias const jlist = firstneigh_skip[i];
      const int jnum = numneigh_skip[i];

      if (THREE) {
        const int jnumhalf = numhalf_skip[ii];
        for (int jj = 0; jj < jnumhalf; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;
          if (!ijskip[itype][type[j]]) neighptr[n++] = joriginal;
        }
        numhalf[my_inum] = n;

        for (int jj = jnumhalf; jj < jnum; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;
          if (!ijskip[itype][type[j]]) neighptr[n++] = joriginal;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;
          if (!ijskip[itype][type[j]]) neighptr[n++] = joriginal;
        }
      }

      ilist[my_inum++] = i;
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

    int last_inum = 0, loop_end;
    _inum_counts[tid] = my_inum;
  }
  int inum = _inum_counts[0];
  for (int tid = 1; tid < packthreads; tid++) {
    for (int i = _inum_starts[tid]; i < _inum_counts[tid]; i++) {
      if (THREE) numhalf[inum] = numhalf[i];
      ilist[inum++] = ilist[i];
    }
  }
  list->inum = inum;

  if (THREE && num_skip > 0) {
    int * const list_start = firstneigh[ilist[0]];
    for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      cnumneigh[ii] = static_cast<int>(firstneigh[i] - list_start);
    }
  }
  if (list->ghost) {
    int num = 0;
    int my_inum = list->inum;
    for (int i = 0; i < my_inum; i++)
      if (ilist[i] < nlocal) num++;
      else break;
    list->inum = num;
    list->gnum = my_inum - num;
  }
}

/* ---------------------------------------------------------------------- */

void NPairSkipIntel::build(NeighList *list)
{
  if (_fix->three_body_neighbor()==0 ||
      _full_props[list->listskip->index] == 0) {
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
      build_t<float,0>(list, nullptr, nullptr, nullptr);
    else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      build_t<double,0>(list, nullptr, nullptr, nullptr);
    else
      build_t<float,0>(list, nullptr, nullptr, nullptr);
  } else {
    int *nhalf, *cnumneigh, *nhalf_skip, *u;
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED) {
      _fix->get_mixed_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_mixed_buffers()->grow_data3(list, nhalf, cnumneigh);
      build_t<float,1>(list, nhalf, cnumneigh, nhalf_skip);
    } else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      _fix->get_double_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_double_buffers()->grow_data3(list, nhalf, cnumneigh);
      build_t<double,1>(list, nhalf, cnumneigh, nhalf_skip);
    } else {
      _fix->get_single_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_single_buffers()->grow_data3(list,nhalf,cnumneigh);
      build_t<float,1>(list, nhalf, cnumneigh, nhalf_skip);
    }
  }
}

/* ---------------------------------------------------------------------- */

NPairSkipTrimIntel::NPairSkipTrimIntel(LAMMPS *lmp) : NPair(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
  _inum_starts = new int[comm->nthreads];
  _inum_counts = new int[comm->nthreads];
  _full_props = nullptr;
}

/* ---------------------------------------------------------------------- */

NPairSkipTrimIntel::~NPairSkipTrimIntel() {
  delete []_inum_starts;
  delete []_inum_counts;
  delete[] _full_props;
}

/* ---------------------------------------------------------------------- */

void NPairSkipTrimIntel::copy_neighbor_info()
{
  NPair::copy_neighbor_info();
  // Only need to set _full_props once; npair object deleted for changes
  if (_full_props) return;
  _full_props = new int[neighbor->nrequest];
  for (int i = 0; i < neighbor->nrequest; i++)
    _full_props[i] = neighbor->requests[i]->full;
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   works for half and full lists
   works for owned (non-ghost) list, also for ghost list
   iskip and ijskip flag which atom types and type pairs to skip
   if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int THREE>
void NPairSkipTrimIntel::build_t(NeighList *list, int *numhalf, int *cnumneigh,
                             int *numhalf_skip, IntelBuffers<flt_t,acc_t> *buffers)
{
  const int nlocal = atom->nlocal;
  const int e_nall = nlocal + atom->nghost;
  const ATOM_T * _noalias const x = buffers->get_x();
  const int * _noalias const type = atom->type;
  int * _noalias const ilist = list->ilist;
  int * _noalias const numneigh = list->numneigh;
  int ** _noalias const firstneigh = (int ** const)list->firstneigh;  // NOLINT
  const int * _noalias const ilist_skip = list->listskip->ilist;
  const int * _noalias const numneigh_skip = list->listskip->numneigh;
  const int ** _noalias const firstneigh_skip = (const int ** const)list->listskip->firstneigh;  // NOLINT
  const int * _noalias const iskip = list->iskip;
  const int **  _noalias const ijskip = (const int ** const)list->ijskip;  // NOLINT

  const flt_t cutsq_custom = cutoff_custom * cutoff_custom;
  int num_skip = list->listskip->inum;
  if (list->ghost) num_skip += list->listskip->gnum;

  int packthreads;
  if (comm->nthreads > INTEL_HTHREADS && THREE==0)
    packthreads = comm->nthreads;
  else
    packthreads = 1;

  #if defined(_OPENMP)
  #pragma omp parallel if (packthreads > 1)
  #endif
  {
    int tid, ifrom, ito;
    IP_PRE_omp_range_id(ifrom, ito, tid, num_skip, packthreads);

    // each thread has its own page allocator
    MyPage<int> &ipage = list->ipage[tid];
    ipage.reset();

    int my_inum = ifrom;
    _inum_starts[tid] = ifrom;

    // loop over parent full list
    for (int ii = ifrom; ii < ito; ii++) {
      const int i = ilist_skip[ii];
      const int itype = type[i];
      if (iskip[itype]) continue;

      const flt_t xtmp = x[i].x;
      const flt_t ytmp = x[i].y;
      const flt_t ztmp = x[i].z;

      int n = 0;
      int *neighptr = ipage.vget();

      // loop over parent non-skip list

      const int * _noalias const jlist = firstneigh_skip[i];
      const int jnum = numneigh_skip[i];

      if (THREE) {
        const int jnumhalf = numhalf_skip[ii];
        for (int jj = 0; jj < jnumhalf; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;

          int addme = 1;
          if (ijskip[itype][type[j]]) addme = 0;

          // trim to shorter cutoff

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutsq_custom) addme = 0;

          if (addme)
            neighptr[n++] = joriginal;
        }
        numhalf[my_inum] = n;

        for (int jj = jnumhalf; jj < jnum; jj++) {
          const int joriginal = jlist[jj];
          const int j = joriginal & NEIGHMASK;

          int addme = 1;
          if (ijskip[itype][type[j]]) addme = 0;

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
          if (ijskip[itype][type[j]]) addme = 0;

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

      ilist[my_inum++] = i;
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

    int last_inum = 0, loop_end;
    _inum_counts[tid] = my_inum;
  }
  int inum = _inum_counts[0];
  for (int tid = 1; tid < packthreads; tid++) {
    for (int i = _inum_starts[tid]; i < _inum_counts[tid]; i++) {
      if (THREE) numhalf[inum] = numhalf[i];
      ilist[inum++] = ilist[i];
    }
  }
  list->inum = inum;

  if (THREE && num_skip > 0) {
    int * const list_start = firstneigh[ilist[0]];
    for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      cnumneigh[ii] = static_cast<int>(firstneigh[i] - list_start);
    }
  }
  if (list->ghost) {
    int num = 0;
    int my_inum = list->inum;
    for (int i = 0; i < my_inum; i++)
      if (ilist[i] < nlocal) num++;
      else break;
    list->inum = num;
    list->gnum = my_inum - num;
  }
}

/* ---------------------------------------------------------------------- */

void NPairSkipTrimIntel::build(NeighList *list)
{
  if (_fix->three_body_neighbor()==0 ||
      _full_props[list->listskip->index] == 0) {
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
      build_t<float,double,0>(list, nullptr, nullptr, nullptr, _fix->get_mixed_buffers());
    else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      build_t<double,double,0>(list, nullptr, nullptr, nullptr, _fix->get_double_buffers());
    else
      build_t<float,float,0>(list, nullptr, nullptr, nullptr, _fix->get_single_buffers());
  } else {
    int *nhalf, *cnumneigh, *nhalf_skip, *u;
    if (_fix->precision() == FixIntel::PREC_MODE_MIXED) {
      _fix->get_mixed_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_mixed_buffers()->grow_data3(list, nhalf, cnumneigh);
      build_t<float,double,1>(list, nhalf, cnumneigh, nhalf_skip, _fix->get_mixed_buffers());
    } else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      _fix->get_double_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_double_buffers()->grow_data3(list, nhalf, cnumneigh);
      build_t<double,double,1>(list, nhalf, cnumneigh, nhalf_skip, _fix->get_double_buffers());
    } else {
      _fix->get_single_buffers()->get_list_data3(list->listskip,nhalf_skip,u);
      _fix->get_single_buffers()->grow_data3(list,nhalf,cnumneigh);
      build_t<float,float,1>(list, nhalf, cnumneigh, nhalf_skip, _fix->get_single_buffers());
    }
  }
}
