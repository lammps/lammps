/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "npair_skip_intel.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_request.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairSkipIntel::NPairSkipIntel(LAMMPS *lmp) : NPair(lmp) {
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  _fix = static_cast<FixIntel *>(modify->fix[ifix]);
  _inum_starts = new int[comm->nthreads];
  _inum_counts = new int[comm->nthreads];
  _full_props = 0;
}

/* ---------------------------------------------------------------------- */

NPairSkipIntel::~NPairSkipIntel() {
  delete []_inum_starts;
  delete []_inum_counts;
  if (_full_props) delete []_full_props;
}

/* ---------------------------------------------------------------------- */

void NPairSkipIntel::copy_neighbor_info()
{
  NPair::copy_neighbor_info();
  if (_full_props) delete []_full_props;
  _full_props = new int[neighbor->nlist];
  for (int i = 0; i < neighbor->nlist; i++)
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
  int ** _noalias const firstneigh = (int ** const)list->firstneigh;
  const int * _noalias const ilist_skip = list->listskip->ilist;
  const int * _noalias const numneigh_skip = list->listskip->numneigh;
  const int ** _noalias const firstneigh_skip =
    (const int ** const)list->listskip->firstneigh;
  const int * _noalias const iskip = list->iskip;
  const int ** _noalias const ijskip = (const int ** const)list->ijskip;

  int num_skip = list->listskip->inum;
  if (list->ghost) num_skip += list->listskip->gnum;

  int packthreads;
  if (comm->nthreads > INTEL_HTHREADS && THREE==0)
    packthreads = comm->nthreads;
  else
    packthreads = 1;

  #if defined(_OPENMP)
  #pragma omp parallel if(packthreads > 1)
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
      build_t<float,0>(list, 0, 0, 0);
    else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      build_t<double,0>(list, 0, 0, 0);
    else
      build_t<float,0>(list, 0, 0, 0);
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
