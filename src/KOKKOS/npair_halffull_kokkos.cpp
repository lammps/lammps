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

#include "npair_halffull_kokkos.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType, int NEWTON>
NPairHalffullKokkos<DeviceType,NEWTON>::NPairHalffullKokkos(LAMMPS *lmp) : NPair(lmp) {
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ----------------------------------------------------------------------
   build half list from full list
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
   works if full list is a skip list
   works for owned (non-ghost) list, also for ghost list
   if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

template<class DeviceType, int NEWTON>
void NPairHalffullKokkos<DeviceType,NEWTON>::build(NeighList *list)
{
  if (NEWTON) {
    x = atomKK->k_x.view<DeviceType>();
    atomKK->sync(execution_space,X_MASK);
  }
  nlocal = atom->nlocal;

  NeighListKokkos<DeviceType>* k_list_full = static_cast<NeighListKokkos<DeviceType>*>(list->listfull);
  d_ilist_full = k_list_full->d_ilist;
  d_numneigh_full = k_list_full->d_numneigh;
  d_neighbors_full = k_list_full->d_neighbors;
  int inum_full = list->listfull->inum;
  if (list->ghost) inum_full += list->listfull->gnum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  k_list->maxneighs = k_list_full->maxneighs; // simple, but could be made more memory efficient
  k_list->grow(atom->nmax);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  // loop over parent full list

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagNPairHalffullCompute>(0,inum_full),*this);

  list->inum = k_list_full->inum;
  list->gnum = k_list_full->gnum;
  k_list->k_ilist.template modify<DeviceType>();

  copymode = 0;
}

template<class DeviceType, int NEWTON>
KOKKOS_INLINE_FUNCTION
void NPairHalffullKokkos<DeviceType,NEWTON>::operator()(TagNPairHalffullCompute, const int &ii) const {
  int n = 0;

  const int i = d_ilist_full(ii);
  F_FLOAT xtmp,ytmp,ztmp;
  if (NEWTON) {
    xtmp = x(i,0);
    ytmp = x(i,1);
    ztmp = x(i,2);
  }

  // loop over full neighbor list

  const int jnum = d_numneigh_full(i);
  const AtomNeighbors neighbors_i = AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                                                  &d_neighbors(i,1)-&d_neighbors(i,0));

  for (int jj = 0; jj < jnum; jj++) {
    const int joriginal = d_neighbors_full(i,jj);
    const int j = joriginal & NEIGHMASK;
    if (NEWTON) {
      if (j < nlocal) {
        if (i > j) continue;
      } else {
        if (x(j,2) < ztmp) continue;
        if (x(j,2) == ztmp) {
          if (x(j,1) < ytmp) continue;
          if (x(j,1) == ytmp && x(j,0) < xtmp) continue;
        }
      }
      neighbors_i(n++) = joriginal;
    } else {
      if (j > i) neighbors_i(n++) = joriginal;
    }
  }

  d_numneigh(i) = n;
  d_ilist(ii) = i;
}

namespace LAMMPS_NS {
template class NPairHalffullKokkos<LMPDeviceType,0>;
template class NPairHalffullKokkos<LMPDeviceType,1>;
#ifdef KOKKOS_ENABLE_CUDA
template class NPairHalffullKokkos<LMPHostType,0>;
template class NPairHalffullKokkos<LMPHostType,1>;
#endif
}
