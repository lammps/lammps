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

#include "npair_halffull_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "neigh_list_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType, int NEWTON, int TRI, int TRIM>
NPairHalffullKokkos<DeviceType,NEWTON,TRI,TRIM>::NPairHalffullKokkos(LAMMPS *lmp) : NPair(lmp) {
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

template<class DeviceType, int NEWTON, int TRI, int TRIM>
void NPairHalffullKokkos<DeviceType,NEWTON,TRI,TRIM>::build(NeighList *list)
{
  if (NEWTON || TRIM) {
    x = atomKK->k_x.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|TYPE_MASK);

    NeighborKokkos* neighborKK = (NeighborKokkos*) neighbor;
    neighborKK->k_cutneighsq.template sync<DeviceType>();
    d_cutneighsq = neighborKK->k_cutneighsq.template view<DeviceType>();
  }

  nlocal = atom->nlocal;

  cutsq_custom = cutoff_custom*cutoff_custom;

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

  delta = 0.01 * force->angstrom;

  // loop over parent full list

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagNPairHalffullCompute>(0,inum_full),*this);
  copymode = 0;

  list->inum = k_list_full->inum;
  list->gnum = k_list_full->gnum;

  k_list->k_ilist.template modify<DeviceType>();
}

template<class DeviceType, int NEWTON, int TRI, int TRIM>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void NPairHalffullKokkos<DeviceType,NEWTON,TRI,TRIM>::operator()(TagNPairHalffullCompute, const int &ii) const {
  int n = 0;

  const int i = d_ilist_full(ii);
  double xtmp,ytmp,ztmp;
  if (NEWTON || TRIM) {
    xtmp = static_cast<double>(x(i,0));
    ytmp = static_cast<double>(x(i,1));
    ztmp = static_cast<double>(x(i,2));
  }

  // loop over full neighbor list
  // use i < j < nlocal to eliminate half the local/local interactions
  // for triclinic, must use delta to eliminate half the local/ghost interactions
  // cannot use I/J exact coord comparision as for orthog
  //   b/c transforming orthog -> lambda -> orthog for ghost atoms
  //   with an added PBC offset can shift all 3 coords by epsilon

  const int jnum = d_numneigh_full(i);
  const AtomNeighbors neighbors_i = AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                                                  &d_neighbors(i,1)-&d_neighbors(i,0));

  for (int jj = 0; jj < jnum; jj++) {
    const int joriginal = d_neighbors_full(i,jj);
    const int j = joriginal & NEIGHMASK;
    if (NEWTON) {
      if (j < nlocal) {
        if (i > j) continue;
      } else if (TRI) {
        if (fabs(static_cast<double>(x(j,2))-ztmp) > delta) {
          if (static_cast<double>(x(j,2)) < ztmp) continue;
        } else if (fabs(static_cast<double>(x(j,1))-ytmp) > delta) {
          if (static_cast<double>(x(j,1)) < ytmp) continue;
        } else {
          if (static_cast<double>(x(j,0)) < xtmp) continue;
        }
      } else {
        if (static_cast<double>(x(j,2)) < ztmp) continue;
        if (static_cast<double>(x(j,2)) == ztmp) {
          if (static_cast<double>(x(j,1)) < ytmp) continue;
          if (static_cast<double>(x(j,1)) == ytmp && static_cast<double>(x(j,0)) < xtmp) continue;
        }
      }

      if (TRIM) {
        const double delx = xtmp - static_cast<double>(x(j,0));
        const double dely = ytmp - static_cast<double>(x(j,1));
        const double delz = ztmp - static_cast<double>(x(j,2));
        const double rsq = delx*delx + dely*dely + delz*delz;

        // a trim list whose own request carries no custom cutoff must fall back
        // to the pairwise neighbour cutoff, as NPairHalffull::build() does
        const double cutsq_trim = (cutsq_custom > 0.0) ? cutsq_custom :
          static_cast<double>(d_cutneighsq(type(i),type(j)));
        if (rsq > cutsq_trim) continue;
      }

      neighbors_i(n++) = joriginal;
    } else if (j > i) {

      if (TRIM) {
        const double delx = xtmp - static_cast<double>(x(j,0));
        const double dely = ytmp - static_cast<double>(x(j,1));
        const double delz = ztmp - static_cast<double>(x(j,2));
        const double rsq = delx*delx + dely*dely + delz*delz;

        // a trim list whose own request carries no custom cutoff must fall back
        // to the pairwise neighbour cutoff, as NPairHalffull::build() does
        const double cutsq_trim = (cutsq_custom > 0.0) ? cutsq_custom :
          static_cast<double>(d_cutneighsq(type(i),type(j)));
        if (rsq > cutsq_trim) continue;
      }

      neighbors_i(n++) = joriginal;
    }
  }

  d_numneigh(i) = n;
  d_ilist(ii) = i;
}

namespace LAMMPS_NS {
template class NPairHalffullKokkos<LMPDeviceType,0,0,0>;
template class NPairHalffullKokkos<LMPDeviceType,0,0,1>;
template class NPairHalffullKokkos<LMPDeviceType,1,0,0>;
template class NPairHalffullKokkos<LMPDeviceType,1,0,1>;
template class NPairHalffullKokkos<LMPDeviceType,1,1,0>;
template class NPairHalffullKokkos<LMPDeviceType,1,1,1>;
#ifdef LMP_KOKKOS_GPU
template class NPairHalffullKokkos<LMPHostType,0,0,0>;
template class NPairHalffullKokkos<LMPHostType,0,0,1>;
template class NPairHalffullKokkos<LMPHostType,1,0,0>;
template class NPairHalffullKokkos<LMPHostType,1,0,1>;
template class NPairHalffullKokkos<LMPHostType,1,1,0>;
template class NPairHalffullKokkos<LMPHostType,1,1,1>;
#endif
}
