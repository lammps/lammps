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

#include "npair_skip_kokkos.h"
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

template<class DeviceType>
NPairSkipKokkos<DeviceType>::NPairSkipKokkos(LAMMPS *lmp) : NPair(lmp) {
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  d_inum = typename AT::t_int_scalar("npair_skip:inum");
}

/* ----------------------------------------------------------------------
   build skip list for subset of types from parent list
   works for half and full lists
   works for owned (non-ghost) list, also for ghost list
   iskip and ijskip flag which atom types and type pairs to skip
   if ghost, also store neighbors of ghost atoms & set inum,gnum correctly
------------------------------------------------------------------------- */

template<class DeviceType>
void NPairSkipKokkos<DeviceType>::build(NeighList *list)
{
  atomKK->sync(execution_space,TYPE_MASK);
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;


  NeighListKokkos<DeviceType>* k_list_skip = static_cast<NeighListKokkos<DeviceType>*>(list->listskip);
  d_ilist_skip = k_list_skip->d_ilist;
  d_numneigh_skip = k_list_skip->d_numneigh;
  d_neighbors_skip = k_list_skip->d_neighbors;

  num_skip = list->listskip->inum;
  if (list->ghost) num_skip += list->listskip->gnum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  k_list->maxneighs = k_list_skip->maxneighs; // simple, but could be made more memory efficient
  k_list->grow(atom->nmax);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  int ntypes = atom->ntypes;

  k_iskip = DAT::tdual_int_1d("npair_skip:iskip",ntypes+1);
  k_ijskip = DAT::tdual_int_2d("npair_skip:ijskip",ntypes+1,ntypes+1);
  d_iskip = k_iskip.view<DeviceType>();
  d_ijskip = k_ijskip.view<DeviceType>();

  for (int itype = 1; itype <= ntypes; itype++) {
    k_iskip.h_view(itype) = list->iskip[itype];
    for (int jtype = 1; jtype <= ntypes; jtype++) {
      k_ijskip.h_view(itype,jtype) = list->ijskip[itype][jtype];
    }
  }
  k_iskip.modify<LMPHostType>();
  k_ijskip.modify<LMPHostType>();

  k_iskip.sync<DeviceType>();
  k_ijskip.sync<DeviceType>();

  // loop over atoms in other list
  // skip I atom entirely if iskip is set for type[I]
  // skip I,J pair if ijskip is set for type[I],type[J]

  copymode = 1;
  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType, TagNPairSkipCompute>(0,num_skip),*this);

  auto h_inum = Kokkos::create_mirror_view(d_inum);
  Kokkos::deep_copy(h_inum,d_inum);
  const int inum = h_inum();
  list->inum = inum;
  if (list->ghost) {
    int num = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagNPairSkipCountLocal>(0,inum),*this,num);
    list->inum = num;
    list->gnum = inum - num;
  }
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NPairSkipKokkos<DeviceType>::operator()(TagNPairSkipCompute, const int &ii, int &inum, const bool &final) const {

  const int i = d_ilist_skip(ii);
  const int itype = type(i);

  if (!d_iskip(itype)) {

    if (final) {

      int n = 0;

      // loop over parent non-skip list

      const int jnum = d_numneigh_skip(i);
      const AtomNeighbors neighbors_i = AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                                                    &d_neighbors(i,1)-&d_neighbors(i,0));

      for (int jj = 0; jj < jnum; jj++) {
        const int joriginal = d_neighbors_skip(i,jj);
        int j = joriginal & NEIGHMASK;
        if (d_ijskip(itype,type(j))) continue;
        neighbors_i(n++) = joriginal;
      }

      d_numneigh(i) = n;
      d_ilist(inum) = i;
    }

    inum++;
  }

  if (final) {
    if (ii == num_skip-1)
      d_inum() = inum;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NPairSkipKokkos<DeviceType>::operator()(TagNPairSkipCountLocal, const int &i, int &num) const {
  if (d_ilist[i] < nlocal) num++;
}


namespace LAMMPS_NS {
template class NPairSkipKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class NPairSkipKokkos<LMPHostType>;
#endif
}
