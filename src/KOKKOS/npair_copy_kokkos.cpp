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

#include "npair_copy_kokkos.h"
#include "neigh_list_kokkos.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NPairCopyKokkos<DeviceType>::NPairCopyKokkos(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
------------------------------------------------------------------------- */

template<class DeviceType>
void NPairCopyKokkos<DeviceType>::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  if (list->kokkos) {
    if (!listcopy->kokkos)
      error->all(FLERR,"Cannot copy non-Kokkos neighbor list to Kokkos neighbor list");
    copy_to_kokkos(list);
  } else {
    if (!listcopy->kokkos)
      error->all(FLERR,"Missing Kokkos neighbor list for copy");
    copy_to_cpu(list);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NPairCopyKokkos<DeviceType>::copy_to_kokkos(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  list->inum = listcopy->inum;
  list->gnum = listcopy->gnum;
  list->ilist = listcopy->ilist;
  list->numneigh = listcopy->numneigh;
  list->ipage = listcopy->ipage;

  NeighListKokkos<DeviceType>* list_kk = (NeighListKokkos<DeviceType>*) list;
  NeighListKokkos<DeviceType>* listcopy_kk = (NeighListKokkos<DeviceType>*) list->listcopy;

  list_kk->d_ilist = listcopy_kk->d_ilist;
  list_kk->d_numneigh = listcopy_kk->d_numneigh;
  list_kk->d_neighbors = listcopy_kk->d_neighbors;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NPairCopyKokkos<DeviceType>::copy_to_cpu(NeighList *list)
{
  NeighList *listcopy = list->listcopy;
  NeighListKokkos<DeviceType>* listcopy_kk = (NeighListKokkos<DeviceType>*) listcopy;

  listcopy_kk->k_ilist.template sync<LMPHostType>();

  int inum = listcopy->inum;
  int gnum = listcopy->gnum;
  int inum_all = inum;
  if (list->ghost) inum_all += gnum;
  auto h_ilist = listcopy_kk->k_ilist.h_view;
  auto h_numneigh = Kokkos::create_mirror_view_and_copy(LMPHostType(),listcopy_kk->d_numneigh);
  auto h_neighbors = Kokkos::create_mirror_view_and_copy(LMPHostType(),listcopy_kk->d_neighbors);

  list->inum = inum;
  list->gnum = gnum;
  auto ilist = list->ilist;
  auto numneigh = list->numneigh;

  // Kokkos neighbor data is stored differently than regular CPU,
  //  must copy element by element

  int *neighptr;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  ipage->reset();

  for (int ii = 0; ii < inum_all; ii++) {
    neighptr = ipage->vget();

    const int i = h_ilist[ii];
    ilist[ii] = i;

    // loop over Kokkos neighbor list

    const int jnum = h_numneigh[i];
    numneigh[i] = jnum;

    for (int jj = 0; jj < jnum; jj++) {
      const int joriginal = h_neighbors(i,jj);
      neighptr[jj] = joriginal;
    }

    firstneigh[i] = neighptr;
    ipage->vgot(jnum);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}

namespace LAMMPS_NS {
template class NPairCopyKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class NPairCopyKokkos<LMPHostType>;
#endif
}
