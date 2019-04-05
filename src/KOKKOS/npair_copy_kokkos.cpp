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

#include "npair_copy_kokkos.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
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

namespace LAMMPS_NS {
template class NPairCopyKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class NPairCopyKokkos<LMPHostType>;
#endif
}
