/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU Gdirectoryeneral Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Megan McCarthy (SNL)
------------------------------------------------------------------------- */

#include "compute_composition_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCompositionAtomKokkos<DeviceType>::ComputeCompositionAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeCompositionAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCompositionAtomKokkos<DeviceType>::~ComputeCompositionAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_result,result);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCompositionAtomKokkos<DeviceType>::init()
{
  ComputeCompositionAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCompositionAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  int size_peratom_cols = 1 + atom->ntypes;
  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_result,result);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_result,result,nmax,size_peratom_cols,"composition/atom:result");
    d_result = k_result.view<DeviceType>();
    array_atom = result;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);
  int inum = list->inum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  atomKK->sync(execution_space,X_MASK|TYPE_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  ntypes = atom->ntypes;

  Kokkos::deep_copy(d_result,0.0);

  copymode = 1;
  typename Kokkos::RangePolicy<DeviceType, TagComputeCompositionAtom> policy(0,inum);
  Kokkos::parallel_for("ComputeLocalComp",policy,*this);
  copymode = 0;

  k_result.modify<DeviceType>();
  k_result.sync_host();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeCompositionAtomKokkos<DeviceType>::operator()(TagComputeCompositionAtom, const int &ii) const
{
  const int i = d_ilist[ii];

  if (mask[i] & groupbit) {

    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    // i atom contribution

    int count = 1.0;
    int itype = type[i];
    d_result(i,itype)++;

    for (int jj = 0; jj < jnum; jj++) {

      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      int jtype = type[j];

      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        count++;
        d_result(i,jtype) += 1.0;
      }
    }

    // total count of atoms found in sampled radius range

    d_result(i,0) = count;

    // local comp fractions per atom type

    double lfac = 1.0 / count;

    for (int n = 1; n < size_peratom_cols; n++) {
      d_result(i,n) *= lfac;
    }
  }
}

namespace LAMMPS_NS {
template class ComputeCompositionAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeCompositionAtomKokkos<LMPHostType>;
#endif
}
