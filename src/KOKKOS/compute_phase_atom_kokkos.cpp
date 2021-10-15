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

#include "compute_phase_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "pair.h"
#include "update.h"
#include "math_const.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputePhaseAtomKokkos<DeviceType>::ComputePhaseAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputePhaseAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputePhaseAtomKokkos<DeviceType>::~ComputePhaseAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_phase,phase);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputePhaseAtomKokkos<DeviceType>::init()
{
  ComputePhaseAtom::init();

  // need an occasional full neighbor list

  // irequest = neigh request made by parent class

  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputePhaseAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow phase array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_phase,phase);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_phase,phase,nmax,2,"phase/atom:phase");
    d_phase = k_phase.view<DeviceType>();
    array_atom = phase;
  }

  // need velocities of ghost atoms

  atomKK->sync(Host,V_MASK);
  comm->forward_comm_compute(this);
  atomKK->modified(Host,V_MASK);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);
  int inum = list->inum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // compute phase for each atom in group
  // use full neighbor list to count atoms less than cutoff

  atomKK->sync(execution_space,X_MASK|V_MASK|TYPE_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  Kokkos::deep_copy(d_phase,0.0);

  copymode = 1;
  typename Kokkos::RangePolicy<DeviceType, TagComputePhaseAtom> policy(0,inum);
  Kokkos::parallel_for("ComputePhaseAtom",policy,*this);
  copymode = 0;

  k_phase.modify<DeviceType>();
  k_phase.sync_host();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePhaseAtomKokkos<DeviceType>::operator()(TagComputePhaseAtom, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    // i atom contribution

    int count = 1;
    double vsum[3];
    vsum[0] = v(i,0);
    vsum[1] = v(i,1);
    vsum[2] = v(i,2);

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        count++;
        vsum[0] += v(j,0);
        vsum[1] += v(j,1);
        vsum[2] += v(j,2);
      }
    }

    double vavg[3];
    vavg[0] = vsum[0]/count;
    vavg[1] = vsum[1]/count;
    vavg[2] = vsum[2]/count;

    // i atom contribution

    count = 1;
    double vnet[3];
    vnet[0] = v(i,0) - vavg[0];
    vnet[1] = v(i,1) - vavg[1];
    vnet[2] = v(i,2) - vavg[2];
    double ke_sum = vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2];

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        count++;
        vnet[0] = v(j,0) - vavg[0];
        vnet[1] = v(j,1) - vavg[1];
        vnet[2] = v(j,2) - vavg[2];
        ke_sum += vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2];
      }
    }
    double density = count/sphere_vol;
    double temp = ke_sum/3.0/count;
    d_phase(i,0) = density;
    d_phase(i,1) = temp;
  }
}

namespace LAMMPS_NS {
template class ComputePhaseAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputePhaseAtomKokkos<LMPHostType>;
#endif
}
