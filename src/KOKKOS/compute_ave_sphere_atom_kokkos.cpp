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

#include "compute_ave_sphere_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeAveSphereAtomKokkos<DeviceType>::ComputeAveSphereAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeAveSphereAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeAveSphereAtomKokkos<DeviceType>::~ComputeAveSphereAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_result,result);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeAveSphereAtomKokkos<DeviceType>::init()
{
  ComputeAveSphereAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeAveSphereAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow result array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_result,result);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_result,result,nmax,2,"ave/sphere/atom:result");
    d_result = k_result.view<DeviceType>();
    array_atom = result;
  }

  // need velocities of ghost atoms

  atomKK->sync(Host,V_MASK);
  comm->forward_comm(this);
  atomKK->modified(Host,V_MASK);

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);
  int inum = list->inum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // compute properties for each atom in group
  // use full neighbor list to count atoms less than cutoff

  atomKK->sync(execution_space,X_MASK|V_MASK|RMASS_MASK|TYPE_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  adof = domain->dimension;
  mvv2e = force->mvv2e;
  mv2d = force->mv2d;
  boltz = force->boltz;

  Kokkos::deep_copy(d_result,0.0);

  copymode = 1;
  typename Kokkos::RangePolicy<DeviceType, TagComputeAveSphereAtom> policy(0,inum);
  Kokkos::parallel_for("ComputeAveSphereAtom",policy,*this);
  copymode = 0;

  k_result.modify<DeviceType>();
  k_result.sync_host();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeAveSphereAtomKokkos<DeviceType>::operator()(TagComputeAveSphereAtom, const int &ii) const
{
  double massone_i,massone_j;

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    if (rmass.data()) massone_i = rmass[i];
    else massone_i = mass[type[i]];

    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    // i atom contribution

    int count = 1;
    double totalmass = massone_i;
    double p[3];
    p[0] = v(i,0)*massone_i;
    p[1] = v(i,1)*massone_i;
    p[2] = v(i,2)*massone_i;

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      if (rmass.data()) massone_j = rmass[j];
      else massone_j = mass[type[j]];

      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        count++;
        totalmass += massone_j;
        p[0] += v(j,0)*massone_j;
        p[1] += v(j,1)*massone_j;
        p[2] += v(j,2)*massone_j;
      }
    }

    double vcom[3];
    vcom[0] = p[0]/totalmass;
    vcom[1] = p[1]/totalmass;
    vcom[2] = p[2]/totalmass;

    // i atom contribution

    double vnet[3];
    vnet[0] = v(i,0) - vcom[0];
    vnet[1] = v(i,1) - vcom[1];
    vnet[2] = v(i,2) - vcom[2];
    double ke_sum = massone_i * (vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2]);

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      if (rmass.data()) massone_j = rmass[j];
      else massone_j = mass[type[j]];

      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        vnet[0] = v(j,0) - vcom[0];
        vnet[1] = v(j,1) - vcom[1];
        vnet[2] = v(j,2) - vcom[2];
        ke_sum += massone_j * (vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2]);
      }
    }
    double density = mv2d*totalmass/volume;
    double temp = mvv2e*ke_sum/(adof*count*boltz);
    d_result(i,0) = density;
    d_result(i,1) = temp;
  }
}

namespace LAMMPS_NS {
template class ComputeAveSphereAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeAveSphereAtomKokkos<LMPHostType>;
#endif
}
