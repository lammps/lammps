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
#include "kokkos.h"
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
  forward_comm_device = 1;

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
  int prev_auto_sync = lmp->kokkos->auto_sync;
  lmp->kokkos->auto_sync = 0;

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

  comm->forward_comm(this);

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
  atomKK->k_v.clear_sync_state();

  lmp->kokkos->auto_sync = prev_auto_sync;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeAveSphereAtomKokkos<DeviceType>::operator()(TagComputeAveSphereAtom, const int &ii) const
{
  KK_FLOAT massone_i,massone_j;

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    if (rmass.data()) massone_i = rmass[i];
    else massone_i = mass[type[i]];

    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    // i atom contribution

    int count = 1;
    KK_ACC_FLOAT totalmass = massone_i;
    KK_ACC_FLOAT p[3];
    p[0] = v(i,0)*massone_i;
    p[1] = v(i,1)*massone_i;
    p[2] = v(i,2)*massone_i;

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      if (rmass.data()) massone_j = rmass[j];
      else massone_j = mass[type[j]];

      const KK_FLOAT delx = x(j,0) - xtmp;
      const KK_FLOAT dely = x(j,1) - ytmp;
      const KK_FLOAT delz = x(j,2) - ztmp;
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        count++;
        totalmass += massone_j;
        p[0] += v(j,0)*massone_j;
        p[1] += v(j,1)*massone_j;
        p[2] += v(j,2)*massone_j;
      }
    }

    KK_FLOAT vcom[3];
    vcom[0] = p[0]/totalmass;
    vcom[1] = p[1]/totalmass;
    vcom[2] = p[2]/totalmass;

    // i atom contribution

    KK_FLOAT vnet[3];
    vnet[0] = v(i,0) - vcom[0];
    vnet[1] = v(i,1) - vcom[1];
    vnet[2] = v(i,2) - vcom[2];
    KK_ACC_FLOAT ke_sum = massone_i * (vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2]);

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      if (rmass.data()) massone_j = rmass[j];
      else massone_j = mass[type[j]];

      const KK_FLOAT delx = x(j,0) - xtmp;
      const KK_FLOAT dely = x(j,1) - ytmp;
      const KK_FLOAT delz = x(j,2) - ztmp;
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        vnet[0] = v(j,0) - vcom[0];
        vnet[1] = v(j,1) - vcom[1];
        vnet[2] = v(j,2) - vcom[2];
        ke_sum += massone_j * (vnet[0]*vnet[0] + vnet[1]*vnet[1] + vnet[2]*vnet[2]);
      }
    }
    KK_FLOAT density = mv2d*totalmass/volume;
    KK_FLOAT temp = mvv2e*ke_sum/(adof*count*boltz);
    d_result(i,0) = density;
    d_result(i,1) = temp;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int ComputeAveSphereAtomKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                         DAT::tdual_double_1d &k_buf,
                                                         int pbc_flag, int* pbc)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  d_buf = k_buf.view<DeviceType>();

  atomKK->sync(execution_space,V_MASK);
  v = atomKK->k_v.view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeAveSphereAtomPackForwardComm>(0,n),*this);
  copymode = 0;

  return n*3;
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeAveSphereAtomKokkos<DeviceType>::operator()(TagComputeAveSphereAtomPackForwardComm, const int &i) const {
  const int j = d_sendlist(i);

  d_buf[3*i] = static_cast<double>(v(j,0));
  d_buf[3*i+1] = static_cast<double>(v(j,1));
  d_buf[3*i+2] = static_cast<double>(v(j,2));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeAveSphereAtomKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_double_1d &buf)
{
  first = first_in;
  d_buf = buf.view<DeviceType>();

  atomKK->sync(execution_space,V_MASK);
  v = atomKK->k_v.view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeAveSphereAtomUnpackForwardComm>(0,n),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeAveSphereAtomKokkos<DeviceType>::operator()(TagComputeAveSphereAtomUnpackForwardComm, const int &i) const {
  v(i + first,0) = static_cast<KK_FLOAT>(d_buf[3*i]);
  v(i + first,1) = static_cast<KK_FLOAT>(d_buf[3*i+1]);
  v(i + first,2) = static_cast<KK_FLOAT>(d_buf[3*i+2]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int ComputeAveSphereAtomKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  atomKK->sync(Host,V_MASK);

  int m = ComputeAveSphereAtom::pack_forward_comm(n,list,buf,pbc_flag,pbc);

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeAveSphereAtomKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  atomKK->sync(Host,V_MASK);

  ComputeAveSphereAtom::unpack_forward_comm(n,first,buf);

  atomKK->modified(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeAveSphereAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeAveSphereAtomKokkos<LMPHostType>;
#endif
}
