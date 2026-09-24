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

#include "compute_hexorder_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeHexOrderAtomKokkos<DeviceType>::ComputeHexOrderAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeHexOrderAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeHexOrderAtomKokkos<DeviceType>::~ComputeHexOrderAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_qnarray,qnarray);
  qnarray = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeHexOrderAtomKokkos<DeviceType>::init()
{
  ComputeHexOrderAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeHexOrderAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_qnarray,qnarray);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_qnarray,qnarray,nmax,ncol,"hexorder/atom:qnarray");
    d_qnarray = k_qnarray.template view<DeviceType>();
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  auto *k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // one distsq/nearest scratch row per atom

  int maxneigh_kk = 0;
  auto h_numneigh = Kokkos::create_mirror_view(d_numneigh);
  Kokkos::deep_copy(h_numneigh,d_numneigh);
  auto h_ilist = Kokkos::create_mirror_view(d_ilist);
  Kokkos::deep_copy(h_ilist,d_ilist);
  for (int ii = 0; ii < inum; ii++)
    maxneigh_kk = MAX(maxneigh_kk,h_numneigh(h_ilist(ii)));
  maxneigh = maxneigh_kk;

  if ((int)d_distsq.extent(0) < inum || (int)d_distsq.extent(1) < maxneigh_kk) {
    d_distsq = typename AT::t_kkfloat_2d("hexorder/atom:distsq",inum,maxneigh_kk);
    d_nearest = typename AT::t_int_2d("hexorder/atom:nearest",inum,maxneigh_kk);
  }

  cutsq_kk = static_cast<KK_FLOAT>(cutsq);
  groupbit_kk = groupbit;

  atomKK->sync(execution_space,datamask_read);
  x = atomKK->k_x.template view<DeviceType>();
  mask = atomKK->k_mask.template view<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for("ComputeHexOrderAtom",
    Kokkos::RangePolicy<DeviceType, TagComputeHexOrderAtom>(0,inum),*this);

  copymode = 0;

  k_qnarray.template modify<DeviceType>();
  k_qnarray.sync_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeHexOrderAtomKokkos<DeviceType>::operator()(TagComputeHexOrderAtom, const int &ii) const
{
  const int i = d_ilist[ii];

  if (!(mask[i] & groupbit_kk)) {
    d_qnarray(i,0) = d_qnarray(i,1) = static_cast<KK_FLOAT>(0.0);
    return;
  }

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  int ncount = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < cutsq_kk) {
      d_distsq(ii,ncount) = rsq;
      d_nearest(ii,ncount++) = j;
    }
  }

  // fewer than nnn neighbors (or no neighbor at all with nnn = NULL):
  // order parameter = 0

  if ((ncount < nnn) || (ncount == 0)) {
    d_qnarray(i,0) = d_qnarray(i,1) = static_cast<KK_FLOAT>(0.0);
    return;
  }

  if (nnn > 0) {
    select2_kk(nnn,ncount,ii);
    ncount = nnn;
  }

  KK_FLOAT usum = static_cast<KK_FLOAT>(0.0);
  KK_FLOAT vsum = static_cast<KK_FLOAT>(0.0);

  for (int jj = 0; jj < ncount; jj++) {
    int j = d_nearest(ii,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);

    // z^ndegree for the unit complex number z = (delx,dely)/|(delx,dely)|,
    // which is what the CPU style's std::pow(complex,int) evaluates

    const KK_FLOAT ntheta = ndegree * Kokkos::atan2(dely,delx);
    usum += Kokkos::cos(ntheta);
    vsum += Kokkos::sin(ntheta);
  }

  // average over the neighbors actually used: ncount, which the branch above
  // has set to nnn when nnn > 0, but stays the in-cutoff count for nnn = NULL

  d_qnarray(i,0) = usum/ncount;
  d_qnarray(i,1) = vsum/ncount;
}

/* ----------------------------------------------------------------------
   device version of ComputeHexOrderAtom::select2(), operating on one row
   of the per-atom distsq and nearest scratch
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeHexOrderAtomKokkos<DeviceType>::select2_kk(int k, int n, int ii) const
{
  int i,ir,j,l,mid,ia,itmp;
  KK_FLOAT a,tmp;

  auto arr = Kokkos::subview(d_distsq, ii, Kokkos::ALL);
  auto iarr = Kokkos::subview(d_nearest, ii, Kokkos::ALL);

  l = 0;
  ir = n-1;
  while (true) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        tmp = arr[l]; arr[l] = arr[ir]; arr[ir] = tmp;
        itmp = iarr[l]; iarr[l] = iarr[ir]; iarr[ir] = itmp;
      }
      return;
    } else {
      mid = ((l+ir+2) >> 1) - 1;
      tmp = arr[mid]; arr[mid] = arr[l+1]; arr[l+1] = tmp;
      itmp = iarr[mid]; iarr[mid] = iarr[l+1]; iarr[l+1] = itmp;
      if (arr[l] > arr[ir]) {
        tmp = arr[l]; arr[l] = arr[ir]; arr[ir] = tmp;
        itmp = iarr[l]; iarr[l] = iarr[ir]; iarr[ir] = itmp;
      }
      if (arr[l+1] > arr[ir]) {
        tmp = arr[l+1]; arr[l+1] = arr[ir]; arr[ir] = tmp;
        itmp = iarr[l+1]; iarr[l+1] = iarr[ir]; iarr[ir] = itmp;
      }
      if (arr[l] > arr[l+1]) {
        tmp = arr[l]; arr[l] = arr[l+1]; arr[l+1] = tmp;
        itmp = iarr[l]; iarr[l] = iarr[l+1]; iarr[l+1] = itmp;
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      while (true) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
        itmp = iarr[i]; iarr[i] = iarr[j]; iarr[j] = itmp;
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j+1 >= k) ir = j-1;
      if (j+1 <= k) l = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeHexOrderAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeHexOrderAtomKokkos<LMPHostType>;
#endif
}
