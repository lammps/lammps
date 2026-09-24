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

#include "compute_centro_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCentroAtomKokkos<DeviceType>::ComputeCentroAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeCentroAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCentroAtomKokkos<DeviceType>::~ComputeCentroAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_centro,centro);
  centro = nullptr;
  if (axes_flag) {
    memoryKK->destroy_kokkos(k_array_atom,array_atom);
    array_atom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCentroAtomKokkos<DeviceType>::init()
{
  ComputeCentroAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCentroAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow centro array if necessary
  // grow array_atom array if axes_flag set

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_centro,centro);
    if (axes_flag) memoryKK->destroy_kokkos(k_array_atom,array_atom);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_centro,centro,nmax,"centro/atom:centro");
    d_centro = k_centro.template view<DeviceType>();
    if (!axes_flag) {
      vector_atom = centro;
    } else {
      memoryKK->create_kokkos(k_array_atom,array_atom,nmax,size_peratom_cols,
                              "centro/atom:array_atom");
      d_array_atom = k_array_atom.template view<DeviceType>();
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  auto *k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  nhalf = nnn / 2;
  npairs = nnn * (nnn - 1) / 2;

  // one distsq/nearest/pairs scratch row per atom

  int maxneigh_kk = 0;
  auto h_numneigh = Kokkos::create_mirror_view(d_numneigh);
  Kokkos::deep_copy(h_numneigh,d_numneigh);
  auto h_ilist = Kokkos::create_mirror_view(d_ilist);
  Kokkos::deep_copy(h_ilist,d_ilist);
  for (int ii = 0; ii < inum; ii++)
    maxneigh_kk = MAX(maxneigh_kk,h_numneigh(h_ilist(ii)));
  maxneigh = maxneigh_kk;

  if ((int)d_distsq.extent(0) < inum || (int)d_distsq.extent(1) < maxneigh_kk) {
    d_distsq = typename AT::t_kkfloat_2d("centro/atom:distsq",inum,maxneigh_kk);
    d_nearest = typename AT::t_int_2d("centro/atom:nearest",inum,maxneigh_kk);
  }
  if ((int)d_pairs.extent(0) < inum || (int)d_pairs.extent(1) < npairs)
    d_pairs = typename AT::t_kkfloat_2d("centro/atom:pairs",inum,npairs);

  cutsq_kk = static_cast<KK_FLOAT>(force->pair->cutforce * force->pair->cutforce);
  groupbit_kk = groupbit;

  atomKK->sync(execution_space,datamask_read);
  x = atomKK->k_x.template view<DeviceType>();
  mask = atomKK->k_mask.template view<DeviceType>();

  copymode = 1;

  if (axes_flag)
    Kokkos::parallel_for("ComputeCentroAtom",
      Kokkos::RangePolicy<DeviceType, TagComputeCentroAtom<1> >(0,inum),*this);
  else
    Kokkos::parallel_for("ComputeCentroAtom",
      Kokkos::RangePolicy<DeviceType, TagComputeCentroAtom<0> >(0,inum),*this);

  copymode = 0;

  k_centro.template modify<DeviceType>();
  k_centro.sync_host();
  if (axes_flag) {
    k_array_atom.template modify<DeviceType>();
    k_array_atom.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int AXES>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeCentroAtomKokkos<DeviceType>::operator()(TagComputeCentroAtom<AXES>, const int &ii) const
{
  const int i = d_ilist[ii];

  if (!(mask[i] & groupbit_kk)) {
    d_centro[i] = static_cast<KK_FLOAT>(0.0);
    if (AXES)
      for (int m = 1; m < 10; m++) d_array_atom(i,m) = static_cast<KK_FLOAT>(0.0);
    return;
  }

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  // distsq[] = squared distance to each neighbor inside the force cutoff
  // nearest[] = the corresponding atom indices

  int n = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < cutsq_kk) {
      d_distsq(ii,n) = rsq;
      d_nearest(ii,n++) = j;
    }
  }

  // fewer than nnn neighbors: centro = 0

  if (n < nnn) {
    d_centro[i] = static_cast<KK_FLOAT>(0.0);
    if (AXES) {
      for (int m = 1; m < 10; m++) d_array_atom(i,m) = static_cast<KK_FLOAT>(0.0);
      d_array_atom(i,0) = static_cast<KK_FLOAT>(0.0);
    }
    return;
  }

  // put the nnn nearest neighbors in the first nnn slots

  select2_kk(nnn,n,ii);

  // R = Ri + Rj for each of npairs i,j pairs among the nnn neighbors

  KK_FLOAT r1[3],r2[3];
  KK_FLOAT rsq1 = cutsq_kk;
  KK_FLOAT rsq2 = cutsq_kk;
  if (AXES) {
    r1[0] = r1[1] = r1[2] = static_cast<KK_FLOAT>(0.0);
    r2[0] = r2[1] = r2[2] = static_cast<KK_FLOAT>(0.0);
  }

  n = 0;
  for (int j = 0; j < nnn; j++) {
    const int jj = d_nearest(ii,j);
    for (int k = j+1; k < nnn; k++) {
      const int kk = d_nearest(ii,k);
      const KK_FLOAT delx = x(jj,0) + x(kk,0) - static_cast<KK_FLOAT>(2.0)*xtmp;
      const KK_FLOAT dely = x(jj,1) + x(kk,1) - static_cast<KK_FLOAT>(2.0)*ytmp;
      const KK_FLOAT delz = x(jj,2) + x(kk,2) - static_cast<KK_FLOAT>(2.0)*ztmp;
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      d_pairs(ii,n++) = rsq;

      if (AXES) {

        // rsq1, rsq2 are the two smallest R^2, r1 and r2 the matching Ri - Rj

        if (rsq < rsq2) {
          if (rsq < rsq1) {
            rsq2 = rsq1;
            r2[0] = r1[0]; r2[1] = r1[1]; r2[2] = r1[2];
            rsq1 = rsq;
            r1[0] = x(jj,0) - x(kk,0);
            r1[1] = x(jj,1) - x(kk,1);
            r1[2] = x(jj,2) - x(kk,2);
          } else {
            rsq2 = rsq;
            r2[0] = x(jj,0) - x(kk,0);
            r2[1] = x(jj,1) - x(kk,1);
            r2[2] = x(jj,2) - x(kk,2);
          }
        }
      }
    }
  }

  if (AXES) {
    KK_FLOAT r3[3];
    r3[0] = r1[1]*r2[2] - r1[2]*r2[1];
    r3[1] = r1[2]*r2[0] - r1[0]*r2[2];
    r3[2] = r1[0]*r2[1] - r1[1]*r2[0];

    KK_FLOAT s = Kokkos::sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
    s = static_cast<KK_FLOAT>(1.0)/s;
    r1[0] *= s; r1[1] *= s; r1[2] *= s;
    s = Kokkos::sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
    s = static_cast<KK_FLOAT>(1.0)/s;
    r2[0] *= s; r2[1] *= s; r2[2] *= s;
    s = Kokkos::sqrt(r3[0]*r3[0] + r3[1]*r3[1] + r3[2]*r3[2]);
    s = static_cast<KK_FLOAT>(1.0)/s;
    r3[0] *= s; r3[1] *= s; r3[2] *= s;

    for (int m = 0; m < 3; m++) {
      d_array_atom(i,1+m) = r1[m];
      d_array_atom(i,4+m) = r2[m];
      d_array_atom(i,7+m) = r3[m];
    }
  }

  // the nhalf smallest pair distances sum to the centrosymmetry parameter

  select_kk(nhalf,npairs,ii);

  KK_FLOAT value = static_cast<KK_FLOAT>(0.0);
  for (int j = 0; j < nhalf; j++) value += d_pairs(ii,j);
  d_centro[i] = value;

  if (AXES) d_array_atom(i,0) = value;
}

/* ----------------------------------------------------------------------
   device version of ComputeCentroAtom::select(), operating on one row of
   the per-atom pairs scratch
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeCentroAtomKokkos<DeviceType>::select_kk(int k, int n, int ii) const
{
  int i,ir,j,l,mid;
  KK_FLOAT a,tmp;

  auto arr = Kokkos::subview(d_pairs, ii, Kokkos::ALL);

  l = 0;
  ir = n-1;
  while (true) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) { tmp = arr[l]; arr[l] = arr[ir]; arr[ir] = tmp; }
      return;
    } else {
      mid = ((l+ir+2) >> 1) - 1;
      tmp = arr[mid]; arr[mid] = arr[l+1]; arr[l+1] = tmp;
      if (arr[l] > arr[ir]) { tmp = arr[l]; arr[l] = arr[ir]; arr[ir] = tmp; }
      if (arr[l+1] > arr[ir]) { tmp = arr[l+1]; arr[l+1] = arr[ir]; arr[ir] = tmp; }
      if (arr[l] > arr[l+1]) { tmp = arr[l]; arr[l] = arr[l+1]; arr[l+1] = tmp; }
      i = l+1;
      j = ir;
      a = arr[l+1];
      while (true) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      if (j+1 >= k) ir = j-1;
      if (j+1 <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   device version of ComputeCentroAtom::select2(), operating on one row of
   the per-atom distsq and nearest scratch
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeCentroAtomKokkos<DeviceType>::select2_kk(int k, int n, int ii) const
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
template class ComputeCentroAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeCentroAtomKokkos<LMPHostType>;
#endif
}
