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

#include "fix_neigh_history_kokkos.h"

#include "atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixNeighHistoryKokkos<DeviceType>::FixNeighHistoryKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNeighHistory(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  memory->destroy(npartner);
  memory->sfree(partner);
  memory->sfree(valuepartner);
  npartner = nullptr;
  partner = nullptr;
  valuepartner = nullptr;

  maxpartner = 8;
  grow_arrays(atom->nmax);

  d_resize = typename ArrayTypes<DeviceType>::t_int_scalar("FixNeighHistoryKokkos::resize");
  h_resize = Kokkos::create_mirror_view(d_resize);
  h_resize() = 1;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixNeighHistoryKokkos<DeviceType>::~FixNeighHistoryKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_npartner, npartner);
  memoryKK->destroy_kokkos(k_partner, partner);
  memoryKK->destroy_kokkos(k_valuepartner, valuepartner);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::init()
{
  if (atomKK->tag_enable == 0)
    error->all(FLERR,"Neighbor history requires atoms have IDs");

  // this fix must come before any fix which migrates atoms in its pre_exchange()
  // b/c this fix's pre_exchange() creates per-atom data structure
  // that data must be current for atom migration to carry it along

  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i] == this) break;
    if (modify->fix[i]->pre_exchange_migrate)
      error->all(FLERR,"Fix neigh_history comes after a fix which "
                 "migrates atoms in pre_exchange");
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::pre_exchange()
{
  copymode = 1;

  k_firstflag.sync<DeviceType>();
  k_firstvalue.sync<DeviceType>();

  h_resize() = 1;
  while (h_resize() > 0) {
    FixNeighHistoryKokkosZeroPartnerCountFunctor<DeviceType> zero(this);
    Kokkos::parallel_for(nlocal_neigh,zero);

    h_resize() = 0;
    Kokkos::deep_copy(d_resize, h_resize);

    FixNeighHistoryKokkosPreExchangeFunctor<DeviceType> f(this);
    Kokkos::parallel_for(nlocal_neigh,f);

    Kokkos::deep_copy(h_resize, d_resize);
    if (h_resize() > 0) {
      maxpartner += 8;
      memoryKK->grow_kokkos(k_partner,partner,atom->nmax,maxpartner,"neighbor_history:partner");
      memoryKK->grow_kokkos(k_valuepartner,valuepartner,atom->nmax,dnum*maxpartner,"neighbor_history:valuepartner");
    }
  }

  copymode = 0;

  maxexchange = (dnum+1)*maxpartner+1;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::zero_partner_count_item(const int &i) const
{
  d_npartner[i] = 0;
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::pre_exchange_item(const int &ii) const
{
  const int i = d_ilist[ii];
  const int jnum = d_numneigh[i];

  for (int jj = 0; jj < jnum; jj++) {
    if (d_firstflag(i,jj)) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      int m = Kokkos::atomic_fetch_add(&d_npartner[i],1);
      if (m < maxpartner) {
        d_partner(i,m) = tag[j];
        for (int k = 0; k < dnum; k++)
          d_valuepartner(i,dnum*m+k) = d_firstvalue(i,dnum*jj+k);
      } else {
        d_resize() = 1;
      }
      if (j < nlocal_neigh) {
        m = Kokkos::atomic_fetch_add(&d_npartner[j],1);
        if (m < maxpartner) {
          d_partner(j,m) = tag[i];
          for (int k = 0; k < dnum; k++)
            d_valuepartner(j,dnum*m+k) = d_firstvalue(i,dnum*jj+k);
        } else {
          d_resize() = 1;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::setup_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::post_neighbor()
{
  tag = atomKK->k_tag.view<DeviceType>();
  atomKK->sync(execution_space,TAG_MASK);

  k_firstflag.sync<DeviceType>();
  k_firstvalue.sync<DeviceType>();

  int inum = pair->list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(pair->list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // store atom counts used for new neighbor list which was just built

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  nlocal_neigh = nlocal;
  nall_neigh = nall;

  // realloc firstflag and firstvalue if needed

  if (maxatom < nlocal || k_list->maxneighs > (int)d_firstflag.extent(1)) {
    maxatom = nall;
    k_firstflag = DAT::tdual_int_2d("neighbor_history:firstflag",maxatom,k_list->maxneighs);
    k_firstvalue = DAT::tdual_float_2d("neighbor_history:firstvalue",maxatom,k_list->maxneighs*dnum);
    d_firstflag = k_firstflag.view<DeviceType>();
    d_firstvalue = k_firstvalue.view<DeviceType>();
  }

  copymode = 1;

  FixNeighHistoryKokkosPostNeighborFunctor<DeviceType> f(this);
  Kokkos::parallel_for(inum,f);

  k_firstflag.modify<DeviceType>();
  k_firstvalue.modify<DeviceType>();

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::post_neighbor_item(const int &ii) const
{
  const int i = d_ilist[ii];
  const int jnum = d_numneigh[i];
  const int np = d_npartner[i];

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const int rflag = j >> SBBITS & 3;
    j &= NEIGHMASK;

    int m;
    if (rflag) {
      int jtag = tag(j);
      for (m = 0; m < np; m++)
        if (d_partner(i, m) == jtag) break;
      if (m < np) {
        d_firstflag(i,jj) = 1;
        for (int k = 0; k < dnum; k++) {
          d_firstvalue(i, dnum*jj+k) = d_valuepartner(i, dnum*m+k);
        }
      } else {
        d_firstflag(i,jj) = 0;
        for (int k = 0; k < dnum; k++) {
          d_firstvalue(i, dnum*jj+k) = 0;
        }
      }
    } else {
      d_firstflag(i,jj) = 0;
      for (int k = 0; k < dnum; k++) {
        d_firstvalue(i, dnum*jj+k) = 0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixNeighHistoryKokkos<DeviceType>::memory_usage()
{
  double bytes = (double)d_firstflag.extent(0)*d_firstflag.extent(1)*sizeof(int);
  bytes += (double)d_firstvalue.extent(0)*d_firstvalue.extent(1)*sizeof(double);
  bytes += (double)2*k_npartner.extent(0)*sizeof(int);
  bytes += (double)2*k_partner.extent(0)*k_partner.extent(1)*sizeof(int);
  bytes += (double)2*k_valuepartner.extent(0)*k_valuepartner.extent(1)*sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_npartner.template sync<LMPHostType>(); // force reallocation on host
  k_partner.template sync<LMPHostType>();
  k_valuepartner.template sync<LMPHostType>();

  memoryKK->grow_kokkos(k_npartner,npartner,nmax,"neighbor_history:npartner");
  memoryKK->grow_kokkos(k_partner,partner,nmax,maxpartner,"neighbor_history:partner");
  memoryKK->grow_kokkos(k_valuepartner,valuepartner,nmax,dnum*maxpartner,"neighbor_history:valuepartner");

  d_npartner = k_npartner.template view<DeviceType>();
  d_partner = k_partner.template view<DeviceType>();
  d_valuepartner = k_valuepartner.template view<DeviceType>();

  k_npartner.template modify<LMPHostType>();
  k_partner.template modify<LMPHostType>();
  k_valuepartner.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_npartner.template sync<LMPHostType>();
  k_partner.template sync<LMPHostType>();
  k_valuepartner.template sync<LMPHostType>();

  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[i]; m++) {
    partner[j][m] = partner[i][m];
    valuepartner[j][m] = valuepartner[i][m];
  }

  k_npartner.template modify<LMPHostType>();
  k_partner.template modify<LMPHostType>();
  k_valuepartner.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_npartner.template sync<LMPHostType>();
  k_partner.template sync<LMPHostType>();
  k_valuepartner.template sync<LMPHostType>();

  int n = 0;
  buf[n++] = npartner[i];
  for (int m = 0; m < npartner[i]; m++) buf[n++] = partner[i][m];
  for (int m = 0; m < dnum*npartner[i]; m++) buf[n++] = valuepartner[i][m];

  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  npartner[nlocal] = static_cast<int>(buf[n++]);
  for (int m = 0; m < npartner[nlocal]; m++) partner[nlocal][m] = static_cast<int>(buf[n++]);
  for (int m = 0; m < dnum*npartner[nlocal]; m++) valuepartner[nlocal][m] = buf[n++];

  k_npartner.template modify<LMPHostType>();
  k_partner.template modify<LMPHostType>();
  k_valuepartner.template modify<LMPHostType>();

  return n;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNeighHistoryKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNeighHistoryKokkos<LMPHostType>;
#endif
}
