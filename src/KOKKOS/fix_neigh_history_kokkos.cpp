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
#include "atom_vec_kokkos.h"
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

  int inum = pair->list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(pair->list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  nlocal = atom->nlocal;

  h_resize() = 1;

  while (h_resize() > 0) {

    Kokkos::deep_copy(d_npartner,0);
    Kokkos::deep_copy(d_resize, 0);

    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryPreExchange>(0,inum),*this);

    Kokkos::deep_copy(h_resize, d_resize);

    if (h_resize()) {
      maxpartner += 8;
      memoryKK->grow_kokkos(k_partner,partner,atom->nmax,maxpartner,"neighbor_history:partner");
      memoryKK->grow_kokkos(k_valuepartner,valuepartner,atom->nmax,dnum*maxpartner,"neighbor_history:valuepartner");
    }
  }

  copymode = 0;

  maxexchange = (dnum+1)*maxpartner+1;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryPreExchange, const int &ii) const
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
      if (j < nlocal) {
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

  nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // realloc firstflag and firstvalue if needed

  if (maxatom < nlocal || k_list->maxneighs > (int)d_firstflag.extent(1)) {
    maxatom = atom->nmax;
    k_firstflag = DAT::tdual_int_2d("neighbor_history:firstflag",maxatom,k_list->maxneighs);
    k_firstvalue = DAT::tdual_float_2d("neighbor_history:firstvalue",maxatom,k_list->maxneighs*dnum);
    d_firstflag = k_firstflag.view<DeviceType>();
    d_firstvalue = k_firstvalue.view<DeviceType>();
  }

  copymode = 1;

  Kokkos::deep_copy(d_firstflag,0);
  Kokkos::deep_copy(d_firstvalue,0);

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryPostNeighbor>(0,inum),*this);

  k_firstflag.modify<DeviceType>();
  k_firstvalue.modify<DeviceType>();

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryPostNeighbor, const int &ii) const
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
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_npartner.sync_host(); // force reallocation on host
  k_partner.sync_host();
  k_valuepartner.sync_host();

  memoryKK->grow_kokkos(k_npartner,npartner,nmax,"neighbor_history:npartner");
  memoryKK->grow_kokkos(k_partner,partner,nmax,maxpartner,"neighbor_history:partner");
  memoryKK->grow_kokkos(k_valuepartner,valuepartner,nmax,dnum*maxpartner,"neighbor_history:valuepartner");

  d_npartner = k_npartner.template view<DeviceType>();
  d_partner = k_partner.template view<DeviceType>();
  d_valuepartner = k_valuepartner.template view<DeviceType>();

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_npartner.sync_host();
  k_partner.sync_host();
  k_valuepartner.sync_host();

  FixNeighHistory::copy_arrays(i,j,delflag);

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_npartner.sync_host();
  k_partner.sync_host();
  k_valuepartner.sync_host();

  return FixNeighHistory::pack_exchange(i,buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryFirstNeigh, const int &i, int &update, const bool &final) const {
  const int n = 1+d_npartner(d_sendlist(i))*(dnum+1);
  if (final) {
    d_firstpartner(i) = d_ubuf(nsend+update).d;
    if (i == nsend - 1)
      d_count() = nsend+update+n;
  }
  update += n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryPackExchange, const int &mysend) const {
  const int i = d_sendlist(mysend);
  const int n = d_npartner(i);
  int m = (int) d_ubuf(d_firstpartner(mysend)).i;
  d_firstpartner(m++) = d_ubuf(n).d;
  for (int p = 0; p < n; p++) {
    d_firstpartner(m++) = d_ubuf(d_partner(i,p)).d;
    for (int v = 0; v < dnum; v++) {
      d_firstpartner(m++) = d_valuepartner(i,dnum*p+v);
    }
  }
  const int j = d_copylist(mysend);
  if (j > -1) {
    const int nj = d_npartner(j);
    d_npartner(i) = nj;
    for (int p = 0; p < nj; p++) {
      d_partner(i,p) = d_partner(j,p);
      for (int v = 0; v < dnum; v++) {
	d_valuepartner(i,dnum*p+v) = d_valuepartner(j,dnum*p+v);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space, int dim, X_FLOAT lo, X_FLOAT hi)
{
  k_npartner.template sync<DeviceType>();
  k_partner.template sync<DeviceType>();
  k_valuepartner.template sync<DeviceType>();

  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();

  d_copylist = k_copylist.view<DeviceType>();
  this->nsend = nsend; 

  typename ArrayTypes<DeviceType>::t_xfloat_1d_um d_firstpartner(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));

  typename ArrayTypes<DeviceType>::tdual_int_scalar k_count("neighbor_history:k_count");

  k_count.h_view() = 0;
  k_count.modify_host();
  k_count.template sync<DeviceType>();

  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryFirstNeigh>(0,nsend),*this); 

  k_count.template modify<DeviceType>();
  k_count.sync_host();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryPackExchange>(0,nsend),*this);   

  return k_count.h_view();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryUnpackExchange, const int &i) const 
{
  int index = d_indices(i);
  if (index > 0) {
    int m = (int) d_ubuf(d_firstpartner(i)).i;
    int n = (int) d_ubuf(d_firstpartner(m++)).i;
    d_npartner(index) = n;
    for (int p = 0; p < n; p++) {
      d_partner(index,p) = (tagint) d_ubuf(d_firstpartner(m++)).i;
      for (int v = 0; v < dnum; v++) {
	d_valuepartner(index,dnum*p+v) = d_firstpartner(m++);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
  ExecutionSpace space)
{
  d_firstpartner = typename AT::t_xfloat_1d_um(k_buf.template view<DeviceType>().data(),k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  d_npartner = k_npartner.template view<DeviceType>();
  d_partner = k_partner.template view<DeviceType>();
  d_valuepartner = k_valuepartner.template view<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryUnpackExchange>(0,
    nrecv/(atom->avec->size_border + atom->avec->size_velocity + 2)),*this);

  k_npartner.template modify<DeviceType>();
  k_partner.template modify<DeviceType>();
  k_valuepartner.template modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int n = FixNeighHistory::unpack_exchange(nlocal,buf);

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();

  return n;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNeighHistoryKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNeighHistoryKokkos<LMPHostType>;
#endif
}
