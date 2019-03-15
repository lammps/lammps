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

#include "fix_neigh_history_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"
#include "comm.h"
#include "atom_vec_kokkos.h"

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
  npartner = NULL;
  partner = NULL;
  valuepartner = NULL;

  maxpartner = 8;
  grow_arrays(atom->nmax);

  d_resize = typename ArrayTypes<DeviceType>::t_int_scalar("FixNeighHistoryKokkos::resize");
#ifndef KOKKOS_USE_CUDA_UVM
  h_resize = Kokkos::create_mirror_view(d_resize);
#else
  h_resize = d_resize;
#endif
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
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::pre_exchange()
{
  copymode = 1;

  h_resize() = 1;
  while (h_resize() > 0) {
    FixNeighHistoryKokkosZeroPartnerCountFunctor<DeviceType> zero(this);
    Kokkos::parallel_for(nlocal_neigh,zero);

    h_resize() = 0;
    deep_copy(d_resize, h_resize);

    FixNeighHistoryKokkosPreExchangeFunctor<DeviceType> f(this);
    Kokkos::parallel_for(nlocal_neigh,f);

    deep_copy(h_resize, d_resize);
    if (h_resize() > 0) {
      maxpartner += 8;
      memoryKK->grow_kokkos(k_partner,partner,atom->nmax,maxpartner,"neighbor_history:partner");
      memoryKK->grow_kokkos(k_valuepartner,valuepartner,atom->nmax,dnum*maxpartner,"neighbor_history:valuepartner");
    }
  }

  copymode = 0;

  comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum+1)*maxpartner+1);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::zero_partner_count_item(const int &i) const
{
  d_npartner[i] = 0;
}

template <class DeviceType>
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

  if (maxatom < nlocal || k_list->maxneighs > d_firstflag.extent(1)) {
    maxatom = nall;
    d_firstflag = Kokkos::View<int**>("neighbor_history:firstflag",maxatom,k_list->maxneighs);
    d_firstvalue = Kokkos::View<LMP_FLOAT**>("neighbor_history:firstvalue",maxatom,k_list->maxneighs*dnum);
  }

  copymode = 1;

  FixNeighHistoryKokkosPostNeighborFunctor<DeviceType> f(this);
  Kokkos::parallel_for(inum,f);

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
  double bytes = d_firstflag.extent(0)*d_firstflag.extent(1)*sizeof(int);
  bytes += d_firstvalue.extent(0)*d_firstvalue.extent(1)*sizeof(double);
  bytes += 2*k_npartner.extent(0)*sizeof(int);
  bytes += 2*k_partner.extent(0)*k_partner.extent(1)*sizeof(int);
  bytes += 2*k_valuepartner.extent(0)*k_valuepartner.extent(1)*sizeof(double);
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
void FixNeighHistoryKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
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

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixNeighHistoryKokkos_ExchangeFirstPartnerFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _npartner;
  typename AT::t_xfloat_1d_um _firstpartner;
  typename AT::t_int_scalar _count;
  const int _nsend;
  const int _dnum;

  FixNeighHistoryKokkos_ExchangeFirstPartnerFunctor(
    const typename AT::tdual_int_1d &sendlist,
    const typename AT::tdual_int_1d &npartner,
    const typename AT::t_xfloat_1d_um &firstpartner,
    const typename AT::tdual_int_scalar &count,
    const int &nsend,
    const int &dnum):
    _sendlist(sendlist.template view<DeviceType>()),
    _npartner(npartner.template view<DeviceType>()),
    _firstpartner(firstpartner),
    _count(count.template view<DeviceType>()),
    _nsend(nsend),
    _dnum(dnum)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, int &update, const bool &final) const {
    const int n = 1+_npartner(_sendlist(i))*(_dnum+1);
    if (final) {
      _firstpartner(i) = d_ubuf(_nsend+update).d;
      if (i == _nsend - 1)
        _count() = _nsend+update+n;
    }
    update += n;
  }
};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixNeighHistoryKokkos_PackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  typename AT::t_int_1d _npartner;
  typename AT::t_tagint_2d _partner;
  typename AT::t_float_2d _valuepartner;
  typename AT::t_xfloat_1d_um _firstpartner;
  typename AT::t_xfloat_1d_um _buf;
  const int _dnum;

  FixNeighHistoryKokkos_PackExchangeFunctor(
    const typename AT::tdual_int_1d &sendlist,
    const typename AT::tdual_int_1d &copylist,
    const typename AT::tdual_int_1d &npartner,
    const typename AT::tdual_tagint_2d &partner,
    const typename AT::tdual_float_2d &valuepartner,
    const typename AT::t_xfloat_1d_um &firstpartner,
    const typename AT::t_xfloat_1d_um &buf,
    const int &dnum):
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _npartner(npartner.template view<DeviceType>()),
    _partner(partner.template view<DeviceType>()),
    _valuepartner(valuepartner.template view<DeviceType>()),
    _firstpartner(firstpartner),
    _buf(buf),
    _dnum(dnum)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &mysend) const {
    const int i = _sendlist(mysend);
    const int n = _npartner(i);
    int m = (int) d_ubuf(_firstpartner(mysend)).i;
    _buf(m++) = d_ubuf(n).d;
    for (int p = 0; p < n; p++) {
      _buf(m++) = d_ubuf(_partner(i,p)).d;
      for (int v = 0; v < _dnum; v++) {
        _buf(m++) = _valuepartner(i,_dnum*p+v);
      }
    }
    const int j = _copylist(mysend);
    if (j > -1) {
      const int nj = _npartner(j);
      _npartner(i) = nj;
      for (int p = 0; p < nj; p++) {
        _partner(i,p) = _partner(j,p);
        for (int v = 0; v < _dnum; v++) {
          _valuepartner(i,_dnum*p+v) = _valuepartner(j,_dnum*p+v);
        }
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend,DAT::tdual_xfloat_2d &buf,
   DAT::tdual_int_1d k_sendlist,
   DAT::tdual_int_1d k_copylist,
   ExecutionSpace space, int dim,
   X_FLOAT lo, X_FLOAT hi)
{
  k_npartner.template sync<DeviceType>();
  k_partner.template sync<DeviceType>();
  k_valuepartner.template sync<DeviceType>();

  typename ArrayTypes<DeviceType>::t_xfloat_1d_um d_firstpartner(
    buf.template view<DeviceType>().data(),
    buf.extent(0)*buf.extent(1));
  typename ArrayTypes<DeviceType>::tdual_int_scalar k_count("neighbor_history:k_count");

  k_count.h_view() = 0;
  if (space == Device) {
    k_count.template modify<LMPHostType>();
    k_count.template sync<LMPDeviceType>();
  }

  Kokkos::parallel_scan(
    nsend,
    FixNeighHistoryKokkos_ExchangeFirstPartnerFunctor<DeviceType>(
      k_sendlist,k_npartner,d_firstpartner,k_count,nsend,dnum));

  if (space == Device) {
    k_count.template modify<LMPDeviceType>();
    k_count.template sync<LMPHostType>();
  }

  Kokkos::parallel_for(
    nsend,
    FixNeighHistoryKokkos_PackExchangeFunctor<DeviceType>(
      k_sendlist,k_copylist,k_npartner,k_partner,k_valuepartner,
      d_firstpartner,d_firstpartner,dnum));

  return k_count.h_view();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixNeighHistoryKokkos_UnpackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_xfloat_1d_um _buf;
  typename AT::t_int_1d _npartner;
  typename AT::t_tagint_2d _partner;
  typename AT::t_float_2d _valuepartner;
  typename AT::t_int_1d _indices;
  const int _dnum;

  FixNeighHistoryKokkos_UnpackExchangeFunctor(
    const typename AT::tdual_xfloat_2d buf,
    const typename AT::tdual_int_1d &npartner,
    const typename AT::tdual_tagint_2d &partner,
    const typename AT::tdual_float_2d &valuepartner,
    const typename AT::tdual_int_1d &indices,
    const int &dnum):
    _npartner(npartner.template view<DeviceType>()),
    _partner(partner.template view<DeviceType>()),
    _valuepartner(valuepartner.template view<DeviceType>()),
    _indices(indices.template view<DeviceType>()),
    _dnum(dnum)
  {
    _buf = typename AT::t_xfloat_1d_um(buf.template view<DeviceType>().data(),buf.extent(0)*buf.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    int index = _indices(i);
    if (index > 0) {
      int m = (int) d_ubuf(_buf(i)).i;
      int n = (int) d_ubuf(_buf(m++)).i;
      _npartner(index) = n;
      for (int p = 0; p < n; p++) {
        _partner(index,p) = (tagint) d_ubuf(_buf(m++)).i;
        for (int v = 0; v < _dnum; v++) {
          _valuepartner(index,_dnum*p+v) = _buf(m++);
        }
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf,DAT::tdual_int_1d &indices,int nrecv,
  int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
  ExecutionSpace space)
{
  Kokkos::parallel_for(
    nrecv/16,
    FixNeighHistoryKokkos_UnpackExchangeFunctor<DeviceType>(
      k_buf,k_npartner,k_partner,k_valuepartner,indices,dnum));

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
#ifdef KOKKOS_HAVE_CUDA
template class FixNeighHistoryKokkos<LMPHostType>;
#endif
}
