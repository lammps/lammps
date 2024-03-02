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

template<class DeviceType>
FixNeighHistoryKokkos<DeviceType>::FixNeighHistoryKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNeighHistory(lmp, narg, arg)
{
  kokkosable = 1;
  exchange_comm_device = sort_device = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(npartner);
  memory->sfree(partner);
  memory->sfree(valuepartner);
  npartner = nullptr;
  partner = nullptr;
  valuepartner = nullptr;

  maxpartner = 8;
  grow_arrays(atom->nmax);

  d_resize = typename AT::t_int_scalar("fix_neigh_history::resize");
  h_resize = Kokkos::create_mirror_view(d_resize);

  d_count = typename AT::t_int_scalar("fix_neigh_history:count");
  h_count = Kokkos::create_mirror_view(d_count);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNeighHistoryKokkos<DeviceType>::~FixNeighHistoryKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_npartner, npartner);
  memoryKK->destroy_kokkos(k_partner, partner);
  memoryKK->destroy_kokkos(k_valuepartner, valuepartner);
}

/* ----------------------------------------------------------------------
   copy partner info from neighbor data structs (NDS) to atom arrays
   should be called whenever NDS store current history info
     and need to transfer the info to owned atoms
   e.g. when atoms migrate to new procs, new neigh list built, or between runs
     when atoms may be added or deleted (NDS becomes out-of-date)
   the next post_neighbor() will put this info back into new NDS
   called during run before atom exchanges, including for restart files
   called at end of run via post_run()
   do not call during setup of run (setup_pre_exchange)
     because there is no guarantee of a current NDS (even on continued run)
   if run command does a 2nd run with pre = no, then no neigh list
     will be built, but old neigh list will still have the info
   onesided and newton on and newton off versions
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::pre_exchange()
{
  if (onesided)
    error->all(FLERR,"Fix neigh/history/kk does not (yet) support onesided exchange communication");

  if (newton_pair)
    error->all(FLERR,"Fix neigh/history/kk requires newton 'off' for exchange communication");

  pre_exchange_no_newton();
}

/* ----------------------------------------------------------------------
   newton OFF version
   do not need partner values from ghost atoms
   assume J values are negative of I values
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::pre_exchange_no_newton()
{
  copymode = 1;

  k_firstflag.sync<DeviceType>();
  k_firstvalue.sync<DeviceType>();

  k_npartner.sync<DeviceType>();
  k_partner.sync<DeviceType>();
  k_valuepartner.sync<DeviceType>();

  // NOTE: all operations until very end are with nlocal_neigh <= current nlocal
  // because previous neigh list was built with nlocal_neigh
  // nlocal can be larger if other fixes added atoms at this pre_exchange()

  int inum = pair->list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(pair->list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

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

  maxexchange = (dnum+1)*maxpartner + 2;

  k_npartner.modify<DeviceType>();
  k_partner.modify<DeviceType>();
  k_valuepartner.modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

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
      if (j < nlocal_neigh) {
        m = Kokkos::atomic_fetch_add(&d_npartner[j],1);
        if (m < maxpartner) {
          d_partner(j,m) = tag[i];
          for (int k = 0; k < dnum; k++)
            d_valuepartner(j,dnum*m+k) = -d_firstvalue(i,dnum*jj+k);
        } else {
          d_resize() = 1;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::post_neighbor()
{
  tag = atomKK->k_tag.view<DeviceType>();
  atomKK->sync(execution_space,TAG_MASK);

  k_firstflag.sync<DeviceType>();
  k_firstvalue.sync<DeviceType>();

  k_npartner.sync<DeviceType>();
  k_partner.sync<DeviceType>();
  k_valuepartner.sync<DeviceType>();

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

  beyond_contact = pair->beyond_contact;

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

    int rflag;
    if (use_bit_flag) {
      rflag = histmask(j) | beyond_contact;
      j &= HISTMASK;
      d_firstflag(i,jj) = j;
    } else {
      rflag = 1;
    }

    // Remove special bond bits
    j &= NEIGHMASK;

    // rflag = 1 if r < radsum in npair_size() method or if pair interactions extend further
    // preserve neigh history info if tag[j] is in old-neigh partner list
    // this test could be more geometrically precise for two sphere/line/tri
    // if use_bit_flag is turned off, always record data since not all npair classes
    // apply a mask for history (and they could use the bits for special bonds)

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
  k_npartner.sync<DeviceType>(); // force reallocation on device
  k_partner.sync<DeviceType>();
  k_valuepartner.sync<DeviceType>();

  memoryKK->grow_kokkos(k_npartner,npartner,nmax,"neighbor_history:npartner");
  memoryKK->grow_kokkos(k_partner,partner,nmax,maxpartner,"neighbor_history:partner");
  memoryKK->grow_kokkos(k_valuepartner,valuepartner,nmax,dnum*maxpartner,"neighbor_history:valuepartner");

  d_npartner = k_npartner.template view<DeviceType>();
  d_partner = k_partner.template view<DeviceType>();
  d_valuepartner = k_valuepartner.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_npartner.sync_host();
  k_partner.sync_host();
  k_valuepartner.sync_host();

  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[i]; m++) partner[j][m] = partner[i][m];
  for (int m = 0; m < dnum*npartner[i]; m++) valuepartner[j][m] = valuepartner[i][m];

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_npartner.sync_device();
  k_partner.sync_device();
  k_valuepartner.sync_device();

  Sorter.sort(LMPDeviceType(), k_npartner.d_view);
  Sorter.sort(LMPDeviceType(), k_partner.d_view);
  Sorter.sort(LMPDeviceType(), k_valuepartner.d_view);

  k_npartner.modify_device();
  k_partner.modify_device();
  k_valuepartner.modify_device();
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

  int n = 0;
  buf[n++] = npartner[i];
  for (int m = 0; m < npartner[i]; m++) buf[n++] = partner[i][m];
  for (int m = 0; m < dnum*npartner[i]; m++) buf[n++] = valuepartner[i][m];

  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryPackExchange, const int &mysend, int &offset, const bool &final) const {

  const int i = d_sendlist(mysend);

  if (!final)
    offset += 1+d_npartner(i)*(dnum+1);
  else {
    int m = nsend + offset;

    d_buf(mysend) = d_ubuf(m).d;
    const int n = d_npartner(i);
    d_buf(m++) = d_ubuf(n).d;
    for (int p = 0; p < n; p++) {
      d_buf(m++) = d_ubuf(d_partner(i,p)).d;
      for (int v = 0; v < dnum; v++) {
        d_buf(m++) = d_valuepartner(i,dnum*p+v);
      }
    }
    if (mysend == nsend-1) d_count() = m;
    offset = m - nsend;

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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixNeighHistoryKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace /*space*/)
{
  k_npartner.template sync<DeviceType>();
  k_partner.template sync<DeviceType>();
  k_valuepartner.template sync<DeviceType>();

  k_buf.sync<DeviceType>();
  k_sendlist.sync<DeviceType>();
  k_copylist.sync<DeviceType>();

  d_sendlist = k_sendlist.view<DeviceType>();
  d_copylist = k_copylist.view<DeviceType>();
  this->nsend = nsend;

  d_buf = typename AT::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));

  Kokkos::deep_copy(d_count,0);

  copymode = 1;

  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryPackExchange>(0,nsend),*this);

  copymode = 0;

  k_npartner.modify<DeviceType>();
  k_partner.modify<DeviceType>();
  k_valuepartner.modify<DeviceType>();

  Kokkos::deep_copy(h_count,d_count);

  return h_count();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<DeviceType>::operator()(TagFixNeighHistoryUnpackExchange, const int &i) const
{
  int index = d_indices(i);

  if (index > -1) {
    int m = (int) d_ubuf(d_buf(i)).i;
    if (i >= nrecv1)
      m = nextrarecv1 + (int) d_ubuf(d_buf(nextrarecv1 + i - nrecv1)).i;

    int n = (int) d_ubuf(d_buf(m++)).i;
    d_npartner(index) = n;
    for (int p = 0; p < n; p++) {
      d_partner(index,p) = (tagint) d_ubuf(d_buf(m++)).i;
      for (int v = 0; v < dnum; v++) {
        d_valuepartner(index,dnum*p+v) = d_buf(m++);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNeighHistoryKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int nrecv1, int nextrarecv1,
  ExecutionSpace /*space*/)
{
  d_buf = typename AT::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  this->nrecv1 = nrecv1;
  this->nextrarecv1 = nextrarecv1;

  d_npartner = k_npartner.template view<DeviceType>();
  d_partner = k_partner.template view<DeviceType>();
  d_valuepartner = k_valuepartner.template view<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNeighHistoryUnpackExchange>(0,nrecv),*this);

  copymode = 0;

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
  k_npartner.sync_host();
  k_partner.sync_host();
  k_valuepartner.sync_host();

  int n = 0;
  npartner[nlocal] = static_cast<int>(buf[n++]);
  for (int m = 0; m < npartner[nlocal]; m++) partner[nlocal][m] = static_cast<tagint>(buf[n++]);
  for (int m = 0; m < dnum*npartner[nlocal]; m++) valuepartner[nlocal][m] = buf[n++];

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();

  return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixNeighHistoryKokkos<DeviceType>::memory_usage()
{
  double bytes = MemKK::memory_usage(d_partner);
  bytes += MemKK::memory_usage(d_valuepartner);
  bytes += MemKK::memory_usage(d_firstflag);
  bytes += MemKK::memory_usage(d_firstvalue);

  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNeighHistoryKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNeighHistoryKokkos<LMPHostType>;
#endif
}
