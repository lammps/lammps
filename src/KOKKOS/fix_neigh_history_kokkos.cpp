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
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
FixNeighHistoryKokkos<Space>::FixNeighHistoryKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNeighHistory(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;

  memory->destroy(npartner);
  memory->sfree(partner);
  memory->sfree(valuepartner);
  npartner = NULL;
  partner = NULL;
  valuepartner = NULL;

  maxpartner = 8;
  grow_arrays(atom->nmax);

  d_resize = typename AT::t_int_scalar("FixNeighHistoryKokkos::resize");
#ifndef KOKKOS_USE_CUDA_UVM
  h_resize = Kokkos::create_mirror_view(d_resize);
#else
  h_resize = d_resize;
#endif
  h_resize() = 1;
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
FixNeighHistoryKokkos<Space>::~FixNeighHistoryKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_npartner, npartner);
  memoryKK->destroy_kokkos(k_partner, partner);
  memoryKK->destroy_kokkos(k_valuepartner, valuepartner);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::init()
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

template <ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::pre_exchange()
{
  copymode = 1;

  h_resize() = 1;
  while (h_resize() > 0) {
    FixNeighHistoryKokkosZeroPartnerCountFunctor<Space> zero(this);
    Kokkos::parallel_for(nlocal_neigh,zero);

    h_resize() = 0;
    deep_copy(d_resize, h_resize);

    FixNeighHistoryKokkosPreExchangeFunctor<Space> f(this);
    Kokkos::parallel_for(nlocal_neigh,f);

    deep_copy(h_resize, d_resize);
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

template <ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<Space>::zero_partner_count_item(const int &i) const
{
  d_npartner[i] = 0;
}

template <ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<Space>::pre_exchange_item(const int &ii) const
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

template <ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::setup_post_neighbor()
{
  post_neighbor();
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::post_neighbor()
{
  tag = DualViewHelper<Space>::view(atomKK->k_tag);

  int inum = pair->list->inum;
  NeighListKokkos<Space>* k_list = static_cast<NeighListKokkos<Space>*>(pair->list);
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
    d_firstvalue = Kokkos::View<typename AT::t_float_1d::data_type*>("neighbor_history:firstvalue",maxatom,k_list->maxneighs*dnum);
  }

  copymode = 1;

  FixNeighHistoryKokkosPostNeighborFunctor<Space> f(this);
  Kokkos::parallel_for(inum,f);

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNeighHistoryKokkos<Space>::post_neighbor_item(const int &ii) const
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

template<ExecutionSpace Space>
double FixNeighHistoryKokkos<Space>::memory_usage()
{
  KK_FLOAT bytes = d_firstflag.extent(0)*d_firstflag.extent(1)*sizeof(int);
  bytes += d_firstvalue.extent(0)*d_firstvalue.extent(1)*sizeof(KK_FLOAT);
  bytes += 2*k_npartner.extent(0)*sizeof(int);
  bytes += 2*k_partner.extent(0)*k_partner.extent(1)*sizeof(int);
  bytes += 2*k_valuepartner.extent(0)*k_valuepartner.extent(1)*sizeof(KK_FLOAT);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::grow_arrays(int nmax)
{
  k_npartner.sync_host(); // force reallocation on host
  k_partner.sync_host();
  k_valuepartner.sync_host();

  memoryKK->grow_kokkos(k_npartner,npartner,nmax,"neighbor_history:npartner");
  memoryKK->grow_kokkos(k_partner,partner,nmax,maxpartner,"neighbor_history:partner");
  memoryKK->grow_kokkos(k_valuepartner,valuepartner,nmax,dnum*maxpartner,"neighbor_history:valuepartner");

  d_npartner = DualViewHelper<Space>::view(k_npartner);
  d_partner = DualViewHelper<Space>::view(k_partner);
  d_valuepartner = DualViewHelper<Space>::view(k_valuepartner);

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNeighHistoryKokkos<Space>::copy_arrays(int i, int j, int delflag)
{
  k_npartner.sync_host();
  k_partner.sync_host();
  k_valuepartner.sync_host();

  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[i]; m++) {
    partner[j][m] = partner[i][m];
    valuepartner[j][m] = valuepartner[i][m];
  }

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
int FixNeighHistoryKokkos<Space>::pack_exchange(int i, double *buf)
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

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
int FixNeighHistoryKokkos<Space>::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  npartner[nlocal] = static_cast<int>(buf[n++]);
  for (int m = 0; m < npartner[nlocal]; m++) partner[nlocal][m] = static_cast<int>(buf[n++]);
  for (int m = 0; m < dnum*npartner[nlocal]; m++) valuepartner[nlocal][m] = buf[n++];

  k_npartner.modify_host();
  k_partner.modify_host();
  k_valuepartner.modify_host();

  return n;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNeighHistoryKokkos<Device>;
template class FixNeighHistoryKokkos<Host>;
}
