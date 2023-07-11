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

/* ----------------------------------------------------------------------
   Contributing authors: Naga Vydyanathan (NVIDIA), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_meam_kokkos.h"
#include "meam_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMEAMKokkos<DeviceType>::PairMEAMKokkos(LAMMPS *lmp) : PairMEAM(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  reverse_comm_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  delete meam_inst;
  meam_inst_kk = new MEAMKokkos<DeviceType>(memory);
  meam_inst = meam_inst_kk;
  myname = "meam/kk";
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMEAMKokkos<DeviceType>::~PairMEAMKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
  delete meam_inst_kk;
  meam_inst = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  // neighbor list info

  int inum_half = listhalf->inum;
  NeighListKokkos<DeviceType>* k_halflist = static_cast<NeighListKokkos<DeviceType>*>(listhalf);
  d_ilist_half = k_halflist->d_ilist;
  d_numneigh_half = k_halflist->d_numneigh;
  d_neighbors_half = k_halflist->d_neighbors;

  NeighListKokkos<DeviceType>* k_fulllist = static_cast<NeighListKokkos<DeviceType>*>(listfull);
  d_numneigh_full = k_fulllist->d_numneigh;
  d_neighbors_full = k_fulllist->d_neighbors;

  EV_FLOAT ev;

  copymode = 1;
  meam_inst_kk->copymode = 1;

  // strip neighbor lists of any special bond flags before using with MEAM
  // necessary before doing neigh_f2c and neigh_c2f conversions each step

  if (neighbor->ago == 0)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMNeighStrip >(0,inum_half),*this);

  // check size of scrfcn based on half neighbor list

  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;

  int n = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMEAMOffsets>(0,inum_half),*this,n);

  meam_inst_kk->meam_dens_setup(atom->nmax, nall, n);
  update_meam_views();

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  atomKK->sync(execution_space,datamask_read);

  int ntype = atom->ntypes;

  // 3 stages of MEAM calculation
  // loop over my atoms followed by communication

  int errorflag = 0;

  d_offset = typename AT::t_int_1d("pair:offset",inum_half+1);
  {
    // local variables for lambda capture

    auto l_ilist_half = d_ilist_half;
    auto l_numneigh_half = d_numneigh_half;
    auto l_offset = d_offset;

    Kokkos::parallel_scan(inum_half, LAMMPS_LAMBDA(int ii, int &m_fill, bool final) {
      int i = l_ilist_half[ii];
      m_fill += l_numneigh_half[i];
      if (final)
        l_offset[ii+1] = m_fill;
    });
  }

  int need_dup = lmp->kokkos->need_dup<DeviceType>();

  meam_inst_kk->meam_dens_init(inum_half,ntype,type,d_map,x,d_numneigh_half,d_numneigh_full,
    d_ilist_half,d_neighbors_half, d_neighbors_full, d_offset, neighflag, need_dup);

  meam_inst_kk->k_rho0.template modify<DeviceType>();
  meam_inst_kk->k_arho2b.template modify<DeviceType>();
  meam_inst_kk->k_arho1.template modify<DeviceType>();
  meam_inst_kk->k_arho2.template modify<DeviceType>();
  meam_inst_kk->k_arho3.template modify<DeviceType>();
  meam_inst_kk->k_arho3b.template modify<DeviceType>();
  meam_inst_kk->k_t_ave.template modify<DeviceType>();
  meam_inst_kk->k_tsq_ave.template modify<DeviceType>();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.template modify<DeviceType>();
    meam_inst_kk->k_arho1m.template modify<DeviceType>();
    meam_inst_kk->k_arho2m.template modify<DeviceType>();
    meam_inst_kk->k_arho3m.template modify<DeviceType>();
    meam_inst_kk->k_arho3mb.template modify<DeviceType>();
  }

  comm->reverse_comm(this);

  meam_inst_kk->k_rho0.template sync<DeviceType>();
  meam_inst_kk->k_arho2b.template sync<DeviceType>();
  meam_inst_kk->k_arho1.template sync<DeviceType>();
  meam_inst_kk->k_arho2.template sync<DeviceType>();
  meam_inst_kk->k_arho3.template sync<DeviceType>();
  meam_inst_kk->k_arho3b.template sync<DeviceType>();
  meam_inst_kk->k_t_ave.template sync<DeviceType>();
  meam_inst_kk->k_tsq_ave.template sync<DeviceType>();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.template sync<DeviceType>();
    meam_inst_kk->k_arho1m.template sync<DeviceType>();
    meam_inst_kk->k_arho2m.template sync<DeviceType>();
    meam_inst_kk->k_arho3m.template sync<DeviceType>();
    meam_inst_kk->k_arho3mb.template sync<DeviceType>();
  }

  meam_inst_kk->meam_dens_final(nlocal,eflag_either,eflag_global,eflag_atom,
                   d_eatom,ntype,type,d_map,d_scale,errorflag,ev);

  if (errorflag)
    error->one(FLERR,"MEAM library error {}",errorflag);

  meam_inst_kk->k_rho0.template modify<DeviceType>();
  meam_inst_kk->k_rho1.template modify<DeviceType>();
  meam_inst_kk->k_rho2.template modify<DeviceType>();
  meam_inst_kk->k_rho3.template modify<DeviceType>();
  meam_inst_kk->k_frhop.template modify<DeviceType>();
  meam_inst_kk->k_gamma.template modify<DeviceType>();
  meam_inst_kk->k_dgamma1.template modify<DeviceType>();
  meam_inst_kk->k_dgamma2.template modify<DeviceType>();
  meam_inst_kk->k_dgamma3.template modify<DeviceType>();
  meam_inst_kk->k_arho2b.template modify<DeviceType>();
  meam_inst_kk->k_arho1.template modify<DeviceType>();
  meam_inst_kk->k_arho2.template modify<DeviceType>();
  meam_inst_kk->k_arho3.template modify<DeviceType>();
  meam_inst_kk->k_arho3b.template modify<DeviceType>();
  meam_inst_kk->k_t_ave.template modify<DeviceType>();
  meam_inst_kk->k_tsq_ave.template modify<DeviceType>();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.template modify<DeviceType>();
    meam_inst_kk->k_arho1m.template modify<DeviceType>();
    meam_inst_kk->k_arho2m.template modify<DeviceType>();
    meam_inst_kk->k_arho3m.template modify<DeviceType>();
    meam_inst_kk->k_arho3mb.template modify<DeviceType>();
  }

  comm->forward_comm(this);

  meam_inst_kk->k_rho0.template sync<DeviceType>();
  meam_inst_kk->k_rho1.template sync<DeviceType>();
  meam_inst_kk->k_rho2.template sync<DeviceType>();
  meam_inst_kk->k_rho3.template sync<DeviceType>();
  meam_inst_kk->k_frhop.template sync<DeviceType>();
  meam_inst_kk->k_gamma.template sync<DeviceType>();
  meam_inst_kk->k_dgamma1.template sync<DeviceType>();
  meam_inst_kk->k_dgamma2.template sync<DeviceType>();
  meam_inst_kk->k_dgamma3.template sync<DeviceType>();
  meam_inst_kk->k_arho2b.template sync<DeviceType>();
  meam_inst_kk->k_arho1.template sync<DeviceType>();
  meam_inst_kk->k_arho2.template sync<DeviceType>();
  meam_inst_kk->k_arho3.template sync<DeviceType>();
  meam_inst_kk->k_arho3b.template sync<DeviceType>();
  meam_inst_kk->k_t_ave.template sync<DeviceType>();
  meam_inst_kk->k_tsq_ave.template sync<DeviceType>();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.template sync<DeviceType>();
    meam_inst_kk->k_arho1m.template sync<DeviceType>();
    meam_inst_kk->k_arho2m.template sync<DeviceType>();
    meam_inst_kk->k_arho3m.template sync<DeviceType>();
    meam_inst_kk->k_arho3mb.template sync<DeviceType>();
  }

  meam_inst_kk->meam_force(inum_half,eflag_global,eflag_atom,vflag_global,
                           vflag_atom,d_eatom,ntype,type,d_map,x,
                           d_numneigh_half, d_numneigh_full,f,d_vatom,
                           d_ilist_half, d_offset, d_neighbors_half, d_neighbors_full,
                           neighflag, need_dup, ev);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  copymode = 0;
  meam_inst_kk->copymode = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairMEAM::coeff(narg,arg);

  // sync map and scale

  int n = atom->ntypes;
  MemKK::realloc_kokkos(d_map,"pair:map",n+1);
  MemKK::realloc_kokkos(d_scale,"pair:scale",n+1,n+1);
  auto h_map = Kokkos::create_mirror_view(d_map);
  auto h_scale = Kokkos::create_mirror_view(d_scale);

  for (int i = 1; i <= n; i++) {
    h_map[i] = map[i];
    for (int j = 1; j <= n; j++)
      h_scale(i,j) = scale[i][j];
  }

  Kokkos::deep_copy(d_map,h_map);
  Kokkos::deep_copy(d_scale,h_scale);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
template<class DeviceType>
void PairMEAMKokkos<DeviceType>::init_style()
{
  PairMEAM::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this,1);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);

  request = neighbor->find_request(this,2);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair meam/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist, int iswap_in, DAT::tdual_xfloat_1d &buf,
                                int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMPackForwardComm>(0,n),*this);
  return n*comm_forward;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  int m = i*comm_forward;
  v_buf[m++] = d_rho0[j];
  v_buf[m++] = d_rho1[j];
  v_buf[m++] = d_rho2[j];
  v_buf[m++] = d_rho3[j];
  v_buf[m++] = d_frhop[j];
  v_buf[m++] = d_gamma[j];
  v_buf[m++] = d_dgamma1[j];
  v_buf[m++] = d_dgamma2[j];
  v_buf[m++] = d_dgamma3[j];
  v_buf[m++] = d_arho2b[j];
  v_buf[m++] = d_arho1(j,0);
  v_buf[m++] = d_arho1(j,1);
  v_buf[m++] = d_arho1(j,2);
  v_buf[m++] = d_arho2(j,0);
  v_buf[m++] = d_arho2(j,1);
  v_buf[m++] = d_arho2(j,2);
  v_buf[m++] = d_arho2(j,3);
  v_buf[m++] = d_arho2(j,4);
  v_buf[m++] = d_arho2(j,5);
  for (int k = 0; k < 10; k++) v_buf[m++] = d_arho3(j,k);
  v_buf[m++] = d_arho3b(j,0);
  v_buf[m++] = d_arho3b(j,1);
  v_buf[m++] = d_arho3b(j,2);
  v_buf[m++] = d_t_ave(j,0);
  v_buf[m++] = d_t_ave(j,1);
  v_buf[m++] = d_t_ave(j,2);
  v_buf[m++] = d_tsq_ave(j,0);
  v_buf[m++] = d_tsq_ave(j,1);
  v_buf[m++] = d_tsq_ave(j,2);
  if (msmeamflag) {
    v_buf[m++] = d_arho2mb[j];
    v_buf[m++] = d_arho1m(j,0);
    v_buf[m++] = d_arho1m(j,1);
    v_buf[m++] = d_arho1m(j,2);
    v_buf[m++] = d_arho2m(j,0);
    v_buf[m++] = d_arho2m(j,1);
    v_buf[m++] = d_arho2m(j,2);
    v_buf[m++] = d_arho2m(j,3);
    v_buf[m++] = d_arho2m(j,4);
    v_buf[m++] = d_arho2m(j,5);
    for (int k = 0; k < 10; k++) v_buf[m++] = d_arho3m(j,k);
    v_buf[m++] = d_arho3mb(j,0);
    v_buf[m++] = d_arho3mb(j,1);
    v_buf[m++] = d_arho3mb(j,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMUnpackForwardComm>(0,n),*this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMUnpackForwardComm, const int &i) const{
  //int m = i*38;
  int m = i*comm_forward;

    d_rho0[i+first] = v_buf[m++];
    d_rho1[i+first] = v_buf[m++];
    d_rho2[i+first] = v_buf[m++];
    d_rho3[i+first] = v_buf[m++];
    d_frhop[i+first] = v_buf[m++];
    d_gamma[i+first] = v_buf[m++];
    d_dgamma1[i+first] = v_buf[m++];
    d_dgamma2[i+first] = v_buf[m++];
    d_dgamma3[i+first] = v_buf[m++];
    d_arho2b[i+first] = v_buf[m++];
    d_arho1(i+first,0) = v_buf[m++];
    d_arho1(i+first,1) = v_buf[m++];
    d_arho1(i+first,2) = v_buf[m++];
    d_arho2(i+first,0) = v_buf[m++];
    d_arho2(i+first,1) = v_buf[m++];
    d_arho2(i+first,2) = v_buf[m++];
    d_arho2(i+first,3) = v_buf[m++];
    d_arho2(i+first,4) = v_buf[m++];
    d_arho2(i+first,5) = v_buf[m++];
    for (int k = 0; k < 10; k++) d_arho3(i+first,k) = v_buf[m++];
    d_arho3b(i+first,0) = v_buf[m++];
    d_arho3b(i+first,1) = v_buf[m++];
    d_arho3b(i+first,2) = v_buf[m++];
    d_t_ave(i+first,0) = v_buf[m++];
    d_t_ave(i+first,1) = v_buf[m++];
    d_t_ave(i+first,2) = v_buf[m++];
    d_tsq_ave(i+first,0) = v_buf[m++];
    d_tsq_ave(i+first,1) = v_buf[m++];
    d_tsq_ave(i+first,2) = v_buf[m++];
    if (msmeamflag) {
      d_arho2mb[i+first] = v_buf[m++];
      d_arho1m(i+first,0) = v_buf[m++];
      d_arho1m(i+first,1) = v_buf[m++];
      d_arho1m(i+first,2) = v_buf[m++];
      d_arho2m(i+first,0) = v_buf[m++];
      d_arho2m(i+first,1) = v_buf[m++];
      d_arho2m(i+first,2) = v_buf[m++];
      d_arho2m(i+first,3) = v_buf[m++];
      d_arho2m(i+first,4) = v_buf[m++];
      d_arho2m(i+first,5) = v_buf[m++];
      for (int k = 0; k < 10; k++) d_arho3m(i+first,k) = v_buf[m++];
      d_arho3mb(i+first,0) = v_buf[m++];
      d_arho3mb(i+first,1) = v_buf[m++];
      d_arho3mb(i+first,2) = v_buf[m++];
    }
 }

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  meam_inst_kk->k_rho0.sync_host();
  meam_inst_kk->k_rho1.sync_host();
  meam_inst_kk->k_rho2.sync_host();
  meam_inst_kk->k_rho3.sync_host();
  meam_inst_kk->k_frhop.sync_host();
  meam_inst_kk->k_gamma.sync_host();
  meam_inst_kk->k_dgamma1.sync_host();
  meam_inst_kk->k_dgamma2.sync_host();
  meam_inst_kk->k_dgamma3.sync_host();
  meam_inst_kk->k_arho2b.sync_host();
  meam_inst_kk->k_arho1.sync_host();
  meam_inst_kk->k_arho2.sync_host();
  meam_inst_kk->k_arho3.sync_host();
  meam_inst_kk->k_arho3b.sync_host();
  meam_inst_kk->k_t_ave.sync_host();
  meam_inst_kk->k_tsq_ave.sync_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.sync_host();
    meam_inst_kk->k_arho1m.sync_host();
    meam_inst_kk->k_arho2m.sync_host();
    meam_inst_kk->k_arho3m.sync_host();
    meam_inst_kk->k_arho3mb.sync_host();
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    buf[m++] = meam_inst_kk->h_rho0[j];
    buf[m++] = meam_inst_kk->h_rho1[j];
    buf[m++] = meam_inst_kk->h_rho2[j];
    buf[m++] = meam_inst_kk->h_rho3[j];
    buf[m++] = meam_inst_kk->h_frhop[j];
    buf[m++] = meam_inst_kk->h_gamma[j];
    buf[m++] = meam_inst_kk->h_dgamma1[j];
    buf[m++] = meam_inst_kk->h_dgamma2[j];
    buf[m++] = meam_inst_kk->h_dgamma3[j];
    buf[m++] = meam_inst_kk->h_arho2b[j];
    buf[m++] = meam_inst_kk->h_arho1(j,0);
    buf[m++] = meam_inst_kk->h_arho1(j,1);
    buf[m++] = meam_inst_kk->h_arho1(j,2);
    buf[m++] = meam_inst_kk->h_arho2(j,0);
    buf[m++] = meam_inst_kk->h_arho2(j,1);
    buf[m++] = meam_inst_kk->h_arho2(j,2);
    buf[m++] = meam_inst_kk->h_arho2(j,3);
    buf[m++] = meam_inst_kk->h_arho2(j,4);
    buf[m++] = meam_inst_kk->h_arho2(j,5);
    for (int k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3(j,k);
    buf[m++] = meam_inst_kk->h_arho3b(j,0);
    buf[m++] = meam_inst_kk->h_arho3b(j,1);
    buf[m++] = meam_inst_kk->h_arho3b(j,2);
    buf[m++] = meam_inst_kk->h_t_ave(j,0);
    buf[m++] = meam_inst_kk->h_t_ave(j,1);
    buf[m++] = meam_inst_kk->h_t_ave(j,2);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,0);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,1);
    buf[m++] = meam_inst_kk->h_tsq_ave(j,2);
    if (msmeamflag) {
      buf[m++] = meam_inst_kk->h_arho2mb[j];
      buf[m++] = meam_inst_kk->h_arho1m(j,0);
      buf[m++] = meam_inst_kk->h_arho1m(j,1);
      buf[m++] = meam_inst_kk->h_arho1m(j,2);
      buf[m++] = meam_inst_kk->h_arho2m(j,0);
      buf[m++] = meam_inst_kk->h_arho2m(j,1);
      buf[m++] = meam_inst_kk->h_arho2m(j,2);
      buf[m++] = meam_inst_kk->h_arho2m(j,3);
      buf[m++] = meam_inst_kk->h_arho2m(j,4);
      buf[m++] = meam_inst_kk->h_arho2m(j,5);
      for (int k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3m(j,k);
      buf[m++] = meam_inst_kk->h_arho3mb(j,0);
      buf[m++] = meam_inst_kk->h_arho3mb(j,1);
      buf[m++] = meam_inst_kk->h_arho3mb(j,2);
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  meam_inst_kk->k_rho0.sync_host();
  meam_inst_kk->k_rho1.sync_host();
  meam_inst_kk->k_rho2.sync_host();
  meam_inst_kk->k_rho3.sync_host();
  meam_inst_kk->k_frhop.sync_host();
  meam_inst_kk->k_gamma.sync_host();
  meam_inst_kk->k_dgamma1.sync_host();
  meam_inst_kk->k_dgamma2.sync_host();
  meam_inst_kk->k_dgamma3.sync_host();
  meam_inst_kk->k_arho2b.sync_host();
  meam_inst_kk->k_arho1.sync_host();
  meam_inst_kk->k_arho2.sync_host();
  meam_inst_kk->k_arho3.sync_host();
  meam_inst_kk->k_arho3b.sync_host();
  meam_inst_kk->k_t_ave.sync_host();
  meam_inst_kk->k_tsq_ave.sync_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.sync_host();
    meam_inst_kk->k_arho1m.sync_host();
    meam_inst_kk->k_arho2m.sync_host();
    meam_inst_kk->k_arho3m.sync_host();
    meam_inst_kk->k_arho3mb.sync_host();
  }

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    meam_inst_kk->h_rho0[i] = buf[m++];
    meam_inst_kk->h_rho1[i] = buf[m++];
    meam_inst_kk->h_rho2[i] = buf[m++];
    meam_inst_kk->h_rho3[i] = buf[m++];
    meam_inst_kk->h_frhop[i] = buf[m++];
    meam_inst_kk->h_gamma[i] = buf[m++];
    meam_inst_kk->h_dgamma1[i] = buf[m++];
    meam_inst_kk->h_dgamma2[i] = buf[m++];
    meam_inst_kk->h_dgamma3[i] = buf[m++];
    meam_inst_kk->h_arho2b[i] = buf[m++];
    meam_inst_kk->h_arho1(i,0) = buf[m++];
    meam_inst_kk->h_arho1(i,1) = buf[m++];
    meam_inst_kk->h_arho1(i,2) = buf[m++];
    meam_inst_kk->h_arho2(i,0) = buf[m++];
    meam_inst_kk->h_arho2(i,1) = buf[m++];
    meam_inst_kk->h_arho2(i,2) = buf[m++];
    meam_inst_kk->h_arho2(i,3) = buf[m++];
    meam_inst_kk->h_arho2(i,4) = buf[m++];
    meam_inst_kk->h_arho2(i,5) = buf[m++];
    for (int k = 0; k < 10; k++) meam_inst_kk->h_arho3(i,k) = buf[m++];
    meam_inst_kk->h_arho3b(i,0) = buf[m++];
    meam_inst_kk->h_arho3b(i,1) = buf[m++];
    meam_inst_kk->h_arho3b(i,2) = buf[m++];
    meam_inst_kk->h_t_ave(i,0) = buf[m++];
    meam_inst_kk->h_t_ave(i,1) = buf[m++];
    meam_inst_kk->h_t_ave(i,2) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,0) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,1) = buf[m++];
    meam_inst_kk->h_tsq_ave(i,2) = buf[m++];
    if (msmeamflag) {
      meam_inst_kk->h_arho2mb[i] = buf[m++];
      meam_inst_kk->h_arho1m(i,0) = buf[m++];
      meam_inst_kk->h_arho1m(i,1) = buf[m++];
      meam_inst_kk->h_arho1m(i,2) = buf[m++];
      meam_inst_kk->h_arho2m(i,0) = buf[m++];
      meam_inst_kk->h_arho2m(i,1) = buf[m++];
      meam_inst_kk->h_arho2m(i,2) = buf[m++];
      meam_inst_kk->h_arho2m(i,3) = buf[m++];
      meam_inst_kk->h_arho2m(i,4) = buf[m++];
      meam_inst_kk->h_arho2m(i,5) = buf[m++];
      for (int k = 0; k < 10; k++) meam_inst_kk->h_arho3m(i,k) = buf[m++];
      meam_inst_kk->h_arho3mb(i,0) = buf[m++];
      meam_inst_kk->h_arho3mb(i,1) = buf[m++];
      meam_inst_kk->h_arho3mb(i,2) = buf[m++];
    }
  }

  meam_inst_kk->k_rho0.modify_host();
  meam_inst_kk->k_rho1.modify_host();
  meam_inst_kk->k_rho2.modify_host();
  meam_inst_kk->k_rho3.modify_host();
  meam_inst_kk->k_frhop.modify_host();
  meam_inst_kk->k_gamma.modify_host();
  meam_inst_kk->k_dgamma1.modify_host();
  meam_inst_kk->k_dgamma2.modify_host();
  meam_inst_kk->k_dgamma3.modify_host();
  meam_inst_kk->k_arho2b.modify_host();
  meam_inst_kk->k_arho1.modify_host();
  meam_inst_kk->k_arho2.modify_host();
  meam_inst_kk->k_arho3.modify_host();
  meam_inst_kk->k_arho3b.modify_host();
  meam_inst_kk->k_t_ave.modify_host();
  meam_inst_kk->k_tsq_ave.modify_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.modify_host();
    meam_inst_kk->k_arho1m.modify_host();
    meam_inst_kk->k_arho2m.modify_host();
    meam_inst_kk->k_arho3m.modify_host();
    meam_inst_kk->k_arho3mb.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMPackReverseComm>(0,n),*this);
  //return n*30;
  return n*comm_reverse;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMPackReverseComm, const int &i) const {
  //int m = i*30;
  int m = i*comm_reverse;

  v_buf[m++] = d_rho0[i+first];
  v_buf[m++] = d_arho2b[i+first];
  v_buf[m++] = d_arho1(i+first,0);
  v_buf[m++] = d_arho1(i+first,1);
  v_buf[m++] = d_arho1(i+first,2);
  v_buf[m++] = d_arho2(i+first,0);
  v_buf[m++] = d_arho2(i+first,1);
  v_buf[m++] = d_arho2(i+first,2);
  v_buf[m++] = d_arho2(i+first,3);
  v_buf[m++] = d_arho2(i+first,4);
  v_buf[m++] = d_arho2(i+first,5);
  for (int k = 0; k < 10; k++) v_buf[m++] = d_arho3(i+first,k);
  v_buf[m++] = d_arho3b(i+first,0);
  v_buf[m++] = d_arho3b(i+first,1);
  v_buf[m++] = d_arho3b(i+first,2);
  v_buf[m++] = d_t_ave(i+first,0);
  v_buf[m++] = d_t_ave(i+first,1);
  v_buf[m++] = d_t_ave(i+first,2);
  v_buf[m++] = d_tsq_ave(i+first,0);
  v_buf[m++] = d_tsq_ave(i+first,1);
  v_buf[m++] = d_tsq_ave(i+first,2);
  if (msmeamflag) {
    v_buf[m++] = d_arho2mb[i+first];
    v_buf[m++] = d_arho1m(i+first,0);
    v_buf[m++] = d_arho1m(i+first,1);
    v_buf[m++] = d_arho1m(i+first,2);
    v_buf[m++] = d_arho2m(i+first,0);
    v_buf[m++] = d_arho2m(i+first,1);
    v_buf[m++] = d_arho2m(i+first,2);
    v_buf[m++] = d_arho2m(i+first,3);
    v_buf[m++] = d_arho2m(i+first,4);
    v_buf[m++] = d_arho2m(i+first,5);
    for (int k = 0; k < 10; k++) v_buf[m++] = d_arho3m(i+first,k);
    v_buf[m++] = d_arho3mb(i+first,0);
    v_buf[m++] = d_arho3mb(i+first,1);
    v_buf[m++] = d_arho3mb(i+first,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMEAMKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  meam_inst_kk->k_rho0.sync_host();
  meam_inst_kk->k_arho2b.sync_host();
  meam_inst_kk->k_arho1.sync_host();
  meam_inst_kk->k_arho2.sync_host();
  meam_inst_kk->k_arho3.sync_host();
  meam_inst_kk->k_arho3b.sync_host();
  meam_inst_kk->k_t_ave.sync_host();
  meam_inst_kk->k_tsq_ave.sync_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.sync_host();
    meam_inst_kk->k_arho1m.sync_host();
    meam_inst_kk->k_arho2m.sync_host();
    meam_inst_kk->k_arho3m.sync_host();
    meam_inst_kk->k_arho3mb.sync_host();
  }

  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    buf[m++] = meam_inst_kk->h_rho0[i];
    buf[m++] = meam_inst_kk->h_arho2b[i];
    buf[m++] = meam_inst_kk->h_arho1(i,0);
    buf[m++] = meam_inst_kk->h_arho1(i,1);
    buf[m++] = meam_inst_kk->h_arho1(i,2);
    buf[m++] = meam_inst_kk->h_arho2(i,0);
    buf[m++] = meam_inst_kk->h_arho2(i,1);
    buf[m++] = meam_inst_kk->h_arho2(i,2);
    buf[m++] = meam_inst_kk->h_arho2(i,3);
    buf[m++] = meam_inst_kk->h_arho2(i,4);
    buf[m++] = meam_inst_kk->h_arho2(i,5);
    for (int k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3(i,k);
    buf[m++] = meam_inst_kk->h_arho3b(i,0);
    buf[m++] = meam_inst_kk->h_arho3b(i,1);
    buf[m++] = meam_inst_kk->h_arho3b(i,2);
    buf[m++] = meam_inst_kk->h_t_ave(i,0);
    buf[m++] = meam_inst_kk->h_t_ave(i,1);
    buf[m++] = meam_inst_kk->h_t_ave(i,2);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,0);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,1);
    buf[m++] = meam_inst_kk->h_tsq_ave(i,2);
    if (msmeamflag) {
      buf[m++] = meam_inst_kk->h_arho2mb[i];
      buf[m++] = meam_inst_kk->h_arho1m(i,0);
      buf[m++] = meam_inst_kk->h_arho1m(i,1);
      buf[m++] = meam_inst_kk->h_arho1m(i,2);
      buf[m++] = meam_inst_kk->h_arho2m(i,0);
      buf[m++] = meam_inst_kk->h_arho2m(i,1);
      buf[m++] = meam_inst_kk->h_arho2m(i,2);
      buf[m++] = meam_inst_kk->h_arho2m(i,3);
      buf[m++] = meam_inst_kk->h_arho2m(i,4);
      buf[m++] = meam_inst_kk->h_arho2m(i,5);
      for (int k = 0; k < 10; k++) buf[m++] = meam_inst_kk->h_arho3m(i,k);
      buf[m++] = meam_inst_kk->h_arho3mb(i,0);
      buf[m++] = meam_inst_kk->h_arho3mb(i,1);
      buf[m++] = meam_inst_kk->h_arho3mb(i,2);
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist, int iswap_in, DAT::tdual_xfloat_1d &buf)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMEAMUnpackReverseComm>(0,n),*this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMUnpackReverseComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  //int m = i*30;
  int m = i*comm_reverse;

  d_rho0[j] += v_buf[m++];
  d_arho2b[j] += v_buf[m++];
  d_arho1(j,0) += v_buf[m++];
  d_arho1(j,1) += v_buf[m++];
  d_arho1(j,2) += v_buf[m++];
  d_arho2(j,0) += v_buf[m++];
  d_arho2(j,1) += v_buf[m++];
  d_arho2(j,2) += v_buf[m++];
  d_arho2(j,3) += v_buf[m++];
  d_arho2(j,4) += v_buf[m++];
  d_arho2(j,5) += v_buf[m++];
  for (int k = 0; k < 10; k++) d_arho3(j,k) += v_buf[m++];
  d_arho3b(j,0) += v_buf[m++];
  d_arho3b(j,1) += v_buf[m++];
  d_arho3b(j,2) += v_buf[m++];
  d_t_ave(j,0) += v_buf[m++];
  d_t_ave(j,1) += v_buf[m++];
  d_t_ave(j,2) += v_buf[m++];
  d_tsq_ave(j,0) += v_buf[m++];
  d_tsq_ave(j,1) += v_buf[m++];
  d_tsq_ave(j,2) += v_buf[m++];
  if (msmeamflag) {
    d_arho2mb[j] += v_buf[m++];
    d_arho1m(j,0) += v_buf[m++];
    d_arho1m(j,1) += v_buf[m++];
    d_arho1m(j,2) += v_buf[m++];
    d_arho2m(j,0) += v_buf[m++];
    d_arho2m(j,1) += v_buf[m++];
    d_arho2m(j,2) += v_buf[m++];
    d_arho2m(j,3) += v_buf[m++];
    d_arho2m(j,4) += v_buf[m++];
    d_arho2m(j,5) += v_buf[m++];
    for (int k = 0; k < 10; k++) d_arho3m(j,k) += v_buf[m++];
    d_arho3mb(j,0) += v_buf[m++];
    d_arho3mb(j,1) += v_buf[m++];
    d_arho3mb(j,2) += v_buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  meam_inst_kk->k_rho0.sync_host();
  meam_inst_kk->k_arho2b.sync_host();
  meam_inst_kk->k_arho1.sync_host();
  meam_inst_kk->k_arho2.sync_host();
  meam_inst_kk->k_arho3.sync_host();
  meam_inst_kk->k_arho3b.sync_host();
  meam_inst_kk->k_t_ave.sync_host();
  meam_inst_kk->k_tsq_ave.sync_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.sync_host();
    meam_inst_kk->k_arho1m.sync_host();
    meam_inst_kk->k_arho2m.sync_host();
    meam_inst_kk->k_arho3m.sync_host();
    meam_inst_kk->k_arho3mb.sync_host();
  }

  int m = 0;
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    meam_inst_kk->h_rho0[j] += buf[m++];
    meam_inst_kk->h_arho2b[j] += buf[m++];
    meam_inst_kk->h_arho1(j,0) += buf[m++];
    meam_inst_kk->h_arho1(j,1) += buf[m++];
    meam_inst_kk->h_arho1(j,2) += buf[m++];
    meam_inst_kk->h_arho2(j,0) += buf[m++];
    meam_inst_kk->h_arho2(j,1) += buf[m++];
    meam_inst_kk->h_arho2(j,2) += buf[m++];
    meam_inst_kk->h_arho2(j,3) += buf[m++];
    meam_inst_kk->h_arho2(j,4) += buf[m++];
    meam_inst_kk->h_arho2(j,5) += buf[m++];
    for (int k = 0; k < 10; k++) meam_inst_kk->h_arho3(j,k) += buf[m++];
    meam_inst_kk->h_arho3b(j,0) += buf[m++];
    meam_inst_kk->h_arho3b(j,1) += buf[m++];
    meam_inst_kk->h_arho3b(j,2) += buf[m++];
    meam_inst_kk->h_t_ave(j,0) += buf[m++];
    meam_inst_kk->h_t_ave(j,1) += buf[m++];
    meam_inst_kk->h_t_ave(j,2) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,0) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,1) += buf[m++];
    meam_inst_kk->h_tsq_ave(j,2) += buf[m++];
    if (msmeamflag) {
      meam_inst_kk->h_arho2mb[j] += buf[m++];
      meam_inst_kk->h_arho1m(j,0) += buf[m++];
      meam_inst_kk->h_arho1m(j,1) += buf[m++];
      meam_inst_kk->h_arho1m(j,2) += buf[m++];
      meam_inst_kk->h_arho2m(j,0) += buf[m++];
      meam_inst_kk->h_arho2m(j,1) += buf[m++];
      meam_inst_kk->h_arho2m(j,2) += buf[m++];
      meam_inst_kk->h_arho2m(j,3) += buf[m++];
      meam_inst_kk->h_arho2m(j,4) += buf[m++];
      meam_inst_kk->h_arho2m(j,5) += buf[m++];
      for (int k = 0; k < 10; k++) meam_inst_kk->h_arho3m(j,k) += buf[m++];
      meam_inst_kk->h_arho3mb(j,0) += buf[m++];
      meam_inst_kk->h_arho3mb(j,1) += buf[m++];
      meam_inst_kk->h_arho3mb(j,2) += buf[m++];
    }
  }

  meam_inst_kk->k_rho0.modify_host();
  meam_inst_kk->k_arho2b.modify_host();
  meam_inst_kk->k_arho1.modify_host();
  meam_inst_kk->k_arho2.modify_host();
  meam_inst_kk->k_arho3.modify_host();
  meam_inst_kk->k_arho3b.modify_host();
  meam_inst_kk->k_t_ave.modify_host();
  meam_inst_kk->k_tsq_ave.modify_host();
  if (msmeamflag) {
    meam_inst_kk->k_arho2mb.modify_host();
    meam_inst_kk->k_arho1m.modify_host();
    meam_inst_kk->k_arho2m.modify_host();
    meam_inst_kk->k_arho3m.modify_host();
    meam_inst_kk->k_arho3mb.modify_host();
  }
}

/* ----------------------------------------------------------------------
   strip special bond flags from neighbor list entries
   are not used with MEAM
   need to do here so Fortran lib doesn't see them
   done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMNeighStrip, const int &ii) const {

  const int i = d_ilist_half[ii];
  const int jnum_half = d_numneigh_half[i];
  const int jnum_full = d_numneigh_full[i];
  for (int jj = 0; jj < jnum_half; jj++)
    d_neighbors_half(i,jj) &= NEIGHMASK;
  for (int jj = 0; jj < jnum_full; jj++)
    d_neighbors_full(i,jj) &= NEIGHMASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMEAMKokkos<DeviceType>::operator()(TagPairMEAMOffsets, const int ii, int &n) const {
  const int i = d_ilist_half[ii];
  n += d_numneigh_half[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMEAMKokkos<DeviceType>::update_meam_views()
{
  d_rho0 = meam_inst_kk->d_rho0;
  d_rho1 = meam_inst_kk->d_rho1;
  d_rho2 = meam_inst_kk->d_rho2;
  d_rho3 = meam_inst_kk->d_rho3;
  d_frhop = meam_inst_kk->d_frhop;
  d_gamma = meam_inst_kk->d_gamma;
  d_dgamma1 = meam_inst_kk->d_dgamma1;
  d_dgamma2 = meam_inst_kk->d_dgamma2;
  d_dgamma3 = meam_inst_kk->d_dgamma3;
  d_arho1 = meam_inst_kk->d_arho1;
  d_arho2 = meam_inst_kk->d_arho2;
  d_arho3 = meam_inst_kk->d_arho3;
  d_arho2b = meam_inst_kk->d_arho2b;
  d_arho3b = meam_inst_kk->d_arho3b;
  d_t_ave = meam_inst_kk->d_t_ave;
  d_tsq_ave = meam_inst_kk->d_tsq_ave;
  // msmeam
  d_arho1m = meam_inst_kk->d_arho1m;
  d_arho2m = meam_inst_kk->d_arho2m;
  d_arho3m = meam_inst_kk->d_arho3m;
  d_arho2mb = meam_inst_kk->d_arho2mb;
  d_arho3mb = meam_inst_kk->d_arho3mb;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMEAMKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairMEAMKokkos<LMPHostType>;
#endif
}

