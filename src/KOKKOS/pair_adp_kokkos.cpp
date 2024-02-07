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
   Contributing authors: Vladislav Galigerov (HSE),  Vsevolod Nikolskiy (HSE)
------------------------------------------------------------------------- */

#include "pair_adp_kokkos.h"

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
#include "pair_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairADPKokkos<DeviceType>::PairADPKokkos(LAMMPS *lmp) : PairADP(lmp)
{
  respa_enable = 0;
  single_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairADPKokkos<DeviceType>::~PairADPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    k_rho = DAT::tdual_ffloat_1d("pair:rho",nmax);
    k_fp = DAT::tdual_ffloat_1d("pair:fp",nmax);
    k_mu = DAT::tdual_f_array("pair:mu", nmax);
    k_lambda = DAT::tdual_virial_array("pair:lambda", nmax);
    d_rho = k_rho.template view<DeviceType>();
    d_fp = k_fp.template view<DeviceType>();
    d_mu = k_mu.template view<DeviceType>();
    d_lambda = k_lambda.template view<DeviceType>();
    h_rho = k_rho.h_view;
    h_fp = k_fp.h_view;
    h_mu = k_mu.h_view;
    h_lambda = k_lambda.h_view;
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  int inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_rho    = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_rho);
    dup_mu     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_mu);
    dup_lambda = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_lambda);
    dup_f      = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_rho    = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_rho);
    ndup_mu     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_mu);
    ndup_lambda = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_lambda);
    ndup_f      = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom  = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  // zero out density

  if (newton_pair)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPInitialize>(0,nall),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPInitialize>(0,nlocal),*this);

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  // compute kernel A

  if (neighflag == HALF || neighflag == HALFTHREAD) {

    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelA<HALF,1> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelA<HALF,0> >(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelA<HALFTHREAD,1> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelA<HALFTHREAD,0> >(0,inum),*this);
      }
    }

    if (need_dup)
    {
      Kokkos::Experimental::contribute(d_rho, dup_rho);
      Kokkos::Experimental::contribute(d_mu, dup_mu);
      Kokkos::Experimental::contribute(d_lambda, dup_lambda);
    }

    // communicate and sum densities (on the host)

    if (newton_pair) {
      k_rho.template modify<DeviceType>();
      k_mu.template modify<DeviceType>();
      k_lambda.template modify<DeviceType>();
      comm->reverse_comm(this);
      k_rho.template sync<DeviceType>();
      k_mu.template sync<DeviceType>();
      k_lambda.template sync<DeviceType>();
    }

    // compute kernel B

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelB<1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelB<0> >(0,inum),*this);

  } else if (neighflag == FULL) {

    // compute kernel AB

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelAB<1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelAB<0> >(0,inum),*this);
  }

  if (eflag) {
    eng_vdwl += ev.evdwl;
    ev.evdwl = 0.0;
  }

  // communicate derivative of embedding function

  k_fp.template modify<DeviceType>();
  k_mu.template modify<DeviceType>();
  k_lambda.template modify<DeviceType>();
  comm->forward_comm(this);
  k_fp.template sync<DeviceType>();
  k_mu.template sync<DeviceType>();
  k_lambda.template sync<DeviceType>();

  // compute kernel C

  if (evflag) {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALF,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALF,0,1> >(0,inum),*this,ev);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALFTHREAD,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALFTHREAD,0,1> >(0,inum),*this,ev);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<FULL,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<FULL,0,1> >(0,inum),*this,ev);
      }
    }
  } else {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALF,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALF,0,0> >(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALFTHREAD,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<HALFTHREAD,0,0> >(0,inum),*this);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<FULL,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPKernelC<FULL,0,0> >(0,inum),*this);
      }
    }
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

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
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_rho      = decltype(dup_rho)();
    dup_mu       = decltype(dup_mu)();
    dup_lambda   = decltype(dup_lambda)();
    dup_f        = decltype(dup_f)();
    dup_eatom    = decltype(dup_eatom)();
    dup_vatom    = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::init_style()
{
  // convert read-in file(s) to arrays and spline them

  PairADP::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::file2array()
{
  PairADP::file2array();

  int i,j;
  int n = atom->ntypes;

  auto k_type2frho = DAT::tdual_int_1d("pair:type2frho",n+1);
  auto k_type2rhor = DAT::tdual_int_2d_dl("pair:type2rhor",n+1,n+1);
  auto k_type2z2r = DAT::tdual_int_2d_dl("pair:type2z2r",n+1,n+1);
  auto k_type2u2r = DAT::tdual_int_2d_dl("pair:type2u2r",n+1,n+1);
  auto k_type2w2r = DAT::tdual_int_2d_dl("pair:type2w2r",n+1,n+1);

  auto h_type2frho =  k_type2frho.h_view;
  auto h_type2rhor = k_type2rhor.h_view;
  auto h_type2z2r = k_type2z2r.h_view;
  auto h_type2u2r = k_type2u2r.h_view;
  auto h_type2w2r = k_type2w2r.h_view;

  for (i = 1; i <= n; i++) {
    h_type2frho[i] = type2frho[i];
    for (j = 1; j <= n; j++) {
      h_type2rhor(i,j) = type2rhor[i][j];
      h_type2z2r(i,j) = type2z2r[i][j];
      h_type2u2r(i,j) = type2u2r[i][j];
      h_type2w2r(i,j) = type2w2r[i][j];
    }
  }
  k_type2frho.template modify<LMPHostType>();
  k_type2frho.template sync<DeviceType>();

  k_type2rhor.template modify<LMPHostType>();
  k_type2rhor.template sync<DeviceType>();

  k_type2z2r.template modify<LMPHostType>();
  k_type2z2r.template sync<DeviceType>();

  k_type2u2r.template modify<LMPHostType>();
  k_type2u2r.template sync<DeviceType>();

  k_type2w2r.template modify<LMPHostType>();
  k_type2w2r.template sync<DeviceType>();


  d_type2frho = k_type2frho.template view<DeviceType>();
  d_type2rhor = k_type2rhor.template view<DeviceType>();
  d_type2z2r = k_type2z2r.template view<DeviceType>();
  d_type2u2r = k_type2u2r.template view<DeviceType>();
  d_type2w2r = k_type2w2r.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  tdual_ffloat_2d_n7 k_frho_spline = tdual_ffloat_2d_n7("pair:frho",nfrho,nrho+1);
  tdual_ffloat_2d_n7 k_rhor_spline = tdual_ffloat_2d_n7("pair:rhor",nrhor,nr+1);
  tdual_ffloat_2d_n7 k_z2r_spline = tdual_ffloat_2d_n7("pair:z2r",nz2r,nr+1);
  tdual_ffloat_2d_n7 k_u2r_spline = tdual_ffloat_2d_n7("pair:z2r",nu2r,nr+1);
  tdual_ffloat_2d_n7 k_w2r_spline = tdual_ffloat_2d_n7("pair:z2r",nw2r,nr+1);

  t_host_ffloat_2d_n7 h_frho_spline = k_frho_spline.h_view;
  t_host_ffloat_2d_n7 h_rhor_spline = k_rhor_spline.h_view;
  t_host_ffloat_2d_n7 h_z2r_spline = k_z2r_spline.h_view;
  t_host_ffloat_2d_n7 h_u2r_spline = k_u2r_spline.h_view;
  t_host_ffloat_2d_n7 h_w2r_spline = k_w2r_spline.h_view;

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],h_frho_spline,i);
  k_frho_spline.template modify<LMPHostType>();
  k_frho_spline.template sync<DeviceType>();

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],h_rhor_spline,i);
  k_rhor_spline.template modify<LMPHostType>();
  k_rhor_spline.template sync<DeviceType>();

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],h_z2r_spline,i);
  k_z2r_spline.template modify<LMPHostType>();
  k_z2r_spline.template sync<DeviceType>();

  for (int i = 0; i < nu2r; i++)
    interpolate(nr,dr,u2r[i],h_u2r_spline,i);
  k_u2r_spline.template modify<LMPHostType>();
  k_u2r_spline.template sync<DeviceType>();

  for (int i = 0; i < nw2r; i++)
    interpolate(nr,dr,w2r[i],h_w2r_spline,i);
  k_w2r_spline.template modify<LMPHostType>();
  k_w2r_spline.template sync<DeviceType>();

  d_frho_spline = k_frho_spline.template view<DeviceType>();
  d_rhor_spline = k_rhor_spline.template view<DeviceType>();
  d_z2r_spline = k_z2r_spline.template view<DeviceType>();
  d_u2r_spline = k_u2r_spline.template view<DeviceType>();
  d_w2r_spline = k_w2r_spline.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::interpolate(int n, double delta, double *f, t_host_ffloat_2d_n7 h_spline, int i)
{
  for (int m = 1; m <= n; m++) h_spline(i,m,6) = f[m];

  h_spline(i,1,5) = h_spline(i,2,6) - h_spline(i,1,6);
  h_spline(i,2,5) = 0.5 * (h_spline(i,3,6)-h_spline(i,1,6));
  h_spline(i,n-1,5) = 0.5 * (h_spline(i,n,6)-h_spline(i,n-2,6));
  h_spline(i,n,5) = h_spline(i,n,6) - h_spline(i,n-1,6);

  for (int m = 3; m <= n-2; m++)
    h_spline(i,m,5) = ((h_spline(i,m-2,6)-h_spline(i,m+2,6)) +
                    8.0*(h_spline(i,m+1,6)-h_spline(i,m-1,6))) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    h_spline(i,m,4) = 3.0*(h_spline(i,m+1,6)-h_spline(i,m,6)) -
      2.0*h_spline(i,m,5) - h_spline(i,m+1,5);
    h_spline(i,m,3) = h_spline(i,m,5) + h_spline(i,m+1,5) -
      2.0*(h_spline(i,m+1,6)-h_spline(i,m,6));
  }

  h_spline(i,n,4) = 0.0;
  h_spline(i,n,3) = 0.0;

  for (int m = 1; m <= n; m++) {
    h_spline(i,m,2) = h_spline(i,m,5)/delta;
    h_spline(i,m,1) = 2.0*h_spline(i,m,4)/delta;
    h_spline(i,m,0) = 3.0*h_spline(i,m,3)/delta;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairADPKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                                        int iswap_in, DAT::tdual_xfloat_1d &buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPPackForwardComm>(0,n),*this);
  return n*10;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  v_buf[10 * i] = d_fp(j);
  v_buf[10 * i + 1] = d_mu(j, 0);
  v_buf[10 * i + 2] = d_mu(j, 1);
  v_buf[10 * i + 3] = d_mu(j, 2);
  v_buf[10 * i + 4] = d_lambda(j, 0);
  v_buf[10 * i + 5] = d_lambda(j, 1);
  v_buf[10 * i + 6] = d_lambda(j, 2);
  v_buf[10 * i + 7] = d_lambda(j, 3);
  v_buf[10 * i + 8] = d_lambda(j, 4);
  v_buf[10 * i + 9] = d_lambda(j, 5);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairADPUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPUnpackForwardComm, const int &i) const {
  d_fp(i + first) = v_buf[10 * i];
  d_mu(i + first, 0) = v_buf[10 * i + 1];
  d_mu(i + first, 1) = v_buf[10 * i + 2];
  d_mu(i + first, 2) = v_buf[10 * i + 3];
  d_lambda(i + first, 0) = v_buf[10 * i + 4];
  d_lambda(i + first, 1) = v_buf[10 * i + 5];
  d_lambda(i + first, 2) = v_buf[10 * i + 6];
  d_lambda(i + first, 3) = v_buf[10 * i + 7];
  d_lambda(i + first, 4) = v_buf[10 * i + 8];
  d_lambda(i + first, 5) = v_buf[10 * i + 9];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairADPKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                 int /*pbc_flag*/, int * /*pbc*/)
{
  k_fp.sync_host();
  k_mu.sync_host();
  k_lambda.sync_host();

  int i,j, m;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    buf[m++] = h_fp(j);
    buf[m++] = h_mu(j, 0);
    buf[m++] = h_mu(j, 1);
    buf[m++] = h_mu(j, 2);
    buf[m++] = h_lambda(j, 0);
    buf[m++] = h_lambda(j, 1);
    buf[m++] = h_lambda(j, 2);
    buf[m++] = h_lambda(j, 3);
    buf[m++] = h_lambda(j, 4);
    buf[m++] = h_lambda(j, 5);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  k_fp.sync_host();
  k_mu.sync_host();
  k_lambda.sync_host();

  int m, last;
  m = 0;
  last = n + first;
  for (int i = first; i < last; i++) {
    h_fp(i) = buf[m++];
    h_mu(i, 0) = buf[m++];
    h_mu(i, 1) = buf[m++];
    h_mu(i, 2) = buf[m++];
    h_lambda(i, 0) = buf[m++];
    h_lambda(i, 1) = buf[m++];
    h_lambda(i, 2) = buf[m++];
    h_lambda(i, 3) = buf[m++];
    h_lambda(i, 4) = buf[m++];
    h_lambda(i, 5) = buf[m++];
  }

  k_fp.modify_host();
  k_mu.modify_host();
  k_lambda.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairADPKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_rho.sync_host();
  k_mu.sync_host();
  k_lambda.sync_host();

  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = h_rho(i);
    buf[m++] = h_mu(i,0);
    buf[m++] = h_mu(i,1);
    buf[m++] = h_mu(i,2);
    buf[m++] = h_lambda(i,0);
    buf[m++] = h_lambda(i,1);
    buf[m++] = h_lambda(i,2);
    buf[m++] = h_lambda(i,3);
    buf[m++] = h_lambda(i,4);
    buf[m++] = h_lambda(i,5);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairADPKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_rho.sync_host();
  k_mu.sync_host();
  k_lambda.sync_host();

  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_rho(j) += buf[m++];
    h_mu(j,0) += buf[m++];
    h_mu(j,1) += buf[m++];
    h_mu(j,2) += buf[m++];
    h_lambda(j,0) += buf[m++];
    h_lambda(j,1) += buf[m++];
    h_lambda(j,2) += buf[m++];
    h_lambda(j,3) += buf[m++];
    h_lambda(j,4) += buf[m++];
    h_lambda(j,5) += buf[m++];
  }

  k_rho.modify_host();
  k_mu.modify_host();
  k_lambda.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPInitialize, const int &i) const {
  d_rho[i] = 0.0;
  d_mu(i, 0) = 0.0;
  d_mu(i, 1) = 0.0;
  d_mu(i, 2) = 0.0;
  d_lambda(i, 0) = 0.0;
  d_lambda(i, 1) = 0.0;
  d_lambda(i, 2) = 0.0;
  d_lambda(i, 3) = 0.0;
  d_lambda(i, 4) = 0.0;
  d_lambda(i, 5) = 0.0;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelA<NEIGHFLAG,NEWTON_PAIR>, const int &ii) const {

  // rho = density at each atom
  // vector mu and tensor lambda are ADP-specific
  // loop over neighbors of my atoms

  // The rho array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_rho = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_rho),decltype(ndup_rho)>::get(dup_rho,ndup_rho);
  auto a_rho = v_rho.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_mu = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_mu),decltype(ndup_mu)>::get(dup_mu,ndup_mu);
  auto a_mu = v_mu.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_lambda = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_lambda),decltype(ndup_lambda)>::get(dup_lambda,ndup_lambda);
  auto a_lambda = v_lambda.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT rhotmp = 0.0;
  F_FLOAT mutmp[3] = {0.0,0.0,0.0};
  F_FLOAT lambdatmp[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  int d_type_ji;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cutforcesq) {
      F_FLOAT p = sqrt(rsq)*rdr + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);

      d_type_ji = d_type2rhor(jtype,itype);
      rhotmp += ((d_rhor_spline(d_type_ji,m,3)*p + d_rhor_spline(d_type_ji,m,4))*p +
                  d_rhor_spline(d_type_ji,m,5))*p + d_rhor_spline(d_type_ji,m,6);

      d_type_ji = d_type2u2r(jtype,itype);
      F_FLOAT u2 = ((d_u2r_spline(d_type_ji,m,3)*p + d_u2r_spline(d_type_ji,m,4))*p +
                  d_u2r_spline(d_type_ji,m,5))*p + d_u2r_spline(d_type_ji,m,6);
      mutmp[0] += u2*delx;
      mutmp[1] += u2*dely;
      mutmp[2] += u2*delz;

      d_type_ji = d_type2w2r(jtype,itype);
      F_FLOAT w2 = ((d_w2r_spline(d_type_ji,m,3)*p + d_w2r_spline(d_type_ji,m,4))*p +
                   d_w2r_spline(d_type_ji,m,5))*p + d_w2r_spline(d_type_ji,m,6);
      lambdatmp[0] += w2*delx*delx;
      lambdatmp[1] += w2*dely*dely;
      lambdatmp[2] += w2*delz*delz;
      lambdatmp[3] += w2*dely*delz;
      lambdatmp[4] += w2*delx*delz;
      lambdatmp[5] += w2*delx*dely;

      if (NEWTON_PAIR || j < nlocal) {
        const int d_type2rhor_ij = d_type2rhor(itype,jtype);
        a_rho[j] += ((d_rhor_spline(d_type2rhor_ij,m,3)*p + d_rhor_spline(d_type2rhor_ij,m,4))*p +
                      d_rhor_spline(d_type2rhor_ij,m,5))*p + d_rhor_spline(d_type2rhor_ij,m,6);

        const int d_type2u2r_ij = d_type2u2r(itype, jtype);
        u2 = ((d_u2r_spline(d_type2u2r_ij,m,3)*p + d_u2r_spline(d_type2u2r_ij,m,4))*p +
                      d_u2r_spline(d_type2u2r_ij,m,5))*p + d_u2r_spline(d_type2u2r_ij,m,6);
        a_mu(j,0) -= u2*delx;
        a_mu(j,1) -= u2*dely;
        a_mu(j,2) -= u2*delz;

        const int d_type2w2r_ij = d_type2w2r(itype, jtype);
        w2 = ((d_w2r_spline(d_type2w2r_ij,m,3)*p + d_w2r_spline(d_type2w2r_ij,m,4))*p +
                      d_w2r_spline(d_type2w2r_ij,m,5))*p + d_w2r_spline(d_type2w2r_ij,m,6);
        a_lambda(j, 0) += w2*delx*delx;
        a_lambda(j, 1) += w2*dely*dely;
        a_lambda(j, 2) += w2*delz*delz;
        a_lambda(j, 3) += w2*dely*delz;
        a_lambda(j, 4) += w2*delx*delz;
        a_lambda(j, 5) += w2*delx*dely;

      }
    }

  }
  a_rho[i] += rhotmp;
  a_mu(i, 0) += mutmp[0];
  a_mu(i, 1) += mutmp[1];
  a_mu(i, 2) += mutmp[2];
  a_lambda(i, 0) += lambdatmp[0];
  a_lambda(i, 1) += lambdatmp[1];
  a_lambda(i, 2) += lambdatmp[2];
  a_lambda(i, 3) += lambdatmp[3];
  a_lambda(i, 4) += lambdatmp[4];
  a_lambda(i, 5) += lambdatmp[5];
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelB<EFLAG>, const int &ii, EV_FLOAT& ev) const {
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  const int i = d_ilist[ii];
  const int itype = type(i);

  F_FLOAT p = d_rho[i]*rdrho + 1.0;
  int m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  p = MIN(p,1.0);
  const int d_type2frho_i = d_type2frho[itype];
  d_fp[i] = (d_frho_spline(d_type2frho_i,m,0)*p + d_frho_spline(d_type2frho_i,m,1))*p + d_frho_spline(d_type2frho_i,m,2);
  if (EFLAG) {
    F_FLOAT phi = ((d_frho_spline(d_type2frho_i,m,3)*p + d_frho_spline(d_type2frho_i,m,4))*p +
                    d_frho_spline(d_type2frho_i,m,5))*p + d_frho_spline(d_type2frho_i,m,6);
    phi += 0.5*(d_mu(i,0)*d_mu(i,0)+d_mu(i,1)*d_mu(i,1)+d_mu(i,2)*d_mu(i,2));
    phi += 0.5*(d_lambda(i,0)*d_lambda(i,0)+d_lambda(i,1)*
                d_lambda(i,1)+d_lambda(i,2)*d_lambda(i,2));
    phi += 1.0*(d_lambda(i,3)*d_lambda(i,3)+d_lambda(i,4)*
                d_lambda(i,4)+d_lambda(i,5)*d_lambda(i,5));
    phi -= 1.0/6.0*(d_lambda(i,0)+d_lambda(i,1)+d_lambda(i,2))*
      (d_lambda(i,0)+d_lambda(i,1)+d_lambda(i,2));
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }
}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairADPKernelB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelAB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // rho = density at each atom
  // loop over neighbors of my atoms

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT rhotmp = 0.0;
  F_FLOAT mutmp[3] = {0.0,0.0,0.0};
  F_FLOAT lambdatmp[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  int d_type_ji;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cutforcesq) {
      F_FLOAT p = sqrt(rsq)*rdr + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);
      d_type_ji = d_type2rhor(jtype,itype);
      rhotmp += ((d_rhor_spline(d_type_ji,m,3)*p + d_rhor_spline(d_type_ji,m,4))*p +
                  d_rhor_spline(d_type_ji,m,5))*p + d_rhor_spline(d_type_ji,m,6);

      d_type_ji = d_type2u2r(jtype,itype);
      F_FLOAT u2 = ((d_u2r_spline(d_type_ji,m,3)*p + d_u2r_spline(d_type_ji,m,4))*p +
                   d_u2r_spline(d_type_ji,m,5))*p + d_u2r_spline(d_type_ji,m,6);
      mutmp[0] += u2*delx;
      mutmp[1] += u2*dely;
      mutmp[2] += u2*delz;

      d_type_ji = d_type2w2r(jtype,itype);
      F_FLOAT w2 = ((d_w2r_spline(d_type_ji,m,3)*p + d_w2r_spline(d_type_ji,m,4))*p +
                   d_w2r_spline(d_type_ji,m,5))*p + d_w2r_spline(d_type_ji,m,6);
      lambdatmp[0] += w2*delx*delx;
      lambdatmp[1] += w2*dely*dely;
      lambdatmp[2] += w2*delz*delz;
      lambdatmp[3] += w2*dely*delz;
      lambdatmp[4] += w2*delx*delz;
      lambdatmp[5] += w2*delx*dely;
    }

  }
  d_rho[i] += rhotmp;

  d_mu(i, 0) += mutmp[0];
  d_mu(i, 1) += mutmp[1];
  d_mu(i, 2) += mutmp[2];

  d_lambda(i, 0) += lambdatmp[0];
  d_lambda(i, 1) += lambdatmp[1];
  d_lambda(i, 2) += lambdatmp[2];
  d_lambda(i, 3) += lambdatmp[3];
  d_lambda(i, 4) += lambdatmp[4];
  d_lambda(i, 5) += lambdatmp[5];

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  F_FLOAT p = d_rho[i]*rdrho + 1.0;
  int m = static_cast<int> (p);
  m = MAX(1,MIN(m,nrho-1));
  p -= m;
  p = MIN(p,1.0);
  const int d_type2frho_i = d_type2frho[itype];
  d_fp[i] = (d_frho_spline(d_type2frho_i,m,0)*p + d_frho_spline(d_type2frho_i,m,1))*p + d_frho_spline(d_type2frho_i,m,2);
  if (EFLAG) {
    F_FLOAT phi = ((d_frho_spline(d_type2frho_i,m,3)*p + d_frho_spline(d_type2frho_i,m,4))*p +
                    d_frho_spline(d_type2frho_i,m,5))*p + d_frho_spline(d_type2frho_i,m,6);

    phi += 0.5*(d_mu(i,0)*d_mu(i,0)+d_mu(i,1)*d_mu(i,1)+d_mu(i,2)*d_mu(i,2));
    phi += 0.5*(d_lambda(i,0)*d_lambda(i,0)+d_lambda(i,1)*
                d_lambda(i,1)+d_lambda(i,2)*d_lambda(i,2));
    phi += 1.0*(d_lambda(i,3)*d_lambda(i,3)+d_lambda(i,4)*
                d_lambda(i,4)+d_lambda(i,5)*d_lambda(i,5));
    phi -= 1.0/6.0*(d_lambda(i,0)+d_lambda(i,1)+d_lambda(i,2))*
     (d_lambda(i,0)+d_lambda(i,1)+d_lambda(i,2));
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }

}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelAB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairADPKernelAB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT fxtmp = 0.0;
  F_FLOAT fytmp = 0.0;
  F_FLOAT fztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cutforcesq) {
      const F_FLOAT r = sqrt(rsq);
      F_FLOAT p = r*rdr + 1.0;
      int m = static_cast<int> (p);
      m = MIN(m,nr-1);
      p -= m;
      p = MIN(p,1.0);

      // rhoip = derivative of (density at atom j due to atom i)
      // rhojp = derivative of (density at atom i due to atom j)
      // phi = pair potential energy
      // phip = phi'
      // z2 = phi * r
      // z2p = (phi * r)' = (phi' r) + phi
      // u2 = u
      // u2p = u'
      // w2 = w
      // w2p = w'
      // psip needs both fp[i] and fp[j] terms since r_ij appears in two
      //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
      //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

      const int d_type2rhor_ij = d_type2rhor(itype,jtype);
      const F_FLOAT rhoip = (d_rhor_spline(d_type2rhor_ij,m,0)*p + d_rhor_spline(d_type2rhor_ij,m,1))*p +
                             d_rhor_spline(d_type2rhor_ij,m,2);
      const int d_type2rhor_ji = d_type2rhor(jtype,itype);
      const F_FLOAT rhojp = (d_rhor_spline(d_type2rhor_ji,m,0)*p + d_rhor_spline(d_type2rhor_ji,m,1))*p +
                             d_rhor_spline(d_type2rhor_ji,m,2);
      const int d_type2z2r_ij = d_type2z2r(itype,jtype);
      const F_FLOAT z2p = (d_z2r_spline(d_type2z2r_ij,m,0)*p+d_z2r_spline(d_type2z2r_ij,m,1))*p + d_z2r_spline(d_type2z2r_ij,m,2);
      const F_FLOAT z2 = ((d_z2r_spline(d_type2z2r_ij,m,3)*p + d_z2r_spline(d_type2z2r_ij,m,4))*p +
          d_z2r_spline(d_type2z2r_ij,m,5))*p+d_z2r_spline(d_type2z2r_ij,m,6);

      const int d_type2u2r_ij = d_type2u2r(itype,jtype);
      const F_FLOAT u2p = (d_u2r_spline(d_type2u2r_ij,m,0)*p + d_u2r_spline(d_type2u2r_ij,m,1))*p +
                     d_u2r_spline(d_type2u2r_ij,m,2);

      const F_FLOAT u2 = ((d_u2r_spline(d_type2u2r_ij,m,3)*p + d_u2r_spline(d_type2u2r_ij,m,4))*p +
                     d_u2r_spline(d_type2u2r_ij,m,5))*p + d_u2r_spline(d_type2u2r_ij,m,6);


      const int d_type2w2r_ij = d_type2w2r(itype,jtype);
      const F_FLOAT w2p = (d_w2r_spline(d_type2w2r_ij,m,0)*p + d_w2r_spline(d_type2w2r_ij,m,1))*p +
                     d_w2r_spline(d_type2w2r_ij,m,2);

      const F_FLOAT w2 = ((d_w2r_spline(d_type2w2r_ij,m,3)*p + d_w2r_spline(d_type2w2r_ij,m,4))*p +
                     d_w2r_spline(d_type2w2r_ij,m,5))*p + d_w2r_spline(d_type2w2r_ij,m,6);



      const F_FLOAT recip = 1.0/r;
      const F_FLOAT phi = z2*recip;
      const F_FLOAT phip = z2p*recip - phi*recip;
      const F_FLOAT psip = d_fp[i]*rhojp + d_fp[j]*rhoip + phip;
      const F_FLOAT fpair = -psip*recip;

      const F_FLOAT delmux = d_mu(i, 0)-d_mu(j,0);
      const F_FLOAT delmuy = d_mu(i, 1)-d_mu(j,1);
      const F_FLOAT delmuz = d_mu(i, 2)-d_mu(j,2);
      const F_FLOAT trdelmu = delmux*delx+delmuy*dely+delmuz*delz;
      const F_FLOAT sumlamxx = d_lambda(i,0)+d_lambda(j,0);
      const F_FLOAT sumlamyy = d_lambda(i,1)+d_lambda(j,1);
      const F_FLOAT sumlamzz = d_lambda(i,2)+d_lambda(j,2);
      const F_FLOAT sumlamyz = d_lambda(i,3)+d_lambda(j,3);
      const F_FLOAT sumlamxz = d_lambda(i,4)+d_lambda(j,4);
      const F_FLOAT sumlamxy = d_lambda(i,5)+d_lambda(j,5);

      const F_FLOAT tradellam = sumlamxx*delx*delx+sumlamyy*dely*dely+
          sumlamzz*delz*delz+2.0*sumlamxy*delx*dely+
          2.0*sumlamxz*delx*delz+2.0*sumlamyz*dely*delz;
      const F_FLOAT nu = sumlamxx+sumlamyy+sumlamzz;

      const F_FLOAT adpx = -1.0*(delmux*u2 + trdelmu*u2p*delx*recip +
          2.0*w2*(sumlamxx*delx+sumlamxy*dely+sumlamxz*delz) +
          w2p*delx*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delx);
      const F_FLOAT adpy = -1.0*(delmuy*u2 + trdelmu*u2p*dely*recip +
          2.0*w2*(sumlamxy*delx+sumlamyy*dely+sumlamyz*delz) +
          w2p*dely*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*dely);
      const F_FLOAT adpz = -1.0*(delmuz*u2 + trdelmu*u2p*delz*recip +
          2.0*w2*(sumlamxz*delx+sumlamyz*dely+sumlamzz*delz) +
          w2p*delz*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delz);

      F_FLOAT fx = delx*fpair + adpx;
      F_FLOAT fy = dely*fpair + adpy;
      F_FLOAT fz = delz*fpair + adpz;

      fxtmp += fx;
      fytmp += fy;
      fztmp += fz;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= fx;
        a_f(j,1) -= fy;
        a_f(j,2) -= fz;
      }

      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*phi;
        }

        if (vflag_either || eflag_atom) this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,phi,fx,fy,fz,delx,dely,delz);
      }

    }
  }

  a_f(i,0) += fxtmp;
  a_f(i,1) += fytmp;
  a_f(i,2) += fztmp;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::operator()(TagPairADPKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairADPKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairADPKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (EFLAG) {
    if (eflag_atom) {
      const E_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const E_FLOAT v0 = delx*fx;
    const E_FLOAT v1 = dely*fy;
    const E_FLOAT v2 = delz*fz;
    const E_FLOAT v3 = delx*fy;
    const E_FLOAT v4 = delx*fz;
    const E_FLOAT v5 = dely*fz;

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
        }
      } else {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          a_vatom(i,0) += 0.5*v0;
          a_vatom(i,1) += 0.5*v1;
          a_vatom(i,2) += 0.5*v2;
          a_vatom(i,3) += 0.5*v3;
          a_vatom(i,4) += 0.5*v4;
          a_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        a_vatom(j,0) += 0.5*v0;
        a_vatom(j,1) += 0.5*v1;
        a_vatom(j,2) += 0.5*v2;
        a_vatom(j,3) += 0.5*v3;
        a_vatom(j,4) += 0.5*v4;
        a_vatom(j,5) += 0.5*v5;
        }
      } else {
        a_vatom(i,0) += 0.5*v0;
        a_vatom(i,1) += 0.5*v1;
        a_vatom(i,2) += 0.5*v2;
        a_vatom(i,3) += 0.5*v3;
        a_vatom(i,4) += 0.5*v4;
        a_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairADPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairADPKokkos<LMPHostType>;
#endif
}
