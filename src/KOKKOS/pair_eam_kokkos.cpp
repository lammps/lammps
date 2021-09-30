// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Stan Moore (SNL), Christian Trott (SNL)
------------------------------------------------------------------------- */

#include "pair_eam_kokkos.h"

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
PairEAMKokkos<DeviceType>::PairEAMKokkos(LAMMPS *lmp) : PairEAM(lmp)
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
PairEAMKokkos<DeviceType>::~PairEAMKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
    d_rho = k_rho.template view<DeviceType>();
    d_fp = k_fp.template view<DeviceType>();
    h_rho = k_rho.h_view;
    h_fp = k_fp.h_view;
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
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
    dup_rho   = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_rho);
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_rho   = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_rho);
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  // zero out density

  if (newton_pair)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMInitialize>(0,nall),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMInitialize>(0,nlocal),*this);

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  // compute kernel A

  if (neighflag == HALF || neighflag == HALFTHREAD) {

    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelA<HALF,1> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelA<HALF,0> >(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelA<HALFTHREAD,1> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelA<HALFTHREAD,0> >(0,inum),*this);
      }
    }

    if (need_dup)
      Kokkos::Experimental::contribute(d_rho, dup_rho);

    // communicate and sum densities (on the host)

    if (newton_pair) {
      k_rho.template modify<DeviceType>();
      comm->reverse_comm_pair(this);
      k_rho.template sync<DeviceType>();
    }

    // compute kernel B

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelB<1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelB<0> >(0,inum),*this);

  } else if (neighflag == FULL) {

    // compute kernel AB

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelAB<1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelAB<0> >(0,inum),*this);
  }

  if (eflag) {
    eng_vdwl += ev.evdwl;
    ev.evdwl = 0.0;
  }

  // communicate derivative of embedding function

  k_fp.template modify<DeviceType>();
  comm->forward_comm_pair(this);
  k_fp.template sync<DeviceType>();

  // compute kernel C

  if (evflag) {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALF,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALF,0,1> >(0,inum),*this,ev);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALFTHREAD,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALFTHREAD,0,1> >(0,inum),*this,ev);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<FULL,1,1> >(0,inum),*this,ev);
      } else {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<FULL,0,1> >(0,inum),*this,ev);
      }
    }
  } else {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALF,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALF,0,0> >(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALFTHREAD,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<HALFTHREAD,0,0> >(0,inum),*this);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<FULL,1,0> >(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMKernelC<FULL,0,0> >(0,inum),*this);
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
    dup_rho   = decltype(dup_rho)();
    dup_f     = decltype(dup_f)();
    dup_eatom = decltype(dup_eatom)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::init_style()
{
  // convert read-in file(s) to arrays and spline them

  PairEAM::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with pair eam/kk");
  }

}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::file2array()
{
  PairEAM::file2array();

  int i,j;
  int n = atom->ntypes;

  auto k_type2frho = DAT::tdual_int_1d("pair:type2frho",n+1);
  auto k_type2rhor = DAT::tdual_int_2d_dl("pair:type2rhor",n+1,n+1);
  auto k_type2z2r = DAT::tdual_int_2d_dl("pair:type2z2r",n+1,n+1);

  auto h_type2frho =  k_type2frho.h_view;
  auto h_type2rhor = k_type2rhor.h_view;
  auto h_type2z2r = k_type2z2r.h_view;

  for (i = 1; i <= n; i++) {
    h_type2frho[i] = type2frho[i];
    for (j = 1; j <= n; j++) {
      h_type2rhor(i,j) = type2rhor[i][j];
      h_type2z2r(i,j) = type2z2r[i][j];
    }
  }
  k_type2frho.template modify<LMPHostType>();
  k_type2frho.template sync<DeviceType>();
  k_type2rhor.template modify<LMPHostType>();
  k_type2rhor.template sync<DeviceType>();
  k_type2z2r.template modify<LMPHostType>();
  k_type2z2r.template sync<DeviceType>();

  d_type2frho = k_type2frho.template view<DeviceType>();
  d_type2rhor = k_type2rhor.template view<DeviceType>();
  d_type2z2r = k_type2z2r.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  tdual_ffloat_2d_n7 k_frho_spline = tdual_ffloat_2d_n7("pair:frho",nfrho,nrho+1);
  tdual_ffloat_2d_n7 k_rhor_spline = tdual_ffloat_2d_n7("pair:rhor",nrhor,nr+1);
  tdual_ffloat_2d_n7 k_z2r_spline = tdual_ffloat_2d_n7("pair:z2r",nz2r,nr+1);

  t_host_ffloat_2d_n7 h_frho_spline = k_frho_spline.h_view;
  t_host_ffloat_2d_n7 h_rhor_spline = k_rhor_spline.h_view;
  t_host_ffloat_2d_n7 h_z2r_spline = k_z2r_spline.h_view;

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

  d_frho_spline = k_frho_spline.template view<DeviceType>();
  d_rhor_spline = k_rhor_spline.template view<DeviceType>();
  d_z2r_spline = k_z2r_spline.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::interpolate(int n, double delta, double *f, t_host_ffloat_2d_n7 h_spline, int i)
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
int PairEAMKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                                        int iswap_in, DAT::tdual_xfloat_1d &buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMPackForwardComm>(0,n),*this);
  return n;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  v_buf[i] = d_fp[j];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairEAMUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMUnpackForwardComm, const int &i) const {
  d_fp[i + first] = v_buf[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairEAMKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                 int /*pbc_flag*/, int * /*pbc*/)
{
  k_fp.sync_host();

  int i,j;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[i] = h_fp[j];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  k_fp.sync_host();

  for (int i = 0; i < n; i++) {
    h_fp[i + first] = buf[i];
  }

  k_fp.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairEAMKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_rho.sync_host();

  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = h_rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairEAMKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_rho.sync_host();

  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_rho[j] += buf[m++];
  }

  k_rho.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMInitialize, const int &i) const {
  d_rho[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelA<NEIGHFLAG,NEWTON_PAIR>, const int &ii) const {

  // rho = density at each atom
  // loop over neighbors of my atoms

  // The rho array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_rho = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_rho),decltype(ndup_rho)>::get(dup_rho,ndup_rho);
  auto a_rho = v_rho.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT rhotmp = 0.0;

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
      const int d_type2rhor_ji = d_type2rhor(jtype,itype);
      rhotmp += ((d_rhor_spline(d_type2rhor_ji,m,3)*p + d_rhor_spline(d_type2rhor_ji,m,4))*p +
                  d_rhor_spline(d_type2rhor_ji,m,5))*p + d_rhor_spline(d_type2rhor_ji,m,6);
      if (NEWTON_PAIR || j < nlocal) {
        const int d_type2rhor_ij = d_type2rhor(itype,jtype);
        a_rho[j] += ((d_rhor_spline(d_type2rhor_ij,m,3)*p + d_rhor_spline(d_type2rhor_ij,m,4))*p +
                      d_rhor_spline(d_type2rhor_ij,m,5))*p + d_rhor_spline(d_type2rhor_ij,m,6);
      }
    }

  }
  a_rho[i] += rhotmp;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelB<EFLAG>, const int &ii, EV_FLOAT& ev) const {
  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

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
    if (d_rho[i] > rhomax) phi += d_fp[i] * (d_rho[i]-rhomax);
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }
}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairEAMKernelB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelAB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // rho = density at each atom
  // loop over neighbors of my atoms

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT rhotmp = 0.0;

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
      const int d_type2rhor_ji = d_type2rhor(jtype,itype);
      rhotmp += ((d_rhor_spline(d_type2rhor_ji,m,3)*p + d_rhor_spline(d_type2rhor_ji,m,4))*p +
                  d_rhor_spline(d_type2rhor_ji,m,5))*p + d_rhor_spline(d_type2rhor_ji,m,6);
    }

  }
  d_rho[i] += rhotmp;

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

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
    if (d_rho[i] > rhomax) phi += d_fp[i] * (d_rho[i]-rhomax);
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }

}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelAB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairEAMKernelAB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

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

      const auto z2r_spline_3 = d_z2r_spline(d_type2z2r_ij,m,3);
      const auto z2r_spline_4 = d_z2r_spline(d_type2z2r_ij,m,4);
      const auto z2r_spline_5 = d_z2r_spline(d_type2z2r_ij,m,5);
      const auto z2r_spline_6 = d_z2r_spline(d_type2z2r_ij,m,6);

      const F_FLOAT z2p = (3.0*rdr*z2r_spline_3*p + 2.0*rdr*z2r_spline_4)*p +
                           rdr*z2r_spline_5; // the rdr and the factors of 3.0 and 2.0 come out of the interpolate function
      const F_FLOAT z2 = ((z2r_spline_3*p + z2r_spline_4)*p +
                           z2r_spline_5)*p + z2r_spline_6;

      const F_FLOAT recip = 1.0/r;
      const F_FLOAT phi = z2*recip;
      const F_FLOAT phip = z2p*recip - phi*recip;
      const F_FLOAT psip = d_fp[i]*rhojp + d_fp[j]*rhoip + phip;
      const F_FLOAT fpair = -psip*recip;

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= delx*fpair;
        a_f(j,1) -= dely*fpair;
        a_f(j,2) -= delz*fpair;
      }

      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*phi;
        }

        if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,phi,fpair,delx,dely,delz);
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
void PairEAMKokkos<DeviceType>::operator()(TagPairEAMKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairEAMKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairEAMKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

  auto v_vatom = ScatterViewHelper<typename NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<typename AtomicDup<NEIGHFLAG,DeviceType>::value>();

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
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

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
template class PairEAMKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairEAMKokkos<LMPHostType>;
#endif
}
