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
   Contributing authors: Peter Chun Pang Li (peterchunpang@gmail.com)
------------------------------------------------------------------------- */

#include "pair_bv_kokkos.h"
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

#define MAX_CACHE_ROWS 500

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBVKokkos<DeviceType>::PairBVKokkos(LAMMPS *lmp) : PairBV(lmp)
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
PairBVKokkos<DeviceType>::~PairBVKokkos()
{
  if (copymode) return;

  if (allocated) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
  
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

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
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  k_energy0.template sync<DeviceType>();
    
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    k_s0 = DAT::tdual_ffloat_1d("pair:s0",nmax);
    k_fp = DAT::tdual_ffloat_1d("pair:fp",nmax);
    d_s0 = k_s0.template view<DeviceType>();
    d_fp = k_fp.template view<DeviceType>();
    h_s0 = k_s0.h_view;
    h_fp = k_fp.h_view;
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
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_s0   = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_s0);
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_s0   = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_s0);
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  // zero out s0

  if (newton_pair)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVInitialize>(0,nall),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVInitialize>(0,nlocal),*this);

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  // compute kernel A

  if (neighflag == HALF || neighflag == HALFTHREAD) {

    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVKernelA<HALF,1>>(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVKernelA<HALF,0>>(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVKernelA<HALFTHREAD,1>>(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVKernelA<HALFTHREAD,0>>(0,inum),*this);
      }
    }

    if (need_dup)
      Kokkos::Experimental::contribute(d_s0, dup_s0);

    // communicate and sum bond valence (on the host)

    if (newton_pair) {
      k_s0.template modify<DeviceType>();
      comm->reverse_comm(this);
      k_s0.template sync<DeviceType>();
    }

    // compute kernel B

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBVKernelB<1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVKernelB<0>>(0,inum),*this);

  } else if (neighflag == FULL) {

    // compute kernel AB

    if (eflag)
      Kokkos::parallel_reduce(
             Kokkos::RangePolicy<DeviceType,TagPairBVKernelAB<1>>(0,inum),
             *this,ev);
    else
      Kokkos::parallel_for(
            policyInstance<TagPairBVKernelAB<0>>::get(inum),
            *this);
  }

  if (eflag) {
    eng_vdwl += ev.evdwl;
    ev.evdwl = 0.0;
  }

  // communicate derivative of embedding function

  k_fp.template modify<DeviceType>();
  comm->forward_comm(this);
  k_fp.template sync<DeviceType>();

  // compute kernel C

  if (evflag) {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<HALF,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<HALF,0,1>>(0,inum),
              *this,ev);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<HALFTHREAD,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<HALFTHREAD,0,1>>(0,inum),
              *this,ev);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<FULL,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBVKernelC<FULL,0,1>>(0,inum),
              *this,ev);
      }
    }
  } else {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<HALF,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<HALF,0,0>>::get(inum),
              *this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<HALFTHREAD,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<HALFTHREAD,0,0>>::get(inum),
              *this);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<FULL,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBVKernelC<FULL,0,0>>::get(inum),
              *this);
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
    dup_s0   = decltype(dup_s0)();
    dup_f     = decltype(dup_f)();
    dup_eatom = decltype(dup_eatom)();
    dup_vatom = decltype(dup_vatom)();
  }
}
/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::allocate()
{
  PairBV::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
    
  k_params = Kokkos::DualView<params_bv**,Kokkos::LayoutRight,DeviceType>("PairBV::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
    
  k_energy0 = DAT::tdual_ffloat_1d("pair:energy0",n+1);
  d_energy0 = k_energy0.template view<DeviceType>();

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 3) error->all(FLERR,"Illegal pair_style command");

  PairBV::settings(2,arg);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::init_style()
{
  PairBV::init_style();

  // adjust neighbor list request for KOKKOS
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairBVKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBV::init_one(i,j);

  k_params.h_view(i,j).r0 = r0[i][j];
  k_params.h_view(i,j).alpha = alpha[i][j];
  k_params.h_view(i,j).sparam = sparam[i][j];
  k_params.h_view(i,j).v0 = v0[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cutsq = cutone*cutone;
  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }
  if (i==j){
    k_energy0.h_view(i) = energy0[i];
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();
  k_energy0.template modify<LMPHostType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBVKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_xfloat_1d &buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVPackForwardComm>(0,n),*this);
  return n;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVPackForwardComm, const int &i) const {
  int j = d_sendlist(i);
  v_buf[i] = d_fp[j];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBVUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVUnpackForwardComm, const int &i) const {
  d_fp[i + first] = v_buf[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBVKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
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
void PairBVKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  k_fp.sync_host();

  for (int i = 0; i < n; i++) {
    h_fp[i + first] = buf[i];
  }

  k_fp.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBVKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_s0.sync_host();

  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = h_s0[i];
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBVKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_s0.sync_host();

  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_s0[j] += buf[m++];
  }

  k_s0.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVInitialize, const int &i) const {
  d_s0[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelA<NEIGHFLAG,NEWTON_PAIR>, const int &ii) const {

  // s0 = bond valence at each atom
  // loop over neighbors of my atoms

  // The s0 array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_s0 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_s0),decltype(ndup_s0)>::get(dup_s0,ndup_s0);
  auto a_s0 = v_s0.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT s0tmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    
    if((params(itype,jtype).alpha)!=0.0){
        if (rsq < (d_cutsq(itype,jtype))) {
          F_FLOAT recip = 1.0/sqrt(rsq);
          F_FLOAT r = sqrt(rsq);
          s0tmp += pow((params(itype,jtype).r0)/r,(params(itype,jtype).alpha));
          if (NEWTON_PAIR || j < nlocal) {
            a_s0[j] += pow((params(jtype,itype).r0)/r,(params(jtype,itype).alpha));
          }
        }
    }

  }
  a_s0[i] += s0tmp;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  const int i = d_ilist[ii];
  const int itype = type(i);

  F_FLOAT s = d_s0[i] - (params(itype,itype).v0);
  F_FLOAT ss = s*s;
  d_fp[i] = (params(itype,itype).sparam)*power_global*s;
  if (EFLAG) {
    F_FLOAT phi = (params(itype,itype).sparam)*ss+(d_energy0(itype));
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }
}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairBVKernelB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelAB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // s0 = bond valence at each atom
  // loop over neighbors of my atoms
  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  F_FLOAT s0tmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if((params(itype,jtype).alpha)!=0.0){
        if (rsq < (d_cutsq(itype,jtype))) {
          F_FLOAT recip = 1.0/sqrt(rsq);
          F_FLOAT r = sqrt(rsq);
          s0tmp += pow((params(itype,jtype).r0)/r,(params(itype,jtype).alpha));
          }
    }
  }
  
  d_s0[i] += s0tmp;

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  F_FLOAT s = d_s0[i]- (params(itype,itype).v0);
  F_FLOAT ss = s*s;
  d_fp[i] = (params(itype,itype).sparam)*power_global*s;
  if (EFLAG) {
    F_FLOAT phi = (params(itype,itype).sparam)*ss+(d_energy0(itype));
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }

}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelAB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairBVKernelAB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

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
    
    if((params(itype,jtype).alpha)!=0.0){
        if (rsq < (d_cutsq(itype,jtype))) {
          const F_FLOAT r = sqrt(rsq);
          const F_FLOAT recip = 1.0/r;
          const F_FLOAT Aij = (params(itype,jtype).alpha)*pow((params(itype,jtype).r0)*recip,(params(itype,jtype).alpha))*recip;
          const F_FLOAT psip = (d_fp[i]+d_fp[j])*Aij;
          const F_FLOAT fpair = psip*recip;

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
              ev.evdwl += 0.0;
            }

            if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,0.0,fpair,delx,dely,delz);
          }

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
void PairBVKokkos<DeviceType>::operator()(TagPairBVKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairBVKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairBVKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
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

/* ---------------------------------------------------------------------- */

template<typename DeviceType>
template<class TAG>
struct PairBVKokkos<DeviceType>::policyInstance {

  static auto get(int inum) {
    auto policy = Kokkos::RangePolicy<DeviceType, TAG>(0,inum);
    return policy;
  }
};

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairBVKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBVKokkos<LMPHostType>;
#endif
}
