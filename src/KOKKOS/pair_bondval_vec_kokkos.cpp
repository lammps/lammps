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

#include "pair_bondval_vec_kokkos.h"

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
PairBondValVecKokkos<DeviceType>::PairBondValVecKokkos(LAMMPS *lmp) : PairBondValVec(lmp)
{
  respa_enable = 0;
  single_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datom = 3;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBondValVecKokkos<DeviceType>::~PairBondValVecKokkos()
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
void PairBondValVecKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    k_s0 = DAT::tdual_kkacc_1d_3("pair:s0",nmax);
    k_Di = DAT::tdual_kkfloat_1d_3("pair:Di",nmax);
    d_s0 = k_s0.template view<DeviceType>();
    d_Di = k_Di.template view<DeviceType>();
    h_s0 = k_s0.view_host();
    h_Di = k_Di.view_host();
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
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecInitialize>(0,nall),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecInitialize>(0,nlocal),*this);

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  // compute kernel A

  if (neighflag == HALF || neighflag == HALFTHREAD) {

    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelA<HALF,1>>(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelA<HALF,0>>(0,inum),*this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelA<HALFTHREAD,1>>(0,inum),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelA<HALFTHREAD,0>>(0,inum),*this);
      }
    }

    if (need_dup)
      Kokkos::Experimental::contribute(d_s0, dup_s0);

    // communicate and sum bond valence vectors (on the host)

    if (newton_pair) {
      k_s0.template modify<DeviceType>();
      comm->reverse_comm(this);
      k_s0.template sync<DeviceType>();
    }

    // compute kernel B

    if (eflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelB<1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecKernelB<0>>(0,inum),*this);

  } else if (neighflag == FULL) {

    // compute kernel AB

    if (eflag)
      Kokkos::parallel_reduce(
             Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelAB<1>>(0,inum),
             *this,ev);
    else
      Kokkos::parallel_for(
            policyInstance<TagPairBondValVecKernelAB<0>>::get(inum),
            *this);
  }

  if (eflag) {
    eng_vdwl += ev.evdwl;
    ev.evdwl = 0.0;
  }

  // communicate derivative of embedding function

  k_Di.template modify<DeviceType>();
  comm->forward_comm(this);
  k_Di.template sync<DeviceType>();

  // compute kernel C

  if (evflag) {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<HALF,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<HALF,0,1>>(0,inum),
              *this,ev);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<HALFTHREAD,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<HALFTHREAD,0,1>>(0,inum),
              *this,ev);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<FULL,1,1>>(0,inum),
              *this,ev);
      } else {
        Kokkos::parallel_reduce(
              Kokkos::RangePolicy<DeviceType,TagPairBondValVecKernelC<FULL,0,1>>(0,inum),
              *this,ev);
      }
    }
  } else {
    if (neighflag == HALF) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<HALF,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<HALF,0,0>>::get(inum),
              *this);
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<HALFTHREAD,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<HALFTHREAD,0,0>>::get(inum),
              *this);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<FULL,1,0>>::get(inum),
              *this);
      } else {
        Kokkos::parallel_for(
              policyInstance<TagPairBondValVecKernelC<FULL,0,0>>::get(inum),
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
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
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
void PairBondValVecKokkos<DeviceType>::allocate()
{
  PairBondValVec::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_bvv**,Kokkos::LayoutRight,DeviceType>("PairBondValVec::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBondValVecKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 3) error->all(FLERR,"Illegal pair_style command");

  PairBondValVec::settings(2,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBondValVecKokkos<DeviceType>::init_style()
{
  PairBondValVec::init_style();

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
double PairBondValVecKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBondValVec::init_one(i,j);

  k_params.view_host()(i,j).r0 = r0[i][j];
  k_params.view_host()(i,j).alpha = alpha[i][j];
  k_params.view_host()(i,j).bvvsparam = bvvsparam[i][j];
  k_params.view_host()(i,j).bvvv0 = bvvv0[i][j];
  k_params.view_host()(i,j).offset = offset[i][j];
  k_params.view_host()(i,j).cutsq = cutone*cutone;
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.template modify<LMPHostType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBondValVecKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_double_1d &buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  int buffer_size = n * datom;
  if (buf.extent(0) < (size_t) buffer_size){
      buf.resize(buffer_size);
  }
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecPackForwardComm>(0,n),*this);
  return buffer_size;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecPackForwardComm, const int &i) const {
  int j = d_sendlist(i);
  v_buf[i*3+0] = d_Di(j,0);
  v_buf[i*3+1] = d_Di(j,1);
  v_buf[i*3+2] = d_Di(j,2);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBondValVecKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_double_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairBondValVecUnpackForwardComm>(0,n),*this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecUnpackForwardComm, const int &i) const {
  d_Di(i + first,0) = v_buf[i*3+0];
  d_Di(i + first,1) = v_buf[i*3+1];
  d_Di(i + first,2) = v_buf[i*3+2];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBondValVecKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                 int /*pbc_flag*/, int * /*pbc*/)
{
  k_Di.sync_host();
  int m = 0;

  int i,j;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = h_Di(j,0);
    buf[m++] = h_Di(j,1);
    buf[m++] = h_Di(j,2);

  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBondValVecKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  k_Di.sync_host();
  int m = 0;
  for (int i = 0; i < n; i++) {
    h_Di(i + first,0) = buf[m++];
    h_Di(i + first,1) = buf[m++];
    h_Di(i + first,2) = buf[m++];
  }

  k_Di.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairBondValVecKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  k_s0.sync_host();

  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++){
      buf[m++]= h_s0(i,0);
      buf[m++]= h_s0(i,1);
      buf[m++]= h_s0(i,2);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBondValVecKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_s0.sync_host();

  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_s0(j,0) += buf[m++];
    h_s0(j,1) += buf[m++];
    h_s0(j,2) += buf[m++];
  }

  k_s0.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecInitialize, const int &i) const {
  d_s0(i,0) = 0.0;
  d_s0(i,1) = 0.0;
  d_s0(i,2) = 0.0;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelA<NEIGHFLAG,NEWTON_PAIR>, const int &ii) const {

  // s0 = bond valence vectors at each atom
  // loop over neighbors of my atoms

  // The s0 array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_s0 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_s0),decltype(ndup_s0)>::get(dup_s0,ndup_s0);
  auto a_s0 = v_s0.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  KK_ACC_FLOAT s0xtmp = 0.0;
  KK_ACC_FLOAT s0ytmp = 0.0;
  KK_ACC_FLOAT s0ztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;


    if (rsq < (d_cutsq(itype,jtype))) {
      KK_FLOAT recip = 1.0/sqrt(rsq);

      const KK_FLOAT Aij = pow((params(itype,jtype).r0)*recip,(params(itype,jtype).alpha))*recip;
      s0xtmp += Aij * (delx);
      s0ytmp += Aij * (dely);
      s0ztmp += Aij * (delz);
      if (NEWTON_PAIR || j < nlocal) {
      a_s0(j,0) -= Aij * (delx);
      a_s0(j,1) -= Aij * (dely);
      a_s0(j,2) -= Aij * (delz);
      }
    }
  }
  a_s0(i,0) += s0xtmp;
  a_s0(i,1) += s0ytmp;
  a_s0(i,2) += s0ztmp;
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // Di = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  const int i = d_ilist[ii];
  const int itype = type(i);

  KK_FLOAT s = d_s0(i,0)*d_s0(i,0)+d_s0(i,1)*d_s0(i,1)+d_s0(i,2)*d_s0(i,2)-(params(itype,itype).bvvv0)*(params(itype,itype).bvvv0);
  KK_FLOAT ss = s*s;
  d_Di(i,0) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,0)*s;
  d_Di(i,1) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,1)*s;
  d_Di(i,2) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,2)*s;
  //printf("i: %d, d_Di(i,0): %f, d_Di(i,1): %f, d_Di(i,2): %f\n", i, d_Di(i,0), d_Di(i,1), d_Di(i,2));

  if (EFLAG) {
    KK_FLOAT phi = (params(itype,itype).bvvsparam)*ss;
    //printf("i: %d, phi: %f", i, phi);
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }
}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairBondValVecKernelB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelAB<EFLAG>, const int &ii, EV_FLOAT& ev) const {

  // s0 = bond valence vectors at each atom
  // loop over neighbors of my atoms

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  KK_ACC_FLOAT s0xtmp = 0.0;
  KK_ACC_FLOAT s0ytmp = 0.0;
  KK_ACC_FLOAT s0ztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < (d_cutsq(itype,jtype))) {
      KK_FLOAT recip = 1.0/sqrt(rsq);
      const KK_FLOAT Aij = pow((params(itype,jtype).r0)*recip,(params(itype,jtype).alpha))*recip;
      s0xtmp += Aij * (delx);
      s0ytmp += Aij * (dely);
      s0ztmp += Aij * (delz);
    }
  }
  d_s0(i,0) += s0xtmp;
  d_s0(i,1) += s0ytmp;
  d_s0(i,2) += s0ztmp;

  KK_FLOAT s = d_s0(i,0)*d_s0(i,0)+d_s0(i,1)*d_s0(i,1)+d_s0(i,2)*d_s0(i,2)-(params(itype,itype).bvvv0)*(params(itype,itype).bvvv0);
  KK_FLOAT ss = s*s;
  d_Di(i,0) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,0)*s;
  d_Di(i,1) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,1)*s;
  d_Di(i,2) = (params(itype,itype).bvvsparam)*power_global*2.0*d_s0(i,2)*s;

  if (EFLAG) {
    KK_FLOAT phi = (params(itype,itype).bvvsparam)*ss;
    if (eflag_global) ev.evdwl += phi;
    if (eflag_atom) d_eatom[i] += phi;
  }

}

template<class DeviceType>
template<int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelAB<EFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<EFLAG>(TagPairBondValVecKernelAB<EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

////Specialisation for Neighborlist types Half, HalfThread, Full
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  const int jnum = d_numneigh[i];

  KK_ACC_FLOAT fxtmp = 0.0;
  KK_ACC_FLOAT fytmp = 0.0;
  KK_ACC_FLOAT fztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;


    if (rsq < (d_cutsq(itype,jtype))) {
      const KK_FLOAT r = sqrt(rsq);
      const KK_FLOAT recip = 1.0/r;
      const KK_FLOAT recip2 = recip*recip;
      const KK_FLOAT Aij = pow((params(itype,jtype).r0)*recip,(params(itype,jtype).alpha))*recip;
      const KK_FLOAT Eij = ((params(itype,jtype).alpha)+1.0)*recip2;

      KK_FLOAT fx = (d_Di(j,0)-d_Di(i,0))*Aij
           + (d_Di(i,0)-d_Di(j,0))*Eij*delx*delx*Aij
           + (d_Di(i,1)-d_Di(j,1))*Eij*delx*dely*Aij
           + (d_Di(i,2)-d_Di(j,2))*Eij*delx*delz*Aij;

      fxtmp += fx;

      KK_FLOAT fy = (d_Di(j,1)-d_Di(i,1))*Aij
           + (d_Di(i,1)-d_Di(j,1))*Eij*dely*dely*Aij
           + (d_Di(i,2)-d_Di(j,2))*Eij*dely*delz*Aij
           + (d_Di(i,0)-d_Di(j,0))*Eij*dely*delx*Aij;

      fytmp += fy;

      KK_FLOAT fz = (d_Di(j,2)-d_Di(i,2))*Aij
           + (d_Di(i,2)-d_Di(j,2))*Eij*delz*delz*Aij
           + (d_Di(i,0)-d_Di(j,0))*Eij*delz*delx*Aij
           + (d_Di(i,1)-d_Di(j,1))*Eij*delz*dely*Aij;

      fztmp += fz;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= fx;
        a_f(j,1) -= fy;
        a_f(j,2) -= fz;
      }

      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += 0.0;
        }

        if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,0.0,fx,fy,fz,delx,dely,delz);
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
void PairBondValVecKokkos<DeviceType>::operator()(TagPairBondValVecKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairBondValVecKernelC<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairBondValVecKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fx, const KK_FLOAT &fy , const KK_FLOAT &fz, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
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
      const KK_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const KK_FLOAT v0 = delx*fx;
    const KK_FLOAT v1 = dely*fy;
    const KK_FLOAT v2 = delz*fz;
    const KK_FLOAT v3 = delx*fy;
    const KK_FLOAT v4 = delx*fz;
    const KK_FLOAT v5 = dely*fz;

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
struct PairBondValVecKokkos<DeviceType>::policyInstance {

  static auto get(int inum) {
    auto policy = Kokkos::RangePolicy<DeviceType, TAG>(0,inum);
    return policy;
  }
};

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairBondValVecKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBondValVecKokkos<LMPHostType>;
#endif
}
