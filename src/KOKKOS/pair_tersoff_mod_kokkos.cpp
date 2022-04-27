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
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_tersoff_mod_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTersoffMODKokkos<DeviceType>::PairTersoffMODKokkos(LAMMPS *lmp) : PairTersoffMOD(lmp)
{
  respa_enable = 0;
  suffix_flag |= Suffix::KOKKOS;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTersoffMODKokkos<DeviceType>::~PairTersoffMODKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffMODKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairTersoffMOD::coeff(narg,arg);

  // sync map

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.h_view;

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.template modify<LMPHostType>();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffMODKokkos<DeviceType>::init_style()
{
  PairTersoffMOD::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair tersoff/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffMODKokkos<DeviceType>::setup_params()
{
  PairTersoffMOD::setup_params();

  // sync elem3param and params

  tdual_int_3d k_elem3param = tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  t_host_int_3d h_elem3param = k_elem3param.h_view;

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.h_view;

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem3param(i,j,k) = elem3param[i][j][k];

  for (int m = 0; m < nparams; m++)
    h_params[m] = params[m];

  k_elem3param.modify_host();
  k_elem3param.template sync<DeviceType>();
  k_params.modify_host();
  k_params.template sync<DeviceType>();

  d_elem3param = k_elem3param.template view<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffMODKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_either) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;

  inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int)d_neighbors_short.extent(1) < max_neighs) ||
     ((int)d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = Kokkos::View<int**,DeviceType>("Tersoff::neighbors_short",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum)
    d_numneigh_short = Kokkos::View<int*,DeviceType>("Tersoff::numneighs_short",ignum*1.2);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTersoffMODComputeShortNeigh>(0,inum), *this);

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffMODCompute<HALF,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffMODCompute<HALF,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffMODCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffMODCompute<HALFTHREAD,0> >(0,inum),*this);
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev_all.evdwl;
  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
  }

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_either) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_eatom = decltype(dup_eatom)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::operator()(TagPairTersoffMODComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const F_FLOAT cutmax_sq = cutmax*cutmax;

    const int jnum = d_numneigh[i];
    int inside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const X_FLOAT delx = xtmp - x(j,0);
      const X_FLOAT dely = ytmp - x(j,1);
      const X_FLOAT delz = ztmp - x(j,2);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutmax_sq) {
        d_neighbors_short(ii,inside) = j;
        inside++;
      }
    }
    d_numneigh_short(ii) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::operator()(TagPairTersoffMODCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  const auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  const auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  if (i >= nlocal) return;
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = d_map(type(i));
  const tagint itag = tag(i);

  F_FLOAT fi[3], fj[3], fk[3];

  //const AtomNeighborsConst d_neighbors_i = k_list.get_neighbors_const(i);
  const int jnum = d_numneigh_short[ii];

  // repulsive

  F_FLOAT f_x = 0.0;
  F_FLOAT f_y = 0.0;
  F_FLOAT f_z = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(ii,jj);
    const int jtype = d_map(type(j));
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const int iparam_ij = d_elem3param(itype,jtype,jtype);
    const F_FLOAT cutsq = d_params(iparam_ij).cutsq;

    if (rsq >= cutsq) continue;

    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT tmp_fce = ters_fc_k(d_params(iparam_ij),r);
    const F_FLOAT tmp_fcd = ters_dfc(d_params(iparam_ij),r);
    const F_FLOAT tmp_exp = exp(-d_params(iparam_ij).lam1 * r);
    const F_FLOAT frep = -d_params(iparam_ij).biga * tmp_exp *
                          (tmp_fcd - tmp_fce*d_params(iparam_ij).lam1) / r;
    const F_FLOAT eng = tmp_fce * d_params(iparam_ij).biga * tmp_exp;

    f_x += delx*frep;
    f_y += dely*frep;
    f_z += delz*frep;
    a_f(j,0) -= delx*frep;
    a_f(j,1) -= dely*frep;
    a_f(j,2) -= delz*frep;

    if (EVFLAG) {
      if (eflag) ev.evdwl += eng;
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,frep,delx,dely,delz);
    }
  }

  // attractive: bond order

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(ii,jj);
    const int jtype = d_map(type(j));

    const F_FLOAT delx1 = xtmp - x(j,0);
    const F_FLOAT dely1 = ytmp - x(j,1);
    const F_FLOAT delz1 = ztmp - x(j,2);
    const F_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    const int iparam_ij = d_elem3param(itype,jtype,jtype);
    const F_FLOAT cutsq1 = d_params(iparam_ij).cutsq;

    F_FLOAT bo_ij = 0.0;
    if (rsq1 > cutsq1) continue;
    const F_FLOAT rij = sqrt(rsq1);

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      int k = d_neighbors_short(ii,kk);
      const int ktype = d_map(type(k));

      const F_FLOAT delx2 = xtmp - x(k,0);
      const F_FLOAT dely2 = ytmp - x(k,1);
      const F_FLOAT delz2 = ztmp - x(k,2);
      const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      const int iparam_ijk = d_elem3param(itype,jtype,ktype);
      const F_FLOAT cutsq2 = d_params(iparam_ijk).cutsq;

      if (rsq2 > cutsq2) continue;
      const F_FLOAT rik = sqrt(rsq2);
      bo_ij += bondorder(d_params(iparam_ijk),rij,delx1,dely1,delz1,rik,delx2,dely2,delz2);
    }

    // attractive: pairwise potential and force

    const F_FLOAT fa = ters_fa_k(d_params(iparam_ij),rij);
    const F_FLOAT dfa = ters_dfa(d_params(iparam_ij),rij);
    const F_FLOAT bij = ters_bij_k(d_params(iparam_ij),bo_ij);
    const F_FLOAT fatt = -0.5*bij * dfa / rij;
    const F_FLOAT prefactor = 0.5*fa * ters_dbij(d_params(iparam_ij),bo_ij);

    f_x += delx1*fatt;
    f_y += dely1*fatt;
    f_z += delz1*fatt;
    F_FLOAT fj_x = -delx1*fatt;
    F_FLOAT fj_y = -dely1*fatt;
    F_FLOAT fj_z = -delz1*fatt;

    if (EVFLAG) {
      const F_FLOAT eng = 0.5*bij * fa;
      if (eflag) ev.evdwl += eng;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,fatt,delx1,dely1,delz1);
    }

    // attractive: three-body force

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      int k = d_neighbors_short(ii,kk);
      const int ktype = d_map(type(k));

      const F_FLOAT delx2 = xtmp - x(k,0);
      const F_FLOAT dely2 = ytmp - x(k,1);
      const F_FLOAT delz2 = ztmp - x(k,2);
      const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      const int iparam_ijk = d_elem3param(itype,jtype,ktype);
      const F_FLOAT cutsq2 = d_params(iparam_ijk).cutsq;

      if (rsq2 > cutsq2) continue;
      const F_FLOAT rik = sqrt(rsq2);
      ters_dthb(d_params(iparam_ijk),prefactor,rij,delx1,dely1,delz1,
                rik,delx2,dely2,delz2,fi,fj,fk);

      f_x += fi[0];
      f_y += fi[1];
      f_z += fi[2];
      fj_x += fj[0];
      fj_y += fj[1];
      fj_z += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (vflag_either) {
        F_FLOAT delrij[3], delrik[3];
        delrij[0] = -delx1; delrij[1] = -dely1; delrij[2] = -delz1;
        delrik[0] = -delx2; delrik[1] = -dely2; delrik[2] = -delz2;
        if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj,fk,delrij,delrik);
      }
    }

    a_f(j,0) += fj_x;
    a_f(j,1) += fj_y;
    a_f(j,2) += fj_z;
  }
  a_f(i,0) += f_x;
  a_f(i,1) += f_y;
  a_f(i,2) += f_z;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::operator()(TagPairTersoffMODCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairTersoffMODCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_fc_k(const Param& param, const F_FLOAT &r) const
{
  const F_FLOAT ters_R = param.bigr;
  const F_FLOAT ters_D = param.bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - 1.125*sin(MY_PI2*(r - ters_R)/ters_D) -
              0.125*sin(3.0*MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_dfc(const Param& param, const F_FLOAT &r) const
{
  const F_FLOAT ters_R = param.bigr;
  const F_FLOAT ters_D = param.bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(0.375*MY_PI4/ters_D) * (3.0*cos(MY_PI2*(r - ters_R)/ters_D) +
                                   cos(3.0*MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::bondorder(const Param& param,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2) const
{
  F_FLOAT arg, ex_delr;

  const F_FLOAT costheta = (dx1*dx2 + dy1*dy2 + dz1*dz2)/(rij*rik);

  const F_FLOAT paramtmp = param.lam3 * (rij-rik);
  if (int(param.powerm) == 3) arg = paramtmp*paramtmp*paramtmp;//pow(param.lam3 * (rij-rik),3.0);
  else arg = paramtmp;

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc_k(param,rik) * ters_gijk(param,costheta) * ex_delr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::
        ters_gijk(const Param& param, const F_FLOAT &cos) const
{
  const F_FLOAT ters_c1 = param.c1;
  const F_FLOAT ters_c2 = param.c2;
  const F_FLOAT ters_c3 = param.c3;
  const F_FLOAT ters_c4 = param.c4;
  const F_FLOAT ters_c5 = param.c5;
  const F_FLOAT tmp_h = (param.h - cos)*(param.h - cos);

  return ters_c1 + (ters_c2*tmp_h/(ters_c3 + tmp_h)) *
      (1.0 + ters_c4*exp(-ters_c5*tmp_h));

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::
        ters_dgijk(const Param& param, const F_FLOAT &cos) const
{
  const F_FLOAT ters_c2 = param.c2;
  const F_FLOAT ters_c3 = param.c3;
  const F_FLOAT ters_c4 = param.c4;
  const F_FLOAT ters_c5 = param.c5;
  const F_FLOAT tmp_h = (param.h - cos)*(param.h - cos);
  const F_FLOAT g1 = (param.h - cos)/(ters_c3 + tmp_h);
  const F_FLOAT g2 = exp(-ters_c5*tmp_h);

  return -2.0*ters_c2*g1*((1 + ters_c4*g2)*(1 + g1*(cos - param.h)) -
                            tmp_h*ters_c4*ters_c5*g2);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_fa_k(const Param& param, const F_FLOAT &r) const
{
  if (r > param.bigr + param.bigd) return 0.0;
  return -param.bigb * exp(-param.lam2 * r)
          * ters_fc_k(param,r);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_dfa(const Param& param, const F_FLOAT &r) const
{
  if (r > param.bigr + param.bigd) return 0.0;
  return param.bigb * exp(-param.lam2 * r) *
    (param.lam2 * ters_fc_k(param,r) - ters_dfc(param,r));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_bij_k(const Param& param, const F_FLOAT &bo) const
{
  const F_FLOAT tmp = param.beta * bo;
  if (tmp > param.ca1)
    return pow(tmp, -param.powern/(2.0*param.powern_del));
  if (tmp < param.ca4)
    return 1.0;
  return pow(1.0 + pow(tmp,param.powern), -1.0/(2.0*param.powern_del));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffMODKokkos<DeviceType>::ters_dbij(const Param& param, const F_FLOAT &bo) const
{
  const F_FLOAT tmp = param.beta * bo;
  if (tmp > param.ca1)
    return -0.5*(param.powern/param.powern_del)*
          pow(tmp,-0.5*(param.powern/param.powern_del)) / bo;
  if (tmp < param.ca4)
    return 0.0;

  const F_FLOAT tmp_n = pow(tmp,param.powern);
  return -0.5 *(param.powern/param.powern_del)*
          pow(1.0+tmp_n, -1.0-(1.0/(2.0*param.powern_del)))*tmp_n / bo;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::ters_dthb(
        const Param& param, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk) const
{
  // from PairTersoffMOD::attractive
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  //rij = sqrt(rsq1);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  //rik = sqrt(rsq2);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  // from PairTersoffMOD::ters_zetaterm_d
  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(param,rik);
  dfc = ters_dfc(param,rik);

  const F_FLOAT paramtmp = param.lam3 * (rij-rik);
  if (int(param.powerm) == 3) tmp = paramtmp*paramtmp*paramtmp;//pow(param.lam3 * (rij-rik),3.0);
  else tmp = paramtmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(param.powerm) == 3)
    dex_delr = 3.0*paramtmp*paramtmp*param.lam3*ex_delr;//pow(rij-rik,2.0)*ex_delr;
  else dex_delr = param.lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(param,cos);
  dgijk = ters_dgijk(param,cos);

  // from PairTersoffMOD::costheta_d
  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,fi);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfi,fi,fi);
  vec3_scaleadd(fc*gijk*dex_delr,rik_hat,fi,fi);
  vec3_scaleadd(-fc*gijk*dex_delr,rij_hat,fi,fi);
  vec3_scale(prefactor,fi,fi);

  vec3_scale(fc*dgijk*ex_delr,dcosfj,fj);
  vec3_scaleadd(fc*gijk*dex_delr,rij_hat,fj,fj);
  vec3_scale(prefactor,fj,fj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::ters_dthbj(
        const Param& param, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fj, F_FLOAT *fk) const
{
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(param,rik);
  dfc = ters_dfc(param,rik);
  const F_FLOAT paramtmp = param.lam3 * (rij-rik);
  if (int(param.powerm) == 3) tmp = paramtmp*paramtmp*paramtmp;//pow(param.lam3 * (rij-rik),3.0);
  else tmp = paramtmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(param.powerm) == 3)
    dex_delr = 3.0*paramtmp*paramtmp*param.lam3*ex_delr;//pow(param.lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else dex_delr = param.lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(param,cos);
  dgijk = ters_dgijk(param,cos);

  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(fc*dgijk*ex_delr,dcosfj,fj);
  vec3_scaleadd(fc*gijk*dex_delr,rij_hat,fj,fj);
  vec3_scale(prefactor,fj,fj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::ters_dthbk(
        const Param& param, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fk) const
{
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(param,rik);
  dfc = ters_dfc(param,rik);
  const F_FLOAT paramtmp = param.lam3 * (rij-rik);
  if (int(param.powerm) == 3) tmp = paramtmp*paramtmp*paramtmp;//pow(param.lam3 * (rij-rik),3.0);
  else tmp = paramtmp;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(param.powerm) == 3)
    dex_delr = 3.0*paramtmp*paramtmp*param.lam3*ex_delr;//pow(param.lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else dex_delr = param.lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(param,cos);
  dgijk = ters_dgijk(param,cos);

  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    a_eatom[i] += epairhalf;
    a_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      ev.v[0] += v0;
      ev.v[1] += v1;
      ev.v[2] += v2;
      ev.v[3] += v3;
      ev.v[4] += v4;
      ev.v[5] += v5;
    }

    if (vflag_atom) {
      a_vatom(i,0) += 0.5*v0;
      a_vatom(i,1) += 0.5*v1;
      a_vatom(i,2) += 0.5*v2;
      a_vatom(i,3) += 0.5*v3;
      a_vatom(i,4) += 0.5*v4;
      a_vatom(i,5) += 0.5*v5;

      a_vatom(j,0) += 0.5*v0;
      a_vatom(j,1) += 0.5*v1;
      a_vatom(j,2) += 0.5*v2;
      a_vatom(j,3) += 0.5*v3;
      a_vatom(j,4) += 0.5*v4;
      a_vatom(j,5) += 0.5*v5;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
        F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT v[6];

  v[0] = (drij[0]*fj[0] + drik[0]*fk[0]);
  v[1] = (drij[1]*fj[1] + drik[1]*fk[1]);
  v[2] = (drij[2]*fj[2] + drik[2]*fk[2]);
  v[3] = (drij[0]*fj[1] + drik[0]*fk[1]);
  v[4] = (drij[0]*fj[2] + drik[0]*fk[2]);
  v[5] = (drij[1]*fj[2] + drik[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    v[0] *= THIRD;
    v[1] *= THIRD;
    v[2] *= THIRD;
    v[3] *= THIRD;
    v[4] *= THIRD;
    v[5] *= THIRD;

    a_vatom(i,0) += v[0]; a_vatom(i,1) += v[1]; a_vatom(i,2) += v[2];
    a_vatom(i,3) += v[3]; a_vatom(i,4) += v[4]; a_vatom(i,5) += v[5];

    a_vatom(j,0) += v[0]; a_vatom(j,1) += v[1]; a_vatom(j,2) += v[2];
    a_vatom(j,3) += v[3]; a_vatom(j,4) += v[4]; a_vatom(j,5) += v[5];

    a_vatom(k,0) += v[0]; a_vatom(k,1) += v[1]; a_vatom(k,2) += v[2];
    a_vatom(k,3) += v[3]; a_vatom(k,4) += v[4]; a_vatom(k,5) += v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffMODKokkos<DeviceType>::v_tally3_atom(EV_FLOAT &ev, const int &i,
                                                     const int & /*j*/, const int & /*k*/,
                                                     F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji,
                                                     F_FLOAT *drjk) const
{
  F_FLOAT v[6];

  v[0] = THIRD * (drji[0]*fj[0] + drjk[0]*fk[0]);
  v[1] = THIRD * (drji[1]*fj[1] + drjk[1]*fk[1]);
  v[2] = THIRD * (drji[2]*fj[2] + drjk[2]*fk[2]);
  v[3] = THIRD * (drji[0]*fj[1] + drjk[0]*fk[1]);
  v[4] = THIRD * (drji[0]*fj[2] + drjk[0]*fk[2]);
  v[5] = THIRD * (drji[1]*fj[2] + drjk[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    d_vatom(i,0) += v[0]; d_vatom(i,1) += v[1]; d_vatom(i,2) += v[2];
    d_vatom(i,3) += v[3]; d_vatom(i,4) += v[4]; d_vatom(i,5) += v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairTersoffMODKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairTersoffMODKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairTersoffMODKokkos<LMPHostType>;
#endif
}
