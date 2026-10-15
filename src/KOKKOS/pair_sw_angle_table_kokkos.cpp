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
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_sw_angle_table_kokkos.h"

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
#include "pair_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

// same values as in pair_sw_angle_table.cpp
enum { LINEAR, SPLINE };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSWAngleTableKokkos<DeviceType>::PairSWAngleTableKokkos(LAMMPS *lmp) : PairSWAngleTable(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairSWAngleTableKokkos<DeviceType>::~PairSWAngleTableKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSWAngleTableKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

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

  int chunk_size = 0;
  if (lmp->kokkos->threads_per_atom_set)
    chunk_size = lmp->kokkos->threads_per_atom;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = typename AT::t_int_2d_dl("SWAngleTable::neighbors_short",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum)
    d_numneigh_short = typename AT::t_int_1d("SWAngleTable::numneighs_short",ignum*1.2);

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairSWAngleTableComputeShortNeigh>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairSWAngleTableComputeShortNeigh>(0,inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWAngleTableCompute<HALFTHREAD,0> >(0,inum),*this);
    }
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += static_cast<double>(ev_all.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev_all.v[0]);
    virial[1] += static_cast<double>(ev_all.v[1]);
    virial[2] += static_cast<double>(ev_all.v[2]);
    virial[3] += static_cast<double>(ev_all.v[3]);
    virial[4] += static_cast<double>(ev_all.v[4]);
    virial[5] += static_cast<double>(ev_all.v[5]);
  }

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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f            = {};
    dup_eatom        = {};
    dup_vatom        = {};
  }
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::operator()(TagPairSWAngleTableComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const int itype = d_map[type[i]];
    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);

    const int jnum = d_numneigh[i];
    int inside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = d_map[type[j]];

      const KK_FLOAT delx = xtmp - x(j,0);
      const KK_FLOAT dely = ytmp - x(j,1);
      const KK_FLOAT delz = ztmp - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      const int ijparam = d_elem3param(itype,jtype,jtype);
      if (rsq < static_cast<KK_FLOAT>(d_params[ijparam].cutsq)) {
        d_neighbors_short(ii,inside) = j;
        inside++;
      }
    }
    d_numneigh_short(ii) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::operator()(TagPairSWAngleTableCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;
  KK_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];
  const tagint itag = tag[i];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  // two-body interactions, skip half of them

  const int jnum = d_numneigh_short[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(ii,jj);
    const tagint jtag = tag[j];

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2) < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1) < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const int jtype = d_map[type[j]];

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const int ijparam = d_elem3param(itype,jtype,jtype);

    twobody(d_params[ijparam],rsq,fpair,eflag,evdwl);

    fxtmpi += static_cast<KK_ACC_FLOAT>(delx*fpair);
    fytmpi += static_cast<KK_ACC_FLOAT>(dely*fpair);
    fztmpi += static_cast<KK_ACC_FLOAT>(delz*fpair);
    a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
    a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
    a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  const int jnumm1 = jnum - 1;

  for (int jj = 0; jj < jnumm1; jj++) {
    int j = d_neighbors_short(ii,jj);
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem3param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    KK_ACC_FLOAT fxtmpj = 0.0;
    KK_ACC_FLOAT fytmpj = 0.0;
    KK_ACC_FLOAT fztmpj = 0.0;

    for (int kk = jj+1; kk < jnum; kk++) {
      int k = d_neighbors_short(ii,kk);
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem3param(itype,ktype,ktype);
      const int ijkparam = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      threebody_table_kk(d_params[ijparam],d_params[ikparam],ijkparam,
                         rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];
      fxtmpj += fj[0];
      fytmpj += fj[1];
      fztmpj += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  a_f(i,0) += fxtmpi;
  a_f(i,1) += fytmpi;
  a_f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::operator()(TagPairSWAngleTableCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSWAngleTableCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSWAngleTableKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairSWAngleTable::coeff(narg,arg);

  // sync map

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.view_host();

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.modify_host();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSWAngleTableKokkos<DeviceType>::init_style()
{
  PairSWAngleTable::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair sw/angle/table/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSWAngleTableKokkos<DeviceType>::setup_params()
{
  PairSWAngleTable::setup_params();

  // sync elem3param and params

  DAT::tdual_int_3d k_elem3param = DAT::tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  HAT::t_int_3d h_elem3param = k_elem3param.view_host();

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.view_host();

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

  // copy the angle tables to the device; the views are zero-initialized and have
  // one extra entry, which the spline lookup reads in the last bin with zero weight

  int maxlength = 1;
  for (int m = 0; m < nparams; m++) maxlength = MAX(maxlength,table_params[m].tablength);

  tdual_tableparam_1d k_table_params = tdual_tableparam_1d("pair:table_params",nparams);
  auto h_table_params = k_table_params.view_host();

  DAT::tdual_kkfloat_2d *tabs[7] = {&k_ang, &k_e, &k_de, &k_f, &k_df, &k_e2, &k_f2};
  const char *names[7] = {"pair:ang", "pair:e", "pair:de", "pair:f", "pair:df", "pair:e2", "pair:f2"};
  for (int n = 0; n < 7; n++)
    *tabs[n] = DAT::tdual_kkfloat_2d(names[n],nparams,maxlength+1);

  for (int m = 0; m < nparams; m++) {
    const Table *tb = table_params[m].angtable;
    const int length = table_params[m].tablength;
    h_table_params[m].tabstyle = table_params[m].tabstyle;
    h_table_params[m].tablength = length;
    h_table_params[m].invdelta = static_cast<KK_FLOAT>(tb->invdelta);
    h_table_params[m].deltasq6 = static_cast<KK_FLOAT>(tb->deltasq6);

    double *src[7] = {tb->ang, tb->e, tb->de, tb->f, tb->df, tb->e2, tb->f2};
    for (int n = 0; n < 7; n++) {
      auto h_tab = tabs[n]->view_host();
      for (int i = 0; i < length; i++) h_tab(m,i) = static_cast<KK_FLOAT>(src[n][i]);
    }
  }

  k_table_params.modify_host();
  k_table_params.template sync<DeviceType>();
  d_table_params = k_table_params.template view<DeviceType>();

  for (int n = 0; n < 7; n++) {
    tabs[n]->modify_host();
    tabs[n]->template sync<DeviceType>();
  }
  d_ang = k_ang.template view<DeviceType>();
  d_e = k_e.template view<DeviceType>();
  d_de = k_de.template view<DeviceType>();
  d_f = k_f.template view<DeviceType>();
  d_df = k_df.template view<DeviceType>();
  d_e2 = k_e2.template view<DeviceType>();
  d_f2 = k_f2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::twobody(const Param& param, const KK_FLOAT& rsq, KK_FLOAT& fforce,
                     const int& eflag, KK_FLOAT& eng) const
{
  KK_FLOAT r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

  r = Kokkos::sqrt(rsq);

  // r can round up to the cutoff in KK_FLOAT; the limit there is zero

  if (r >= static_cast<KK_FLOAT>(param.cut)) {
    fforce = 0.0;
    if (eflag) eng = 0.0;
    return;
  }

  rinvsq = static_cast<KK_FLOAT>(1.0)/rsq;
  rp = Kokkos::pow(r,static_cast<KK_FLOAT>(-param.powerp));
  rq = Kokkos::pow(r,static_cast<KK_FLOAT>(-param.powerq));
  rainv = static_cast<KK_FLOAT>(1.0) / (r - static_cast<KK_FLOAT>(param.cut));
  rainvsq = rainv*rainv*r;
  expsrainv = Kokkos::exp(static_cast<KK_FLOAT>(param.sigma) * rainv);
  fforce = (static_cast<KK_FLOAT>(param.c1)*rp - static_cast<KK_FLOAT>(param.c2)*rq +
            (static_cast<KK_FLOAT>(param.c3)*rp -static_cast<KK_FLOAT>(param.c4)*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (static_cast<KK_FLOAT>(param.c5)*rp - static_cast<KK_FLOAT>(param.c6)*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::threebody_table_kk(const Param& paramij, const Param& paramik, const int& ijkparam,
                       const KK_FLOAT& rsq1, const KK_FLOAT& rsq2,
                       KK_FLOAT *delr1, KK_FLOAT *delr2,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, const int& eflag, KK_FLOAT& eng) const
{
  KK_FLOAT r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  KK_FLOAT r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  KK_FLOAT rinv12,cs,facexp;
  KK_FLOAT ftheta,facradtable,frad1table,frad2table,var;
  KK_FLOAT acosprime,gradj1,gradj2,gradk1,gradk2,fprimetheta;

  r1 = Kokkos::sqrt(rsq1);
  r2 = Kokkos::sqrt(rsq2);

  // either separation can round up to the cutoff, see the note in twobody()

  if ((r1 >= static_cast<KK_FLOAT>(paramij.cut)) || (r2 >= static_cast<KK_FLOAT>(paramik.cut))) {
    fj[0] = fj[1] = fj[2] = 0.0;
    fk[0] = fk[1] = fk[2] = 0.0;
    if (eflag) eng = 0.0;
    return;
  }

  rinvsq1 = static_cast<KK_FLOAT>(1.0)/rsq1;
  rainv1 = static_cast<KK_FLOAT>(1.0)/(r1 - static_cast<KK_FLOAT>(paramij.cut));
  gsrainv1 = static_cast<KK_FLOAT>(paramij.sigma_gamma) * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = Kokkos::exp(gsrainv1);

  rinvsq2 = static_cast<KK_FLOAT>(1.0)/rsq2;
  rainv2 = static_cast<KK_FLOAT>(1.0)/(r2 - static_cast<KK_FLOAT>(paramik.cut));
  gsrainv2 = static_cast<KK_FLOAT>(paramik.sigma_gamma) * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = Kokkos::exp(gsrainv2);

  facexp = expgsrainv1*expgsrainv2;

  rinv12 = static_cast<KK_FLOAT>(1.0)/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  cs = Kokkos::fmax(static_cast<KK_FLOAT>(-1.0),Kokkos::fmin(cs,static_cast<KK_FLOAT>(1.0)));

  // look up energy (f(theta), ftheta) and force (df(theta)/dtheta, fprimetheta) at
  // angle theta (var) in angle table belonging to parameter set ijkparam

  var = Kokkos::acos(cs);
  uf_lookup_kk(ijkparam, var, ftheta, fprimetheta);

  if ((cs*cs - static_cast<KK_FLOAT>(1.0)) != static_cast<KK_FLOAT>(0.0))
    acosprime = static_cast<KK_FLOAT>(1.0) / (Kokkos::sqrt(static_cast<KK_FLOAT>(1.0) - cs*cs));
  else
    acosprime = static_cast<KK_FLOAT>(0.0);

  facradtable = facexp*ftheta;
  frad1table = facradtable*gsrainvsq1;
  frad2table = facradtable*gsrainvsq2;
  gradj1 = acosprime * cs * rinvsq1 * facexp * fprimetheta;
  gradj2 = acosprime * rinv12 * facexp * fprimetheta;
  gradk1 = acosprime * cs * rinvsq2 * facexp * fprimetheta;
  gradk2 = acosprime * rinv12 * facexp * fprimetheta;

  fj[0] = static_cast<KK_ACC_FLOAT>(delr1[0]*(frad1table+gradj1)-delr2[0]*gradj2);
  fj[1] = static_cast<KK_ACC_FLOAT>(delr1[1]*(frad1table+gradj1)-delr2[1]*gradj2);
  fj[2] = static_cast<KK_ACC_FLOAT>(delr1[2]*(frad1table+gradj1)-delr2[2]*gradj2);

  fk[0] = static_cast<KK_ACC_FLOAT>(delr2[0]*(frad2table+gradk1)-delr1[0]*gradk2);
  fk[1] = static_cast<KK_ACC_FLOAT>(delr2[1]*(frad2table+gradk1)-delr1[1]*gradk2);
  fk[2] = static_cast<KK_ACC_FLOAT>(delr2[2]*(frad2table+gradk1)-delr1[2]*gradk2);

  if (eflag) eng = facradtable;
}

/* ----------------------------------------------------------------------
   interpolate the angle table of parameter set m at angle x
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::uf_lookup_kk(const int &m, const KK_FLOAT &x, KK_FLOAT &u, KK_FLOAT &f) const
{
  const TableParamKokkos &pm = d_table_params[m];

  // invdelta is based on tablength-1
  int itable = static_cast<int>(x * pm.invdelta);
  if (itable < 0) itable = 0;
  if (itable >= pm.tablength) itable = pm.tablength - 1;

  if (pm.tabstyle == LINEAR) {
    const KK_FLOAT fraction = (x - d_ang(m,itable)) * pm.invdelta;
    u = d_e(m,itable) + fraction*d_de(m,itable);
    f = d_f(m,itable) + fraction*d_df(m,itable);
  } else if (pm.tabstyle == SPLINE) {
    const KK_FLOAT b = (x - d_ang(m,itable)) * pm.invdelta;
    const KK_FLOAT a = static_cast<KK_FLOAT>(1.0) - b;
    u = a * d_e(m,itable) + b * d_e(m,itable+1) +
      ((a * a * a - a) * d_e2(m,itable) + (b * b * b - b) * d_e2(m,itable+1)) * pm.deltasq6;
    f = a * d_f(m,itable) + b * d_f(m,itable+1) +
      ((a * a * a - a) * d_f2(m,itable) + (b * b * b - b) * d_f2(m,itable+1)) * pm.deltasq6;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    const KK_FLOAT epairhalf = static_cast<KK_FLOAT>(0.5) * epair;
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairhalf);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairhalf);
  }

  if (vflag_either) {
    const KK_FLOAT v0 = delx*delx*fpair;
    const KK_FLOAT v1 = dely*dely*fpair;
    const KK_FLOAT v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair;
    const KK_FLOAT v4 = delx*delz*fpair;
    const KK_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v0);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v1);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v2);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v3);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v4);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v5);
    }

    if (vflag_atom) {
      a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v0);
      a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v1);
      a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v2);
      a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v3);
      a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v4);
      a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v5);

      a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v0);
      a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v1);
      a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v2);
      a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v3);
      a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v4);
      a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v5);
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairSWAngleTableKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
          const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                     KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const
{
  KK_FLOAT epairthird,v[6];

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    epairthird = static_cast<KK_FLOAT>(THIRD) * (evdwl + ecoul);
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairthird);
    a_eatom[k] += static_cast<KK_ACC_FLOAT>(epairthird);
  }

  if (vflag_either) {
    v[0] = drji[0]*static_cast<KK_FLOAT>(fj[0]) + drki[0]*static_cast<KK_FLOAT>(fk[0]);
    v[1] = drji[1]*static_cast<KK_FLOAT>(fj[1]) + drki[1]*static_cast<KK_FLOAT>(fk[1]);
    v[2] = drji[2]*static_cast<KK_FLOAT>(fj[2]) + drki[2]*static_cast<KK_FLOAT>(fk[2]);
    v[3] = drji[0]*static_cast<KK_FLOAT>(fj[1]) + drki[0]*static_cast<KK_FLOAT>(fk[1]);
    v[4] = drji[0]*static_cast<KK_FLOAT>(fj[2]) + drki[0]*static_cast<KK_FLOAT>(fk[2]);
    v[5] = drji[1]*static_cast<KK_FLOAT>(fj[2]) + drki[1]*static_cast<KK_FLOAT>(fk[2]);

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v[0]);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v[1]);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v[2]);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v[3]);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v[4]);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v[5]);
    }

    if (vflag_atom) {
      a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[0]); a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[1]);
      a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[2]); a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[3]);
      a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[4]); a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[5]);

      a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[0]); a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[1]);
      a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[2]); a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[3]);
      a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[4]); a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[5]);

      a_vatom(k,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[0]); a_vatom(k,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[1]);
      a_vatom(k,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[2]); a_vatom(k,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[3]);
      a_vatom(k,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[4]); a_vatom(k,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*v[5]);
    }
  }
}

namespace LAMMPS_NS {
template class PairSWAngleTableKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairSWAngleTableKokkos<LMPHostType>;
#endif
}

