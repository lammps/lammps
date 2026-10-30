/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(gw/kk,PairGWKokkos<LMPDeviceType>);
PairStyle(gw/kk/device,PairGWKokkos<LMPDeviceType>);
PairStyle(gw/kk/host,PairGWKokkos<LMPHostType>);
// clang-format on
#else

/* ----------------------------------------------------------------------
   PairGWKokkosBase<DeviceType,PairCPUBase> holds the device kernels of
   pair styles gw/kk and gw/zbl/kk.  It derives from the CPU class (PairGW
   or PairGWZBL), which reads the potential file, and mirrors
   PairGW::compute().  The ZBL variant replaces the repulsive term and
   multiplies the attractive term with a Fermi-like switching function,
   selected through "zblflag".
------------------------------------------------------------------------- */

// clang-format off
#ifndef LMP_PAIR_GW_KOKKOS_H
#define LMP_PAIR_GW_KOKKOS_H

#include "pair_gw.h"
#include "pair_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairGWCompute{};

struct TagPairGWComputeShortNeigh{};

namespace LAMMPS_NS {

template<class DeviceType, class PairCPUBase>
class PairGWKokkosBase : public PairCPUBase {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  using Param = typename PairCPUBase::Param;

  PairGWKokkosBase(class LAMMPS *lmp) : PairCPUBase(lmp)
  {
    this->respa_enable = 0;

    this->kokkosable = 1;
    this->atomKK = (AtomKokkos *) this->atom;
    this->execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
    this->datamask_read = X_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
    this->datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

    zblflag = 0;
    zbl_a0 = zbl_epsilon0 = zbl_e = 0.0;
  }

  ~PairGWKokkosBase() override
  {
    if (this->copymode) return;
    this->memoryKK->destroy_kokkos(k_eatom,this->eatom);
    this->memoryKK->destroy_kokkos(k_vatom,this->vatom);
    this->eatom = nullptr;
    this->vatom = nullptr;
  }

  void compute(int eflag_in, int vflag_in) override;
  void coeff(int narg, char **arg) override;
  void init_style() override;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGWCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGWCompute<NEIGHFLAG,EVFLAG>, const int &ii) const {
    EV_FLOAT ev;
    this->template operator()<NEIGHFLAG,EVFLAG>(TagPairGWCompute<NEIGHFLAG,EVFLAG>(), ii, ev);
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairGWComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
                KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drij, KK_FLOAT *drik) const;

 protected:
  void setup_params() override;

  // device versions of the PairGW (and PairGWZBL) functions

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_fc_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_fc_d_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT F_fermi_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT F_fermi_d_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_fa_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_fa_d_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_gijk_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT gw_gijk_d_k(const Param &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ex_delr_k(const Param &, const KK_FLOAT &, const KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void repulsive_k(const Param &, const KK_FLOAT &, KK_FLOAT &, const int &, KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT zeta_k(const Param &, const KK_FLOAT &, const KK_FLOAT &, KK_FLOAT *, KK_FLOAT *) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void force_zeta_k(const Param &, const KK_FLOAT &, const KK_FLOAT &, KK_FLOAT &, KK_FLOAT &,
                    const int &, KK_FLOAT &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void attractive_k(const Param &, const KK_FLOAT &, const KK_FLOAT &, const KK_FLOAT &,
                    KK_FLOAT *, KK_FLOAT *, KK_ACC_FLOAT *, KK_ACC_FLOAT *, KK_ACC_FLOAT *) const;

  int zblflag;
  double zbl_a0, zbl_epsilon0, zbl_e;    // physical constants of the ZBL variant

  typename AT::t_int_3d_randomread d_elem3param;
  typename AT::t_int_1d_randomread d_map;

  typedef Kokkos::DualView<Param*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  KK_FLOAT cutmax_sq;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  int inum;
  typename AT::t_int_2d_dl d_neighbors_short;
  typename AT::t_int_1d d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairGWKokkosBase>(PairGWKokkosBase*);
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
class PairGWKokkos : public PairGWKokkosBase<DeviceType,PairGW> {
 public:
  PairGWKokkos(class LAMMPS *lmp) : PairGWKokkosBase<DeviceType,PairGW>(lmp) {}
};

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
void PairGWKokkosBase<DeviceType,PairCPUBase>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  this->ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (this->eflag_atom) {
    this->memoryKK->destroy_kokkos(k_eatom,this->eatom);
    this->memoryKK->create_kokkos(k_eatom,this->eatom,this->maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (this->vflag_atom) {
    this->memoryKK->destroy_kokkos(k_vatom,this->vatom);
    this->memoryKK->create_kokkos(k_vatom,this->vatom,this->maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  this->atomKK->sync(this->execution_space,this->datamask_read);
  if (eflag || vflag) this->atomKK->modified(this->execution_space,this->datamask_modify);
  else this->atomKK->modified(this->execution_space,F_MASK);

  x = this->atomKK->k_x.template view<DeviceType>();
  f = this->atomKK->k_f.template view<DeviceType>();
  tag = this->atomKK->k_tag.template view<DeviceType>();
  type = this->atomKK->k_type.template view<DeviceType>();
  nlocal = this->atom->nlocal;
  nall = this->atom->nlocal + this->atom->nghost;

  inum = this->list->inum;
  const int ignum = inum + this->list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  need_dup = this->lmp->kokkos->template need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  cutmax_sq = static_cast<KK_FLOAT>(this->cutmax * this->cutmax);

  this->copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  int chunk_size = 0;
  if (this->lmp->kokkos->threads_per_atom_set)
    chunk_size = this->lmp->kokkos->threads_per_atom;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int) d_neighbors_short.extent(1) < max_neighs) ||
      ((int) d_neighbors_short.extent(0) < ignum)) {
    d_neighbors_short = typename AT::t_int_2d_dl("GW::neighbors_short",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short.extent(0) < ignum)
    d_numneigh_short = typename AT::t_int_1d("GW::numneighs_short",ignum*1.2);

  if (chunk_size)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairGWComputeShortNeigh>(0,inum,Kokkos::ChunkSize(chunk_size)), *this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairGWComputeShortNeigh>(0,inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (this->evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALF,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALF,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALF,0> >(0,inum),*this);
    }
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (this->evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALFTHREAD,1> >(0,inum),*this,ev);
    else {
      if (chunk_size)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALFTHREAD,0> >(0,inum,Kokkos::ChunkSize(chunk_size)),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGWCompute<HALFTHREAD,0> >(0,inum),*this);
    }
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (this->eflag_global) this->eng_vdwl += static_cast<double>(ev_all.evdwl);
  if (this->vflag_global) {
    for (int n = 0; n < 6; n++) this->virial[n] += static_cast<double>(ev_all.v[n]);
  }

  if (this->eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (this->vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (this->vflag_fdotr) pair_virial_fdotr_compute(this);

  this->copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f            = {};
    dup_eatom        = {};
    dup_vatom        = {};
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::operator()(TagPairGWComputeShortNeigh, const int& ii) const {
  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  const int jnum = d_numneigh[i];
  int inside = 0;
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq <= cutmax_sq) {
      d_neighbors_short(ii,inside) = j;
      inside++;
    }
  }
  d_numneigh_short(ii) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::operator()(TagPairGWCompute<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fi[3],fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;
  KK_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];
  const tagint itag = tag[i];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  const int jnum = d_numneigh_short[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  // two-body interactions, skip half of them

  for (int jj = 0; jj < jnum; jj++) {
    const int j = d_neighbors_short(ii,jj);
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

    const int iparam_ij = d_elem3param(itype,jtype,jtype);
    if (rsq > static_cast<KK_FLOAT>(d_params[iparam_ij].cutsq)) continue;

    repulsive_k(d_params[iparam_ij],rsq,fpair,eflag,evdwl);

    fxtmpi += static_cast<KK_ACC_FLOAT>(delx*fpair);
    fytmpi += static_cast<KK_ACC_FLOAT>(dely*fpair);
    fztmpi += static_cast<KK_ACC_FLOAT>(delz*fpair);
    a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
    a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
    a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (this->vflag_either || this->eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  // three-body interactions
  // skip immediately if I-J is not within cutoff

  for (int jj = 0; jj < jnum; jj++) {
    const int j = d_neighbors_short(ii,jj);
    const int jtype = d_map[type[j]];
    const int iparam_ij = d_elem3param(itype,jtype,jtype);

    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
    if (rsq1 > static_cast<KK_FLOAT>(d_params[iparam_ij].cutsq)) continue;

    // accumulate bondorder zeta for each i-j interaction via loop over k

    KK_FLOAT zeta_ij = 1.0;

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      const int k = d_neighbors_short(ii,kk);
      const int ktype = d_map[type[k]];
      const int iparam_ijk = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
      if (rsq2 > static_cast<KK_FLOAT>(d_params[iparam_ijk].cutsq)) continue;

      zeta_ij += zeta_k(d_params[iparam_ijk],rsq1,rsq2,delr1,delr2);
    }

    // pairwise force due to zeta

    KK_FLOAT prefactor;
    force_zeta_k(d_params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

    fxtmpi += static_cast<KK_ACC_FLOAT>(delr1[0]*fpair);
    fytmpi += static_cast<KK_ACC_FLOAT>(delr1[1]*fpair);
    fztmpi += static_cast<KK_ACC_FLOAT>(delr1[2]*fpair);
    KK_ACC_FLOAT fxtmpj = -static_cast<KK_ACC_FLOAT>(delr1[0]*fpair);
    KK_ACC_FLOAT fytmpj = -static_cast<KK_ACC_FLOAT>(delr1[1]*fpair);
    KK_ACC_FLOAT fztmpj = -static_cast<KK_ACC_FLOAT>(delr1[2]*fpair);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (this->vflag_either || this->eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,-fpair,-delr1[0],-delr1[1],-delr1[2]);
    }

    // attractive term via loop over k

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      const int k = d_neighbors_short(ii,kk);
      const int ktype = d_map[type[k]];
      const int iparam_ijk = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
      if (rsq2 > static_cast<KK_FLOAT>(d_params[iparam_ijk].cutsq)) continue;

      attractive_k(d_params[iparam_ijk],prefactor,rsq1,rsq2,delr1,delr2,fi,fj,fk);

      fxtmpi += fi[0];
      fytmpi += fi[1];
      fztmpi += fi[2];
      fxtmpj += fj[0];
      fytmpj += fj[1];
      fztmpj += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (EVFLAG) {
        if (this->vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj,fk,delr1,delr2);
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

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
void PairGWKokkosBase<DeviceType,PairCPUBase>::coeff(int narg, char **arg)
{
  PairCPUBase::coeff(narg,arg);

  // sync map

  int n = this->atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.view_host();

  for (int i = 1; i <= n; i++)
    h_map[i] = this->map[i];

  k_map.modify_host();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
void PairGWKokkosBase<DeviceType,PairCPUBase>::init_style()
{
  PairCPUBase::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = this->lmp->kokkos->neighflag;
  auto request = this->neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);

  if (neighflag == FULL)
    this->error->all(FLERR,"Must use half neighbor list style with pair {}/kk",
                     zblflag ? "gw/zbl" : "gw");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
void PairGWKokkosBase<DeviceType,PairCPUBase>::setup_params()
{
  PairCPUBase::setup_params();

  // sync elem3param and params

  const int nelements = this->nelements;
  const int nparams = this->nparams;

  DAT::tdual_int_3d k_elem3param = DAT::tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  HAT::t_int_3d h_elem3param = k_elem3param.view_host();

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.view_host();

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem3param(i,j,k) = this->elem3param[i][j][k];

  for (int m = 0; m < nparams; m++)
    h_params[m] = this->params[m];

  k_elem3param.modify_host();
  k_elem3param.template sync<DeviceType>();
  k_params.modify_host();
  k_params.template sync<DeviceType>();

  d_elem3param = k_elem3param.template view<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   cutoff function and its derivative
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_fc_k(const Param &param, const KK_FLOAT &r) const
{
  const KK_FLOAT gw_R = static_cast<KK_FLOAT>(param.bigr);
  const KK_FLOAT gw_D = static_cast<KK_FLOAT>(param.bigd);

  if (r < gw_R-gw_D) return static_cast<KK_FLOAT>(1.0);
  if (r > gw_R+gw_D) return static_cast<KK_FLOAT>(0.0);
  return static_cast<KK_FLOAT>(0.5)*(static_cast<KK_FLOAT>(1.0) -
                                     Kokkos::sin(static_cast<KK_FLOAT>(MathConst::MY_PI2)*(r - gw_R)/gw_D));
}

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_fc_d_k(const Param &param, const KK_FLOAT &r) const
{
  const KK_FLOAT gw_R = static_cast<KK_FLOAT>(param.bigr);
  const KK_FLOAT gw_D = static_cast<KK_FLOAT>(param.bigd);

  if (r < gw_R-gw_D) return static_cast<KK_FLOAT>(0.0);
  if (r > gw_R+gw_D) return static_cast<KK_FLOAT>(0.0);
  return -(static_cast<KK_FLOAT>(MathConst::MY_PI4)/gw_D) *
    Kokkos::cos(static_cast<KK_FLOAT>(MathConst::MY_PI2)*(r - gw_R)/gw_D);
}

/* ----------------------------------------------------------------------
   Fermi-like smoothing function of the ZBL variant and its derivative
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::F_fermi_k(const Param &param, const KK_FLOAT &r) const
{
  return static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + Kokkos::exp(-static_cast<KK_FLOAT>(param.ZBLexpscale) *
                                              (r - static_cast<KK_FLOAT>(param.ZBLcut))));
}

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::F_fermi_d_k(const Param &param, const KK_FLOAT &r) const
{
  const KK_FLOAT expscale = static_cast<KK_FLOAT>(param.ZBLexpscale);
  const KK_FLOAT zblcut = static_cast<KK_FLOAT>(param.ZBLcut);
  return expscale*Kokkos::exp(-expscale*(r-zblcut)) /
    Kokkos::pow(static_cast<KK_FLOAT>(1.0) + Kokkos::exp(-expscale*(r-zblcut)),static_cast<KK_FLOAT>(2.0));
}

/* ----------------------------------------------------------------------
   attractive pair term and its derivative
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_fa_k(const Param &param, const KK_FLOAT &r) const
{
  if (r > static_cast<KK_FLOAT>(param.bigr + param.bigd)) return static_cast<KK_FLOAT>(0.0);
  const KK_FLOAT fa = -static_cast<KK_FLOAT>(param.bigb) * Kokkos::exp(-static_cast<KK_FLOAT>(param.lam2) * r) *
    gw_fc_k(param,r);
  if (zblflag) return fa * F_fermi_k(param,r);
  return fa;
}

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_fa_d_k(const Param &param, const KK_FLOAT &r) const
{
  if (r > static_cast<KK_FLOAT>(param.bigr + param.bigd)) return static_cast<KK_FLOAT>(0.0);
  const KK_FLOAT lam2 = static_cast<KK_FLOAT>(param.lam2);
  if (zblflag)
    return static_cast<KK_FLOAT>(param.bigb) * Kokkos::exp(-lam2 * r) *
      (lam2 * gw_fc_k(param,r) * F_fermi_k(param,r) -
       gw_fc_d_k(param,r) * F_fermi_k(param,r) - gw_fc_k(param,r) *
       F_fermi_d_k(param,r));
  return static_cast<KK_FLOAT>(param.bigb) * Kokkos::exp(-lam2 * r) *
    (lam2 * gw_fc_k(param,r) - gw_fc_d_k(param,r));
}

/* ----------------------------------------------------------------------
   angular function and its derivative
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_gijk_k(const Param &param, const KK_FLOAT &costheta) const
{
  const KK_FLOAT gw_c = static_cast<KK_FLOAT>(param.c * param.c);
  const KK_FLOAT gw_d = static_cast<KK_FLOAT>(param.d * param.d);
  const KK_FLOAT hcth = static_cast<KK_FLOAT>(param.h) - costheta;

  return static_cast<KK_FLOAT>(param.gamma) *
    (static_cast<KK_FLOAT>(1.0) + gw_c / gw_d - gw_c / (gw_d + hcth * hcth));
}

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::gw_gijk_d_k(const Param &param, const KK_FLOAT &costheta) const
{
  const KK_FLOAT gw_c = static_cast<KK_FLOAT>(param.c * param.c);
  const KK_FLOAT gw_d = static_cast<KK_FLOAT>(param.d * param.d);
  const KK_FLOAT hcth = static_cast<KK_FLOAT>(param.h) - costheta;
  const KK_FLOAT numerator = static_cast<KK_FLOAT>(-2.0) * gw_c * hcth;
  const KK_FLOAT denominator = static_cast<KK_FLOAT>(1.0) / (gw_d + hcth * hcth);
  return static_cast<KK_FLOAT>(param.gamma) * numerator * denominator * denominator;
}

/* ----------------------------------------------------------------------
   exponential of the radial difference, limited as in PairGW::zeta()
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::ex_delr_k(const Param &param, const KK_FLOAT &rij,
                                                             const KK_FLOAT &rik) const
{
  KK_FLOAT arg;
  if (param.powermint == 3)
    arg = Kokkos::pow(static_cast<KK_FLOAT>(param.lam3) * (rij-rik),static_cast<KK_FLOAT>(3.0));
  else arg = static_cast<KK_FLOAT>(param.lam3) * (rij-rik);

  if (arg > static_cast<KK_FLOAT>(69.0776)) return static_cast<KK_FLOAT>(1.e30);
  else if (arg < static_cast<KK_FLOAT>(-69.0776)) return static_cast<KK_FLOAT>(0.0);
  return Kokkos::exp(arg);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::repulsive_k(const Param &param, const KK_FLOAT &rsq,
                                                           KK_FLOAT &fforce, const int &eflag, KK_FLOAT &eng) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT tmp_fc = gw_fc_k(param,r);
  const KK_FLOAT tmp_fc_d = gw_fc_d_k(param,r);
  const KK_FLOAT lam1 = static_cast<KK_FLOAT>(param.lam1);
  const KK_FLOAT biga = static_cast<KK_FLOAT>(param.biga);
  const KK_FLOAT tmp_exp = Kokkos::exp(-lam1 * r);

  if (!zblflag) {
    fforce = -biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1) / r;
    if (eflag) eng = tmp_fc * biga * tmp_exp;
    return;
  }

  // GW repulsive portion

  const KK_FLOAT fforce_gw = biga * tmp_exp * (tmp_fc_d - tmp_fc*lam1);
  const KK_FLOAT eng_gw = tmp_fc * biga * tmp_exp;

  // ZBL repulsive portion

  const KK_FLOAT Z_i = static_cast<KK_FLOAT>(param.Z_i);
  const KK_FLOAT Z_j = static_cast<KK_FLOAT>(param.Z_j);
  const KK_FLOAT esq = static_cast<KK_FLOAT>(zbl_e*zbl_e);
  const KK_FLOAT a_ij = (static_cast<KK_FLOAT>(0.8854)*static_cast<KK_FLOAT>(zbl_a0)) /
    (Kokkos::pow(Z_i,static_cast<KK_FLOAT>(0.23)) + Kokkos::pow(Z_j,static_cast<KK_FLOAT>(0.23)));
  const KK_FLOAT premult = (Z_i * Z_j * esq)/(static_cast<KK_FLOAT>(4.0*MathConst::MY_PI)*
                                              static_cast<KK_FLOAT>(zbl_epsilon0));
  const KK_FLOAT r_ov_a = r/a_ij;
  const KK_FLOAT e1 = Kokkos::exp(static_cast<KK_FLOAT>(-3.2)*r_ov_a);
  const KK_FLOAT e2 = Kokkos::exp(static_cast<KK_FLOAT>(-0.9423)*r_ov_a);
  const KK_FLOAT e3 = Kokkos::exp(static_cast<KK_FLOAT>(-0.4029)*r_ov_a);
  const KK_FLOAT e4 = Kokkos::exp(static_cast<KK_FLOAT>(-0.2016)*r_ov_a);
  const KK_FLOAT phi = static_cast<KK_FLOAT>(0.1818)*e1 + static_cast<KK_FLOAT>(0.5099)*e2 +
    static_cast<KK_FLOAT>(0.2802)*e3 + static_cast<KK_FLOAT>(0.02817)*e4;
  const KK_FLOAT dphi = (static_cast<KK_FLOAT>(1.0)/a_ij) *
    (static_cast<KK_FLOAT>(-3.2*0.1818)*e1 - static_cast<KK_FLOAT>(0.9423*0.5099)*e2 -
     static_cast<KK_FLOAT>(0.4029*0.2802)*e3 - static_cast<KK_FLOAT>(0.2016*0.02817)*e4);
  const KK_FLOAT fforce_ZBL = premult*-phi/rsq + premult*dphi/r;
  const KK_FLOAT eng_ZBL = premult*(static_cast<KK_FLOAT>(1.0)/r)*phi;

  // combine two parts with smoothing by Fermi-like function

  const KK_FLOAT fermi = F_fermi_k(param,r);
  const KK_FLOAT fermi_d = F_fermi_d_k(param,r);

  fforce = -(-fermi_d * eng_ZBL + (static_cast<KK_FLOAT>(1.0) - fermi)*fforce_ZBL +
             fermi_d*eng_gw + fermi*fforce_gw) / r;

  if (eflag) eng = (static_cast<KK_FLOAT>(1.0) - fermi)*eng_ZBL + fermi*eng_gw;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairGWKokkosBase<DeviceType,PairCPUBase>::zeta_k(const Param &param, const KK_FLOAT &rsqij,
                                                          const KK_FLOAT &rsqik, KK_FLOAT *delrij,
                                                          KK_FLOAT *delrik) const
{
  const KK_FLOAT rij = Kokkos::sqrt(rsqij);
  const KK_FLOAT rik = Kokkos::sqrt(rsqik);
  const KK_FLOAT costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
                             delrij[2]*delrik[2]) / (rij*rik);

  return gw_fc_k(param,rik) * gw_gijk_k(param,costheta) * ex_delr_k(param,rij,rik);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::force_zeta_k(const Param &param, const KK_FLOAT &rsq,
                                                            const KK_FLOAT &zeta_ij, KK_FLOAT &fforce,
                                                            KK_FLOAT &prefactor, const int &eflag,
                                                            KK_FLOAT &eng) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT fa = gw_fa_k(param,r);
  const KK_FLOAT fa_d = gw_fa_d_k(param,r);

  // bond order term and its derivative, gw_bij() and gw_bij_d()

  const KK_FLOAT powern = static_cast<KK_FLOAT>(param.powern);
  const KK_FLOAT tmp = static_cast<KK_FLOAT>(param.beta) * zeta_ij;
  const KK_FLOAT bij = Kokkos::pow(tmp,-powern);
  const KK_FLOAT bij_d = - powern * Kokkos::pow(tmp,-powern-static_cast<KK_FLOAT>(1.0))*tmp / zeta_ij;

  fforce = static_cast<KK_FLOAT>(0.5)*bij*fa_d / r;
  prefactor = static_cast<KK_FLOAT>(-0.5)*fa * bij_d;
  if (eflag) eng = static_cast<KK_FLOAT>(0.5)*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term, combines PairGW::attractive(), gw_zetaterm_d() and
   costheta_d() with the same order of operations
------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::attractive_k(const Param &param, const KK_FLOAT &prefactor,
                                                            const KK_FLOAT &rsqij, const KK_FLOAT &rsqik,
                                                            KK_FLOAT *delrij, KK_FLOAT *delrik,
                                                            KK_ACC_FLOAT *fi, KK_ACC_FLOAT *fj,
                                                            KK_ACC_FLOAT *fk) const
{
  KK_FLOAT rij_hat[3],rik_hat[3];
  KK_FLOAT dri[3],drj[3],drk[3];
  KK_FLOAT dcosdri[3],dcosdrj[3],dcosdrk[3];

  const KK_FLOAT rij = Kokkos::sqrt(rsqij);
  const KK_FLOAT rijinv = static_cast<KK_FLOAT>(1.0)/rij;
  for (int n = 0; n < 3; n++) rij_hat[n] = rijinv*delrij[n];

  const KK_FLOAT rik = Kokkos::sqrt(rsqik);
  const KK_FLOAT rikinv = static_cast<KK_FLOAT>(1.0)/rik;
  for (int n = 0; n < 3; n++) rik_hat[n] = rikinv*delrik[n];

  const KK_FLOAT lam3 = static_cast<KK_FLOAT>(param.lam3);
  const KK_FLOAT fc = gw_fc_k(param,rik);
  const KK_FLOAT dfc = gw_fc_d_k(param,rik);
  const KK_FLOAT ex_delr = ex_delr_k(param,rij,rik);
  KK_FLOAT ex_delr_d;
  if (param.powermint == 3)
    ex_delr_d = static_cast<KK_FLOAT>(3.0)*Kokkos::pow(lam3,static_cast<KK_FLOAT>(3.0)) *
      Kokkos::pow(rij-rik,static_cast<KK_FLOAT>(2.0))*ex_delr;
  else ex_delr_d = lam3 * ex_delr;

  const KK_FLOAT cos_theta = rij_hat[0]*rik_hat[0] + rij_hat[1]*rik_hat[1] + rij_hat[2]*rik_hat[2];
  const KK_FLOAT gijk = gw_gijk_k(param,cos_theta);
  const KK_FLOAT gijk_d = gw_gijk_d_k(param,cos_theta);

  // costheta_d(): derivatives of cos(theta) wrt Ri, Rj, Rk

  const KK_FLOAT rijinv2 = static_cast<KK_FLOAT>(1.0)/rij;
  const KK_FLOAT rikinv2 = static_cast<KK_FLOAT>(1.0)/rik;
  for (int n = 0; n < 3; n++) {
    dcosdrj[n] = (-cos_theta*rij_hat[n] + rik_hat[n]) * rijinv2;
    dcosdrk[n] = (-cos_theta*rik_hat[n] + rij_hat[n]) * rikinv2;
    dcosdri[n] = -(dcosdrj[n] + dcosdrk[n]);
  }

  // derivative wrt Ri

  for (int n = 0; n < 3; n++) {
    dri[n] = -dfc*gijk*ex_delr*rik_hat[n];
    dri[n] = fc*gijk_d*ex_delr*dcosdri[n] + dri[n];
    dri[n] = fc*gijk*ex_delr_d*rik_hat[n] + dri[n];
    dri[n] = -fc*gijk*ex_delr_d*rij_hat[n] + dri[n];
    fi[n] = static_cast<KK_ACC_FLOAT>(dri[n]*prefactor);
  }

  // derivative wrt Rj

  for (int n = 0; n < 3; n++) {
    drj[n] = fc*gijk_d*ex_delr*dcosdrj[n];
    drj[n] = fc*gijk*ex_delr_d*rij_hat[n] + drj[n];
    fj[n] = static_cast<KK_ACC_FLOAT>(drj[n]*prefactor);
  }

  // derivative wrt Rk

  for (int n = 0; n < 3; n++) {
    drk[n] = dfc*gijk*ex_delr*rik_hat[n];
    drk[n] = fc*gijk_d*ex_delr*dcosdrk[n] + drk[n];
    drk[n] = -fc*gijk*ex_delr_d*rik_hat[n] + drk[n];
    fk[n] = static_cast<KK_ACC_FLOAT>(drk[n]*prefactor);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (this->eflag_atom) {
    const KK_FLOAT epairhalf = static_cast<KK_FLOAT>(0.5) * epair;
    a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairhalf);
    a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairhalf);
  }

  if (this->vflag_either) {
    KK_FLOAT v[6];
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (this->vflag_global) {
      for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(v[n]);
    }

    if (this->vflag_atom) {
      for (int n = 0; n < 6; n++) {
        const KK_ACC_FLOAT vhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*v[n]);
        a_vatom(i,n) += vhalf;
        a_vatom(j,n) += vhalf;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drij*fj + drik*fk
 ------------------------------------------------------------------------- */

template<class DeviceType, class PairCPUBase>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairGWKokkosBase<DeviceType,PairCPUBase>::v_tally3(EV_FLOAT &ev, const int &i, const int &j,
        const int &k, KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drij, KK_FLOAT *drik) const
{
  // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  KK_ACC_FLOAT v[6];

  v[0] = (static_cast<KK_ACC_FLOAT>(drij[0])*fj[0] + static_cast<KK_ACC_FLOAT>(drik[0])*fk[0]);
  v[1] = (static_cast<KK_ACC_FLOAT>(drij[1])*fj[1] + static_cast<KK_ACC_FLOAT>(drik[1])*fk[1]);
  v[2] = (static_cast<KK_ACC_FLOAT>(drij[2])*fj[2] + static_cast<KK_ACC_FLOAT>(drik[2])*fk[2]);
  v[3] = (static_cast<KK_ACC_FLOAT>(drij[0])*fj[1] + static_cast<KK_ACC_FLOAT>(drik[0])*fk[1]);
  v[4] = (static_cast<KK_ACC_FLOAT>(drij[0])*fj[2] + static_cast<KK_ACC_FLOAT>(drik[0])*fk[2]);
  v[5] = (static_cast<KK_ACC_FLOAT>(drij[1])*fj[2] + static_cast<KK_ACC_FLOAT>(drik[1])*fk[2]);

  if (this->vflag_global) {
    for (int n = 0; n < 6; n++) ev.v[n] += v[n];
  }

  if (this->vflag_atom) {
    for (int n = 0; n < 6; n++) {
      v[n] *= static_cast<KK_ACC_FLOAT>(MathConst::THIRD);
      a_vatom(i,n) += v[n];
      a_vatom(j,n) += v[n];
      a_vatom(k,n) += v[n];
    }
  }
}

}

#endif
#endif
