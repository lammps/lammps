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

#include "pair_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#if defined(LMP_FEP)
  #include "pair_lj_cut_tip4p_long_soft.h"
  #include "pair_tip4p_long_soft.h"
#endif
#include "pair_coul_long.h"
#include "pair_tip4p_long.h"
#if defined(LMP_KSPACE)
  #include "pair_lj_long_tip4p_long.h"
  #include "pair_lj_cut_tip4p_long.h"
#endif
#if defined(LMP_MOLECULE)
  #include "pair_lj_cut_tip4p_cut.h"
  #include "pair_tip4p_cut.h"
#endif


#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::PairKokkos(LAMMPS *lmp):
  PairBase(lmp)
{
  Pair::respa_enable = 0;
  Pair::kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

}


/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::~PairKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    if constexpr(LJ) memoryKK->destroy_kokkos(k_cut_ljsq,PairBase::cut_ljsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
void PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::compute(int eflag_in, int vflag_in)
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
  if constexpr(LJ) k_cut_ljsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  // DEPRECATED c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;

  g_ewald_kk = static_cast<KK_FLOAT>(PairBase::g_ewald);

  // loop over neighbors of my atoms

  EV_FLOAT ev;
  if (Pair::ncoultablebits)
    ev = pair_compute<PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>,CoulLongTable<1> >
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>,CoulLongTable<0> >
      (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += static_cast<double>(ev.evdwl);
    eng_coul += static_cast<double>(ev.ecoul);
  }

  if (Pair::vflag_global) {
    Pair::virial[0] += static_cast<double>(ev.v[0]);
    Pair::virial[1] += static_cast<double>(ev.v[1]);
    Pair::virial[2] += static_cast<double>(ev.v[2]);
    Pair::virial[3] += static_cast<double>(ev.v[3]);
    Pair::virial[4] += static_cast<double>(ev.v[4]);
    Pair::virial[5] += static_cast<double>(ev.v[5]);
  }

  if (Pair::eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (Pair::vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (Pair::vflag_fdotr) pair_virial_fdotr_compute(this);

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
void PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::allocate()
{
  PairBase::allocate();

  const int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  if constexpr(LJ) {
    memory->destroy(PairBase::cut_ljsq);
    memoryKK->create_kokkos(k_cut_ljsq,PairBase::cut_ljsq,n+1,n+1,"pair:cut_ljsq");
    d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  }

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType>("PairKokkos::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
void PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  int ntable = 1;
  for (int i = 0; i < Pair::ncoultablebits; i++) ntable *= 2;

  tabinnersq_kk = static_cast<KK_FLOAT>(Pair::tabinnersq);

  auto transform_copy = [&](auto& d_view, const auto& h_array) {
    HAT::t_kkfloat_1d h_table("HostTable",ntable);
    typename AT::t_kkfloat_1d d_table("DeviceTable",ntable);
    for (int i = 0; i < ntable; i++)
      h_table(i) = static_cast<KK_FLOAT>(h_array[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_view = d_table;
  };
  transform_copy(d_rtable, rtable);
  transform_copy(d_drtable, drtable);
  transform_copy(d_ftable, ftable);
  transform_copy(d_dftable, dftable);
  transform_copy(d_ctable, ctable);
  transform_copy(d_dctable, dctable);
  transform_copy(d_etable, etable);
  transform_copy(d_detable, detable);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
void PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::init_style()
{
  PairBase::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

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

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
double PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::init_one(int i, int j)
{
  double cutone = PairBase::init_one(i,j);
  double cut_ljsqm = PairBase::cut_ljsq[i][j];

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).offset = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(PairBase::cut_ljsqm);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsqm);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cut_ljsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}











/*
template<class DeviceType, class PairBase, bool TIP4P, bool SOFT>
KOKKOS_INLINE_FUNCTION void tip4p_newsite(const KK_FLOAT *xO, const KK_FLOAT *xH1, const KK_FLOAT *xH2,
                                          KK_FLOAT *xM, const KK_FLOAT alpha)
{
  const KK_FLOAT delx1 = xH1[0] - xO[0];
  const KK_FLOAT dely1 = xH1[1] - xO[1];
  const KK_FLOAT delz1 = xH1[2] - xO[2];
  const KK_FLOAT delx2 = xH2[0] - xO[0];
  const KK_FLOAT dely2 = xH2[1] - xO[1];
  const KK_FLOAT delz2 = xH2[2] - xO[2];
  xM[0] = xO[0] + alpha * static_cast<KK_FLOAT>(0.5) * (delx1 + delx2);
  xM[1] = xO[1] + alpha * static_cast<KK_FLOAT>(0.5) * (dely1 + dely2);
  xM[2] = xO[2] + alpha * static_cast<KK_FLOAT>(0.5) * (delz1 + delz2);
}

template<class DeviceType>
struct PairTIP4PPreprocess {
  int nall, ago_zero, typeO, typeH;
  Kokkos::View<int **, Kokkos::LayoutRight, DeviceType> hneigh;
  Kokkos::View<KK_FLOAT **, Kokkos::LayoutRight, DeviceType> newsite;
  typename ArrayTypes<DeviceType>::t_kkfloat_1d_3_lr_randomread x;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread tag;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread sametag;
  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type map_hash;

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    if (i >= nall) return;
    if (ago_zero) hneigh(i, 0) = -1;
    hneigh(i, 2) = 0;
    if (type(i) != typeO) return;

    const tagint ti = tag(i);
    int iH1 = AtomKokkos::map_kokkos<DeviceType>(ti + 1, map_style, k_map_array, map_hash);
    int iH2 = AtomKokkos::map_kokkos<DeviceType>(ti + 2, map_style, k_map_array, map_hash);
    if (iH1 < 0 || iH2 < 0) return;
    if (type(iH1) != typeH || type(iH2) != typeH) return;

    iH1 = tip4p_closest_image<DeviceType>(i, iH1, x, sametag);
    iH2 = tip4p_closest_image<DeviceType>(i, iH2, x, sametag);

    hneigh(i, 0) = iH1;
    hneigh(i, 1) = iH2;
    hneigh(i, 2) = 1;

    KK_FLOAT xO[3] = {x(i, 0), x(i, 1), x(i, 2)};
    KK_FLOAT xH1[3] = {x(iH1, 0), x(iH1, 1), x(iH1, 2)};
    KK_FLOAT xH2[3] = {x(iH2, 0), x(iH2, 1), x(iH2, 2)};
    KK_FLOAT xM[3];
    tip4p_newsite<DeviceType>(xO, xH1, xH2, xM, alpha);
    newsite(i, 0) = xM[0];
    newsite(i, 1) = xM[1];
    newsite(i, 2) = xM[2];
  }
};

template<typename EatAccess, typename VatAccess>
KOKKOS_INLINE_FUNCTION void tip4p_ev_tally_tip4p(
    EV_FLOAT &ev, const int key, const int *vlist, const KK_FLOAT v[6], const KK_FLOAT ecoul,
    const KK_FLOAT alpha, const EatAccess &a_eatom, const VatAccess &a_vatom, const int eflag_atom,
    const int vflag_global, const int vflag_atom, const int eflag_global, const KK_FLOAT scale)
{
  const KK_ACC_FLOAT z = static_cast<KK_ACC_FLOAT>(scale);
  const KK_FLOAT a = alpha;
  const KK_ACC_FLOAT half = static_cast<KK_ACC_FLOAT>(0.5);
  const KK_ACC_FLOAT fourth = static_cast<KK_ACC_FLOAT>(0.25);

  if (eflag_global) ev.ecoul += static_cast<KK_ACC_FLOAT>(ecoul) * z;

  if (eflag_atom) {
    if (key == 0) {
      a_eatom[vlist[0]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul);
      a_eatom[vlist[1]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul);
    } else if (key == 1) {
      a_eatom[vlist[0]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul * (1.0 - a));
      a_eatom[vlist[1]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[2]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[3]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul);
    } else if (key == 2) {
      a_eatom[vlist[0]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul);
      a_eatom[vlist[1]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul * (1.0 - a));
      a_eatom[vlist[2]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[3]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
    } else {
      a_eatom[vlist[0]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul * (1.0 - a));
      a_eatom[vlist[1]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[2]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[3]] += z * half * static_cast<KK_ACC_FLOAT>(ecoul * (1.0 - a));
      a_eatom[vlist[4]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
      a_eatom[vlist[5]] += z * fourth * static_cast<KK_ACC_FLOAT>(ecoul * a);
    }
  }

  if (vflag_global) {
    for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(v[n]) * z;
  }

  if (vflag_atom) {
    if (key == 0) {
      for (int n = 0; n < 6; n++) {
        const KK_ACC_FLOAT t = z * half * static_cast<KK_ACC_FLOAT>(v[n]);
        a_vatom(vlist[0], n) += t;
        a_vatom(vlist[1], n) += t;
      }
    } else if (key == 1) {
      for (int n = 0; n < 6; n++) {
        a_vatom(vlist[0], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n] * (1.0 - a));
        a_vatom(vlist[1], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[2], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[3], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n]);
      }
    } else if (key == 2) {
      for (int n = 0; n < 6; n++) {
        a_vatom(vlist[0], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n]);
        a_vatom(vlist[1], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n] * (1.0 - a));
        a_vatom(vlist[2], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[3], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
      }
    } else {
      for (int n = 0; n < 6; n++) {
        a_vatom(vlist[0], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n] * (1.0 - a));
        a_vatom(vlist[1], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[2], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[3], n) += z * half * static_cast<KK_ACC_FLOAT>(v[n] * (1.0 - a));
        a_vatom(vlist[4], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
        a_vatom(vlist[5], n) += z * fourth * static_cast<KK_ACC_FLOAT>(v[n] * a);
      }
    }
  }
}

template<class DeviceType, unsigned NEIGHFLAG, int ZEROFLAG>
struct PairTIP4PLongComputeFunctor {
  typedef typename PairLJCutTIP4PLongKokkos<DeviceType>::AT AT;
  typedef EV_FLOAT value_type;
  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  using DUP = NeedDup_v<NEIGHFLAG, DeviceType>;

  PairLJCutTIP4PLongKokkos<DeviceType> c;
  NeighListKokkos<DeviceType> list;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  int inum;

  KKScatterView<KK_ACC_FLOAT *[3], typename DAT::t_kkacc_1d_3::array_layout, KKDeviceType, KKScatterSum,
                DUP>
      dup_f;
  KKScatterView<KK_ACC_FLOAT *, typename DAT::t_kkacc_1d::array_layout, KKDeviceType, KKScatterSum, DUP>
      dup_eatom;
  KKScatterView<KK_ACC_FLOAT *[6], typename DAT::t_kkacc_1d_6::array_layout, KKDeviceType, KKScatterSum,
                DUP>
      dup_vatom;

  PairTIP4PLongComputeFunctor(PairLJCutTIP4PLongKokkos<DeviceType> *c_ptr,
                              NeighListKokkos<DeviceType> *list_ptr) :
      c(*c_ptr), list(*list_ptr)
  {
    f = c.f;
    d_eatom = c.d_eatom;
    d_vatom = c.d_vatom;
    dup_f = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.d_vatom);
    inum = list.inum;
  }

  ~PairTIP4PLongComputeFunctor()
  {
    c.copymode = 1;
    list.copymode = 1;
  }

  KOKKOS_INLINE_FUNCTION int sbmask(const int &j) const { return j >> SBBITS & 3; }

  void contribute()
  {
    const int need_dup = std::is_same_v<DUP, Kokkos::Experimental::ScatterDuplicated>;
    if (need_dup) {
      Kokkos::Experimental::contribute(c.f, dup_f);
      if (c.eflag_atom) Kokkos::Experimental::contribute(c.d_eatom, dup_eatom);
      if (c.vflag_atom) Kokkos::Experimental::contribute(c.d_vatom, dup_vatom);
    }
  }

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION void ev_tally_lj(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &evdwl,
                                   const KK_FLOAT &fpair, const KK_FLOAT &delx, const KK_FLOAT &dely,
                                   const KK_FLOAT &delz) const
  {
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG, DeviceType>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG, DeviceType>::value>();

    if (c.eflag_either) {
      if (c.eflag_atom) {
        const KK_ACC_FLOAT epairhalf =
            static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * evdwl);
        if (NEWTON_PAIR || i < c.nlocal) a_eatom[i] += epairhalf;
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) a_eatom[j] += epairhalf;
      }
    }

    if (c.vflag_either) {
      const KK_FLOAT v0 = delx * delx * fpair;
      const KK_FLOAT v1 = dely * dely * fpair;
      const KK_FLOAT v2 = delz * delz * fpair;
      const KK_FLOAT v3 = delx * dely * fpair;
      const KK_FLOAT v4 = delx * delz * fpair;
      const KK_FLOAT v5 = dely * delz * fpair;
      const auto one_half = static_cast<KK_FLOAT>(0.5);
      const KK_ACC_FLOAT v_acc[6] = {
          static_cast<KK_ACC_FLOAT>(one_half * v0), static_cast<KK_ACC_FLOAT>(one_half * v1),
          static_cast<KK_ACC_FLOAT>(one_half * v2), static_cast<KK_ACC_FLOAT>(one_half * v3),
          static_cast<KK_ACC_FLOAT>(one_half * v4), static_cast<KK_ACC_FLOAT>(one_half * v5)};

      if (c.vflag_global) {
        if (NEIGHFLAG != FULL) {
          if (NEWTON_PAIR) {
            for (int n = 0; n < 6; n++) ev.v[n] += static_cast<KK_ACC_FLOAT>(2) * v_acc[n];
          } else {
            if (i < c.nlocal) {
              for (int n = 0; n < 6; n++) ev.v[n] += v_acc[n];
            }
            if (j < c.nlocal) {
              for (int n = 0; n < 6; n++) ev.v[n] += v_acc[n];
            }
          }
        } else {
          for (int n = 0; n < 6; n++) ev.v[n] += v_acc[n];
        }
      }

      if (c.vflag_atom) {
        if (NEWTON_PAIR || i < c.nlocal) {
          for (int n = 0; n < 6; n++) a_vatom(i, n) += v_acc[n];
        }
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) {
          for (int n = 0; n < 6; n++) a_vatom(j, n) += v_acc[n];
        }
      }
    }
  }

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION EV_FLOAT compute_item(const int &ii, const NeighListKokkos<DeviceType> &list) const
  {
    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG, DeviceType>::value>();
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG, DeviceType>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG, DeviceType>::value>();

    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const KK_FLOAT xtmp = c.x(i, 0);
    const KK_FLOAT ytmp = c.x(i, 1);
    const KK_FLOAT ztmp = c.x(i, 2);
    const int itype = c.type(i);
    const KK_FLOAT qtmp = c.q(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    KK_ACC_FLOAT fxtmp = 0;
    KK_ACC_FLOAT fytmp = 0;
    KK_ACC_FLOAT fztmp = 0;

    if (NEIGHFLAG == FULL && ZEROFLAG) {
      f(i, 0) = 0;
      f(i, 1) = 0;
      f(i, 2) = 0;
    }

    const bool stack = c.kokkos_ntypes <= MAX_TYPES_STACKPARAMS;
    const int typeO = c.tip4p_typeO;
    const int typeH = c.tip4p_typeH;
    const KK_FLOAT alpha = c.tip4p_alpha;

    int iH1 = 0, iH2 = 0;
    KK_FLOAT x1[3];
    if (itype == typeO) {
      iH1 = c.d_tip4p_hneigh(i, 0);
      iH2 = c.d_tip4p_hneigh(i, 1);
      x1[0] = c.d_tip4p_newsite(i, 0);
      x1[1] = c.d_tip4p_newsite(i, 1);
      x1[2] = c.d_tip4p_newsite(i, 2);
    } else {
      x1[0] = xtmp;
      x1[1] = ytmp;
      x1[2] = ztmp;
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
      const KK_FLOAT factor_coul = c.special_coul[sbmask(j)];
      j &= NEIGHMASK;
      const int jhalf =
          ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || (j < c.nlocal)));

      const KK_FLOAT delx_lj = xtmp - c.x(j, 0);
      const KK_FLOAT dely_lj = ytmp - c.x(j, 1);
      const KK_FLOAT delz_lj = ztmp - c.x(j, 2);
      const int jtype = c.type(j);
      KK_FLOAT rsq = delx_lj * delx_lj + dely_lj * dely_lj + delz_lj * delz_lj;

      const KK_FLOAT cutsq_ij =
          stack ? c.m_cutsq[itype][jtype] : c.d_cutsq(itype, jtype);
      const KK_FLOAT cut_ljsq_ij =
          stack ? c.m_cut_ljsq[itype][jtype] : c.d_cut_ljsq(itype, jtype);

      if (rsq < cutsq_ij) {
        KK_FLOAT fpair = KK_FLOAT();
        if (rsq < cut_ljsq_ij) {
          fpair += factor_lj * c.eval_fpair(stack, rsq, itype, jtype);
          fxtmp += static_cast<KK_ACC_FLOAT>(delx_lj * fpair);
          fytmp += static_cast<KK_ACC_FLOAT>(dely_lj * fpair);
          fztmp += static_cast<KK_ACC_FLOAT>(delz_lj * fpair);
          if (jhalf) {
            a_f(j, 0) -= static_cast<KK_ACC_FLOAT>(delx_lj * fpair);
            a_f(j, 1) -= static_cast<KK_ACC_FLOAT>(dely_lj * fpair);
            a_f(j, 2) -= static_cast<KK_ACC_FLOAT>(delz_lj * fpair);
          }
        }

        if (EVFLAG) {
          KK_FLOAT evdwl = 0.0;
          if (c.eflag_either && rsq < cut_ljsq_ij) {
            evdwl = factor_lj * c.eval_evdwl(stack, rsq, itype, jtype);
            const auto scale =
                (((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || (j < c.nlocal)))
                     ? static_cast<KK_FLOAT>(1.0)
                     : static_cast<KK_FLOAT>(0.5));
            ev.evdwl += static_cast<KK_ACC_FLOAT>(scale * evdwl);
          }
          if ((c.vflag_either || c.eflag_atom) && rsq < cut_ljsq_ij)
            ev_tally_lj<EVFLAG, NEWTON_PAIR>(ev, i, j, evdwl, fpair, delx_lj, dely_lj, delz_lj);
        }
      }

      KK_FLOAT delx = delx_lj;
      KK_FLOAT dely = dely_lj;
      KK_FLOAT delz = delz_lj;
      if (rsq < c.tip4p_cut_coulsqplus) {
        int jH1 = 0, jH2 = 0;
        if (itype == typeO || jtype == typeO) {
          KK_FLOAT x2[3];
          if (jtype == typeO) {
            jH1 = c.d_tip4p_hneigh(j, 0);
            jH2 = c.d_tip4p_hneigh(j, 1);
            x2[0] = c.d_tip4p_newsite(j, 0);
            x2[1] = c.d_tip4p_newsite(j, 1);
            x2[2] = c.d_tip4p_newsite(j, 2);
          } else {
            x2[0] = c.x(j, 0);
            x2[1] = c.x(j, 1);
            x2[2] = c.x(j, 2);
          }
          delx = x1[0] - x2[0];
          dely = x1[1] - x2[1];
          delz = x1[2] - x2[2];
          rsq = delx * delx + dely * dely + delz * delz;
        }

        const KK_FLOAT cut_coulsq_ij =
            stack ? c.m_cut_coulsq[itype][jtype] : c.d_cut_coulsq(itype, jtype);
        if (rsq < cut_coulsq_ij) {
          const KK_FLOAT cforce = c.eval_fcoul(stack, rsq, j, factor_coul, qtmp);

          int key = 0;
          int vlist[6];
          int n = 0;
          KK_FLOAT v[6] = {0, 0, 0, 0, 0, 0};
          KK_FLOAT fO[3], fH[3], fd[3];

          if (itype != typeO) {
            fxtmp += static_cast<KK_ACC_FLOAT>(delx * cforce);
            fytmp += static_cast<KK_ACC_FLOAT>(dely * cforce);
            fztmp += static_cast<KK_ACC_FLOAT>(delz * cforce);
            if (EVFLAG && c.vflag_either) {
              v[0] = c.x(i, 0) * delx * cforce;
              v[1] = c.x(i, 1) * dely * cforce;
              v[2] = c.x(i, 2) * delz * cforce;
              v[3] = c.x(i, 0) * dely * cforce;
              v[4] = c.x(i, 0) * delz * cforce;
              v[5] = c.x(i, 1) * delz * cforce;
            }
            vlist[n++] = i;
          } else {
            key++;
            fd[0] = delx * cforce;
            fd[1] = dely * cforce;
            fd[2] = delz * cforce;
            fO[0] = fd[0] * (1 - alpha);
            fO[1] = fd[1] * (1 - alpha);
            fO[2] = fd[2] * (1 - alpha);
            fH[0] = static_cast<KK_FLOAT>(0.5) * alpha * fd[0];
            fH[1] = static_cast<KK_FLOAT>(0.5) * alpha * fd[1];
            fH[2] = static_cast<KK_FLOAT>(0.5) * alpha * fd[2];
            fxtmp += static_cast<KK_ACC_FLOAT>(fO[0]);
            fytmp += static_cast<KK_ACC_FLOAT>(fO[1]);
            fztmp += static_cast<KK_ACC_FLOAT>(fO[2]);
            a_f(iH1, 0) += static_cast<KK_ACC_FLOAT>(fH[0]);
            a_f(iH1, 1) += static_cast<KK_ACC_FLOAT>(fH[1]);
            a_f(iH1, 2) += static_cast<KK_ACC_FLOAT>(fH[2]);
            a_f(iH2, 0) += static_cast<KK_ACC_FLOAT>(fH[0]);
            a_f(iH2, 1) += static_cast<KK_ACC_FLOAT>(fH[1]);
            a_f(iH2, 2) += static_cast<KK_ACC_FLOAT>(fH[2]);
            if (EVFLAG && c.vflag_either) {
              v[0] = c.x(i, 0) * fO[0] + c.x(iH1, 0) * fH[0] + c.x(iH2, 0) * fH[0];
              v[1] = c.x(i, 1) * fO[1] + c.x(iH1, 1) * fH[1] + c.x(iH2, 1) * fH[1];
              v[2] = c.x(i, 2) * fO[2] + c.x(iH1, 2) * fH[2] + c.x(iH2, 2) * fH[2];
              v[3] = c.x(i, 0) * fO[1] + c.x(iH1, 0) * fH[1] + c.x(iH2, 0) * fH[1];
              v[4] = c.x(i, 0) * fO[2] + c.x(iH1, 0) * fH[2] + c.x(iH2, 0) * fH[2];
              v[5] = c.x(i, 1) * fO[2] + c.x(iH1, 1) * fH[2] + c.x(iH2, 1) * fH[2];
            }
            vlist[n++] = i;
            vlist[n++] = iH1;
            vlist[n++] = iH2;
          }

          if (jtype != typeO) {
            if (jhalf) {
              a_f(j, 0) -= static_cast<KK_ACC_FLOAT>(delx * cforce);
              a_f(j, 1) -= static_cast<KK_ACC_FLOAT>(dely * cforce);
              a_f(j, 2) -= static_cast<KK_ACC_FLOAT>(delz * cforce);
            }
            if (EVFLAG && c.vflag_either) {
              v[0] -= c.x(j, 0) * delx * cforce;
              v[1] -= c.x(j, 1) * dely * cforce;
              v[2] -= c.x(j, 2) * delz * cforce;
              v[3] -= c.x(j, 0) * dely * cforce;
              v[4] -= c.x(j, 0) * delz * cforce;
              v[5] -= c.x(j, 1) * delz * cforce;
            }
            vlist[n++] = j;
          } else {
            key += 2;
            fd[0] = -delx * cforce;
            fd[1] = -dely * cforce;
            fd[2] = -delz * cforce;
            fO[0] = fd[0] * (1 - alpha);
            fO[1] = fd[1] * (1 - alpha);
            fO[2] = fd[2] * (1 - alpha);
            fH[0] = static_cast<KK_FLOAT>(0.5) * alpha * fd[0];
            fH[1] = static_cast<KK_FLOAT>(0.5) * alpha * fd[1];
            fH[2] = static_cast<KK_FLOAT>(0.5) * alpha * fd[2];
            if (jhalf) {
              a_f(j, 0) += static_cast<KK_ACC_FLOAT>(fO[0]);
              a_f(j, 1) += static_cast<KK_ACC_FLOAT>(fO[1]);
              a_f(j, 2) += static_cast<KK_ACC_FLOAT>(fO[2]);
              a_f(jH1, 0) += static_cast<KK_ACC_FLOAT>(fH[0]);
              a_f(jH1, 1) += static_cast<KK_ACC_FLOAT>(fH[1]);
              a_f(jH1, 2) += static_cast<KK_ACC_FLOAT>(fH[2]);
              a_f(jH2, 0) += static_cast<KK_ACC_FLOAT>(fH[0]);
              a_f(jH2, 1) += static_cast<KK_ACC_FLOAT>(fH[1]);
              a_f(jH2, 2) += static_cast<KK_ACC_FLOAT>(fH[2]);
            }
            if (EVFLAG && c.vflag_either) {
              v[0] += c.x(j, 0) * fO[0] + c.x(jH1, 0) * fH[0] + c.x(jH2, 0) * fH[0];
              v[1] += c.x(j, 1) * fO[1] + c.x(jH1, 1) * fH[1] + c.x(jH2, 1) * fH[1];
              v[2] += c.x(j, 2) * fO[2] + c.x(jH1, 2) * fH[2] + c.x(jH2, 2) * fH[2];
              v[3] += c.x(j, 0) * fO[1] + c.x(jH1, 0) * fH[1] + c.x(jH2, 0) * fH[1];
              v[4] += c.x(j, 0) * fO[2] + c.x(jH1, 0) * fH[2] + c.x(jH2, 0) * fH[2];
              v[5] += c.x(j, 1) * fO[2] + c.x(jH1, 1) * fH[2] + c.x(jH2, 1) * fH[2];
            }
            vlist[n++] = j;
            vlist[n++] = jH1;
            vlist[n++] = jH2;
          }

          KK_FLOAT ecoul = 0.0;
          if (EVFLAG && c.eflag_either)
            ecoul = c.eval_ecoul(stack, rsq, j, factor_coul, qtmp);

          if (EVFLAG && c.evflag) {
            const KK_FLOAT escale = jhalf ? static_cast<KK_FLOAT>(1.0) : static_cast<KK_FLOAT>(0.5);
            tip4p_ev_tally_tip4p(ev, key, vlist, v, ecoul, alpha, a_eatom, a_vatom, c.eflag_atom,
                                 c.vflag_global, c.vflag_atom, c.eflag_global, escale);
          }
        }
      }
    }

    a_f(i, 0) += static_cast<KK_ACC_FLOAT>(fxtmp);
    a_f(i, 1) += static_cast<KK_ACC_FLOAT>(fytmp);
    a_f(i, 2) += static_cast<KK_ACC_FLOAT>(fztmp);

    return ev;
  }

  KOKKOS_FUNCTION void operator()(const int i) const
  {
    if (c.newton_pair)
      compute_item<0, 1>(i, list);
    else
      compute_item<0, 0>(i, list);
  }

  KOKKOS_FUNCTION void operator()(const int i, EV_FLOAT &energy_virial) const
  {
    if (c.newton_pair)
      energy_virial += compute_item<1, 1>(i, list);
    else
      energy_virial += compute_item<1, 0>(i, list);
  }
};

template<class DeviceType>
EV_FLOAT pair_compute_tip4p(PairLJCutTIP4PLongKokkos<DeviceType> *fpair,
                            NeighListKokkos<DeviceType> *list)
{
  EV_FLOAT ev;
  const int inum = list->inum;

  // RangePolicy path only (no TeamPolicy): functor implements the same neighbor
  // loop pattern as PairComputeFunctor::operator()(int).
  if (fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (fpair->neighflag == FULL) {
      if (fpair->fuse_force_clear_flag) {
        PairTIP4PLongComputeFunctor<DeviceType, FULL, 1> ff(fpair, list);
        if (fpair->eflag || fpair->vflag)
          Kokkos::parallel_reduce(inum, ff, ev);
        else
          Kokkos::parallel_for(inum, ff);
        ff.contribute();
      } else {
        PairTIP4PLongComputeFunctor<DeviceType, FULL, 0> ff(fpair, list);
        if (fpair->eflag || fpair->vflag)
          Kokkos::parallel_reduce(inum, ff, ev);
        else
          Kokkos::parallel_for(inum, ff);
        ff.contribute();
      }
    } else if (fpair->neighflag == HALFTHREAD) {
      PairTIP4PLongComputeFunctor<DeviceType, HALFTHREAD, 0> ff(fpair, list);
      if (fpair->eflag || fpair->vflag)
        Kokkos::parallel_reduce(inum, ff, ev);
      else
        Kokkos::parallel_for(inum, ff);
      ff.contribute();
    } else {
      PairTIP4PLongComputeFunctor<DeviceType, HALF, 0> ff(fpair, list);
      if (fpair->eflag || fpair->vflag)
        Kokkos::parallel_reduce(inum, ff, ev);
      else
        Kokkos::parallel_for(inum, ff);
      ff.contribute();
    }
  } else {
    if (fpair->neighflag == FULL) {
      if (fpair->fuse_force_clear_flag) {
        PairTIP4PLongComputeFunctor<DeviceType, FULL, 1> ff(fpair, list);
        if (fpair->eflag || fpair->vflag)
          Kokkos::parallel_reduce(inum, ff, ev);
        else
          Kokkos::parallel_for(inum, ff);
        ff.contribute();
      } else {
        PairTIP4PLongComputeFunctor<DeviceType, FULL, 0> ff(fpair, list);
        if (fpair->eflag || fpair->vflag)
          Kokkos::parallel_reduce(inum, ff, ev);
        else
          Kokkos::parallel_for(inum, ff);
        ff.contribute();
      }
    } else if (fpair->neighflag == HALFTHREAD) {
      PairTIP4PLongComputeFunctor<DeviceType, HALFTHREAD, 0> ff(fpair, list);
      if (fpair->eflag || fpair->vflag)
        Kokkos::parallel_reduce(inum, ff, ev);
      else
        Kokkos::parallel_for(inum, ff);
      ff.contribute();
    } else {
      PairTIP4PLongComputeFunctor<DeviceType, HALF, 0> ff(fpair, list);
      if (fpair->eflag || fpair->vflag)
        Kokkos::parallel_reduce(inum, ff, ev);
      else
        Kokkos::parallel_for(inum, ff);
      ff.contribute();
    }
  }
  return ev;
}

}    // namespace LAMMPS_NS



template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag, vflag, 0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space, datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_ljsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag)
    atomKK->modified(execution_space, datamask_modify);
  else
    atomKK->modified(execution_space, F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  kokkos_ntypes = atom->ntypes;

  if (atom->nmax > nmax_tip4p) {
    nmax_tip4p = atom->nmax;
    k_tip4p_hneigh = Kokkos::DualView<int **, Kokkos::LayoutRight, DeviceType>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair:tip4p_hneigh"), nmax_tip4p, 3);
    k_tip4p_newsite = Kokkos::DualView<KK_FLOAT **, Kokkos::LayoutRight, DeviceType>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair:tip4p_newsite"), nmax_tip4p, 3);
  }
  d_tip4p_hneigh = k_tip4p_hneigh.template view<DeviceType>();
  d_tip4p_newsite = k_tip4p_newsite.template view<DeviceType>();

  atomKK->map_set();

  const int ago_zero = (neighbor->ago == 0) ? 1 : 0;
  PairTIP4PPreprocess<DeviceType> pre{};
  pre.nall = nall;
  pre.ago_zero = ago_zero;
  pre.typeO = typeO;
  pre.typeH = typeH;
  pre.hneigh = d_tip4p_hneigh;
  pre.newsite = d_tip4p_newsite;
  pre.x = x;
  pre.type = type;
  pre.tag = atomKK->k_tag.view<DeviceType>();
  pre.sametag = atomKK->k_sametag.view<DeviceType>();
  pre.map_style = atom->map_style;
  pre.k_map_array = atomKK->k_map_array;
  pre.map_hash = atomKK->k_map_hash;
  pre.alpha = static_cast<KK_FLOAT>(alpha);

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nall), pre);

  atomKK->sync(Host, TYPE_MASK | TAG_MASK);
  for (int ii = 0; ii < nall; ii++) {
    if (atom->type[ii] != typeO) continue;
    const tagint ti = atom->tag[ii];
    const int ih1 = atom->map(ti + 1);
    const int ih2 = atom->map(ti + 2);
    if (ih1 < 0 || ih2 < 0)
      error->one(FLERR, "TIP4P hydrogen is missing");
    if (atom->type[ih1] != typeH || atom->type[ih2] != typeH)
      error->one(FLERR, "TIP4P hydrogen has incorrect atom type");
  }

  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;

  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  tip4p_cut_coulsqplus =
      static_cast<KK_FLOAT>((cut_coul + 2.0 * qdist) * (cut_coul + 2.0 * qdist));
  tip4p_alpha = static_cast<KK_FLOAT>(alpha);
  tip4p_typeO = typeO;
  tip4p_typeH = typeH;

  EV_FLOAT ev = pair_compute_tip4p(this, (NeighListKokkos<DeviceType> *) list);

  if (eflag) {
    eng_vdwl += static_cast<double>(ev.evdwl);
    eng_coul += static_cast<double>(ev.ecoul);
  }

  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION KK_FLOAT PairLJCutTIP4PLongKokkos<DeviceType>::eval_fpair(bool stack,
                                                                                 const KK_FLOAT &rsq,
                                                                                 const int &itype,
                                                                                 const int &jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv * r2inv * r2inv;
  KK_FLOAT forcelj =
      r6inv * ((stack ? m_params[itype][jtype].lj1 : params(itype, jtype).lj1) * r6inv -
               (stack ? m_params[itype][jtype].lj2 : params(itype, jtype).lj2));
  return forcelj * r2inv;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION KK_FLOAT PairLJCutTIP4PLongKokkos<DeviceType>::eval_fcoul(
    bool stack, const KK_FLOAT &rsq, const int &j, const KK_FLOAT &factor_coul,
    const KK_FLOAT &qtmp) const
{
  if (ncoultablebits && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction =
        ((KK_FLOAT) rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_ftable[itable] + fraction * d_dftable[itable];
    KK_FLOAT forcecoul = qtmp * q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table2 = d_ctable[itable] + fraction * d_dctable[itable];
      const KK_FLOAT prefactor = qtmp * q[j] * table2;
      forcecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;
    }
    return forcecoul / rsq;
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT grij = g_ewald_kk * r;
    const KK_FLOAT expm2 = exp(-grij * grij);
    const KK_FLOAT t =
        static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P) * grij);
    const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
    const KK_FLOAT erfc =
        t * (static_cast<KK_FLOAT>(A1) +
             t * (static_cast<KK_FLOAT>(A2) +
                  t * (static_cast<KK_FLOAT>(A3) +
                       t * (static_cast<KK_FLOAT>(A4) + t * static_cast<KK_FLOAT>(A5))))) *
        expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp * q[j] * rinv;
    KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F) * grij * expm2);
    if (factor_coul < static_cast<KK_FLOAT>(1.0))
      forcecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;
    return forcecoul * rinv * rinv;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION KK_FLOAT PairLJCutTIP4PLongKokkos<DeviceType>::eval_evdwl(
    bool stack, const KK_FLOAT &rsq, const int &itype, const int &jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv * r2inv * r2inv;
  return r6inv * ((stack ? m_params[itype][jtype].lj3 : params(itype, jtype).lj3) * r6inv -
                  (stack ? m_params[itype][jtype].lj4 : params(itype, jtype).lj4)) -
      (stack ? m_params[itype][jtype].offset : params(itype, jtype).offset);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION KK_FLOAT PairLJCutTIP4PLongKokkos<DeviceType>::eval_ecoul(
    bool stack, const KK_FLOAT &rsq, const int &j, const KK_FLOAT &factor_coul,
    const KK_FLOAT &qtmp) const
{
  (void) stack;
  if (ncoultablebits && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction =
        ((KK_FLOAT) rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction * d_detable[itable];
    KK_FLOAT ecoul = qtmp * q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table2 = d_ctable[itable] + fraction * d_dctable[itable];
      const KK_FLOAT prefactor = qtmp * q[j] * table2;
      ecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;
    }
    return ecoul;
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT grij = g_ewald_kk * r;
    const KK_FLOAT expm2 = exp(-grij * grij);
    const KK_FLOAT t =
        static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P) * grij);
    const KK_FLOAT erfc =
        t * (static_cast<KK_FLOAT>(A1) +
             t * (static_cast<KK_FLOAT>(A2) +
                  t * (static_cast<KK_FLOAT>(A3) +
                       t * (static_cast<KK_FLOAT>(A4) + t * static_cast<KK_FLOAT>(A5))))) *
        expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp * q[j] / r;
    KK_FLOAT ecoul = prefactor * erfc;
    if (factor_coul < static_cast<KK_FLOAT>(1.0))
      ecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;
    return ecoul;
  }
}

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::allocate()
{
  PairLJCutTIP4PLong::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq, cutsq, n + 1, n + 1, "pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq, cut_ljsq, n + 1, n + 1, "pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq", n + 1, n + 1);

  k_params = Kokkos::DualView<params_lj_coul **, Kokkos::LayoutRight, DeviceType>(
      "PairLJCutTIP4PLongKokkos::params", n + 1, n + 1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul, cut_respa);

  typedef typename AT::t_kkfloat_1d table_type;
  typedef HAT::t_kkfloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);

  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(rtable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_rtable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(drtable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_drtable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(ftable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_ftable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(dftable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_dftable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(ctable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_ctable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(dctable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_dctable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(etable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_etable = d_table;
  }
  {
    host_table_type h_table("HostTable", ntable);
    table_type d_table("DeviceTable", ntable);
    for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(detable[i]);
    Kokkos::deep_copy(d_table, h_table);
    d_detable = d_table;
  }
}

template<class DeviceType>
void PairLJCutTIP4PLongKokkos<DeviceType>::init_style()
{
  PairLJCutTIP4PLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq, static_cast<KK_FLOAT>(cut_coulsq));

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa) error->all(FLERR, "Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

template<class DeviceType>
double PairLJCutTIP4PLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCutTIP4PLong::init_one(i, j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i, j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i, j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i, j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i, j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i, j).offset = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i, j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsqm);
  k_params.view_host()(i, j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);

  k_params.view_host()(j, i) = k_params.view_host()(i, j);
  if (i < MAX_TYPES_STACKPARAMS + 1 && j < MAX_TYPES_STACKPARAMS + 1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i, j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone * cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsqm);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i, j) = k_cutsq.view_host()(j, i) = cutone * cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i, j) = k_cut_ljsq.view_host()(j, i) = cut_ljsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}

*/


namespace LAMMPS_NS {

template class PairKokkos<LMPDeviceType,PairCoulLong,false,false,false>;

//template class PairKokkos<LMPDeviceType,PairTIP4PLong,false,true,false>;


#ifdef LMP_KOKKOS_GPU
//template class PairLJCutTIP4PLongKokkos<LMPHostType>;
#endif

}    // namespace LAMMPS_NS

