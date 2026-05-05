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
#include "domain_kokkos.h"
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

#include "pair_coul_long.h"
#include "pair_tip4p_long.h"
#include "pair_lj_cut_coul_long.h"
#include "pair_lj_cut_tip4p_long.h"

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

/* ----------------------------------------------------------------------
   PAIR METHODS
------------------------------------------------------------------------- */

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

  if constexpr (TIP4P) tip4p_precompute();

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

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(PairBase::cut_coulsq));

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

  // TIP4P

  if constexpr (TIP4P) {
    const double cut_coulplus = PairBase::cut_coul + 2.0 * PairBase::qdist;
    tip4p_kk.cut_coulsqplus = static_cast<KK_FLOAT>(cut_coulplus * cut_coulplus);
    tip4p_kk.typeO = PairBase::typeO;
    tip4p_kk.typeH = PairBase::typeH;
    tip4p_kk.half_alpha = static_cast<KK_FLOAT>(0.5 * PairBase::alpha);
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
double PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::init_one(int i, int j)
{
  double cutone = PairBase::init_one(i,j);

  // LJ
  if constexpr(LJ) {
    double cut_ljsqm = PairBase::cut_ljsq[i][j];
    k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(PairBase::lj1[i][j]);
    k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(PairBase::lj2[i][j]);
    k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(PairBase::lj3[i][j]);
    k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(PairBase::lj4[i][j]);
    k_params.view_host()(i,j).offset = static_cast<KK_FLOAT>(PairBase::offset[i][j]);
    k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsqm);
    if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1)
      m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsqm);
    k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cut_ljsqm;
    k_cut_ljsq.modify_host();
  }

  // COUL
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(PairBase::cut_coulsq);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(PairBase::cut_coulsq);
  }
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();

  k_params.modify_host();
  return cutone;
}





















/* ----------------------------------------------------------------------
   PROTECTED METHODS
------------------------------------------------------------------------- */

template<class DeviceType, class PairBase, bool LJ, bool TIP4P, bool SOFT>
void PairKokkos<DeviceType,PairBase,LJ,TIP4P,SOFT>::tip4p_precompute() requires (TIP4P)
{

  if (atom->nmax > tip4p_kk.nmax) {
    tip4p_kk.nmax = atom->nmax;
    tip4p_kk.k_hneigh = DAT::tdual_int_2d("pair:tip4p_hneigh", tip4p_kk.nmax, 3);
    tip4p_kk.k_newsite = DAT::tdual_kkfloat_2d("pair:tip4p_newsite", tip4p_kk.nmax, 3);
  }
  tip4p_kk.d_hneigh = tip4p_kk.k_hneigh.template view<DeviceType>();
  tip4p_kk.d_newsite = tip4p_kk.k_newsite.template view<DeviceType>();

  atomKK->map_set();
  auto l_map_style = atom->map_style;
  auto l_map_array = atomKK->k_map_array;
  auto l_map_hash = atomKK->k_map_hash;

  auto l_tag = atomKK->k_tag.template view<DeviceType>();
  auto l_x = atomKK->k_x.template view<DeviceType>();
  auto l_type = atomKK->k_type.template view<DeviceType>();
  auto l_sametag = atomKK->k_sametag.template view<DeviceType>();

  const int l_ago_zero = (neighbor->ago == 0) ? 1 : 0;
  int flag_missing, flag_incorrect;

  auto l_tip4p_kk = tip4p_kk;

  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<DeviceType>(0, atom->nlocal + atom->nghost),
    KOKKOS_LAMBDA (const int i, int& l_flag_missing, int& l_flag_incorrect) {

      if (l_ago_zero) l_tip4p_kk.d_hneigh(i,0) = -1;
      l_tip4p_kk.d_hneigh(i,2) = 0;
      if (l_type(i) != l_tip4p_kk.typeO) return;

      const tagint ti = l_tag(i);
      int iH1 = AtomKokkos::map_kokkos<DeviceType>(ti + 1, l_map_style, l_map_array, l_map_hash);
      int iH2 = AtomKokkos::map_kokkos<DeviceType>(ti + 2, l_map_style, l_map_array, l_map_hash);

      l_flag_missing = l_flag_incorrect = 0;
      if (iH1 < 0 || iH2 < 0) {
        l_flag_missing = 1;
        return;
      }
      if (l_type(iH1) != l_tip4p_kk.typeH || l_type(iH2) != l_tip4p_kk.typeH) {
        l_flag_incorrect = 1;
        return;
      }

      iH1 = DomainKokkos::closest_image(l_x, l_sametag, i, iH1);
      iH2 = DomainKokkos::closest_image(l_x, l_sametag, i, iH2);

      l_tip4p_kk.d_hneigh(i,0) = iH1;
      l_tip4p_kk.d_hneigh(i,1) = iH2;
      l_tip4p_kk.d_hneigh(i,2) = 1;

      const KK_FLOAT xO = l_x(i,0);
      const KK_FLOAT yO = l_x(i,1);
      const KK_FLOAT zO = l_x(i,2);

      const KK_FLOAT delx1 = l_x(iH1,0) - xO;
      const KK_FLOAT dely1 = l_x(iH1,1) - yO;
      const KK_FLOAT delz1 = l_x(iH1,2) - zO;
      const KK_FLOAT delx2 = l_x(iH2,0) - xO;
      const KK_FLOAT dely2 = l_x(iH2,1) - yO;
      const KK_FLOAT delz2 = l_x(iH2,2) - zO;

      l_tip4p_kk.d_newsite(i,0) = xO + l_tip4p_kk.half_alpha * (delx1 + delx2);
      l_tip4p_kk.d_newsite(i,1) = yO + l_tip4p_kk.half_alpha * (dely1 + dely2);
      l_tip4p_kk.d_newsite(i,2) = zO + l_tip4p_kk.half_alpha * (delz1 + delz2);

    }, Kokkos::Max<int>(flag_missing), Kokkos::Max<int>(flag_incorrect)
  );

  if (flag_missing) error->one(FLERR, "TIP4P hydrogen is missing");
  if (flag_incorrect) error->one(FLERR, "TIP4P hydrogen has incorrect atom type");

}



namespace LAMMPS_NS {

  // coul/long/kk
  template class PairKokkos<LMPDeviceType,PairCoulLong,false,false,false>;

  // coul/tip4p/kk
  template class PairKokkos<LMPDeviceType,PairTIP4PLong,false,true,false>;

  // lj/cut/coul/long/kk
  template class PairKokkos<LMPDeviceType,PairLJCutCoulLong,true,false,false>;

  // lj/cut/tip4p/long/kk
  template class PairKokkos<LMPDeviceType,PairLJCutTIP4PLong,true,true,false>;


  #ifdef LMP_KOKKOS_GPU
    //template class PairLJCutTIP4PLongKokkos<LMPHostType>;
  #endif

}    // namespace LAMMPS_NS

