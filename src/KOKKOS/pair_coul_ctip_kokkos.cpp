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

#include "pair_coul_ctip_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PIS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulCTIPKokkos<DeviceType>::PairCoulCTIPKokkos(LAMMPS *lmp) : PairCoulCTIP(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulCTIPKokkos<DeviceType>::~PairCoulCTIPKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulCTIPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cut_ljsq.template sync<DeviceType>();
  k_cut_coulsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  alpha_kk = static_cast<KK_FLOAT>(alpha);
  cut_coul_kk = static_cast<KK_FLOAT>(cut_coul);
  const double erfcd_cut = exp(-alpha * alpha * cut_coulsq);
  const double t_cut = 1.0 / (1.0 + EWALD_P * alpha * cut_coul);
  erfcc_cut_kk = static_cast<KK_FLOAT>(t_cut * (A1 + t_cut * (A2 + t_cut * (A3 + t_cut * (A4 + t_cut * A5)))) * erfcd_cut);
  k_params_self.template sync<DeviceType>();

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairCoulCTIPKokkos<DeviceType>,void >
    (this,(NeighListKokkos<DeviceType>*)list);

  // add the self energy of the local atoms. it does not contribute to the
  // forces, so it is only needed when energies are requested

  if (eflag) {
    d_ilist = ((NeighListKokkos<DeviceType>*)list)->d_ilist;
    KK_ACC_FLOAT eself = 0;
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairCoulCTIPSelf>(0,list->inum),*this,eself);
    copymode = 0;
    ev.ecoul += eself;
  }

  if (eflag_global) eng_coul += static_cast<double>(ev.ecoul);
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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulCTIPKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const params_coul &p = STACKPARAMS ? m_params[itype][jtype] : params(itype,jtype);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT reff = Kokkos::cbrt(rsq * r + p.shieldcu_inv);
  const KK_FLOAT reffsq = reff * reff;
  const KK_FLOAT reff4 = reffsq * reffsq;
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) / r;
  const KK_FLOAT erfcd = Kokkos::exp(-alpha_kk * alpha_kk * rsq);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P) * alpha_kk * r);
  const KK_FLOAT erfcc = t * (static_cast<KK_FLOAT>(A1) + t * (static_cast<KK_FLOAT>(A2) +
    t * (static_cast<KK_FLOAT>(A3) + t * (static_cast<KK_FLOAT>(A4) + t * static_cast<KK_FLOAT>(A5))))) * erfcd;

  KK_FLOAT forcecoul = prefactor *
    (erfcc / r + static_cast<KK_FLOAT>(2.0) * alpha_kk / static_cast<KK_FLOAT>(MY_PIS) * erfcd +
     rsq * r / reff4 - static_cast<KK_FLOAT>(1.0) / r - r * p.f_shift +
     r * p.s2d_shift * (r - cut_coul_kk)) * r;
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;

  return forcecoul / rsq;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulCTIPKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const params_coul &p = STACKPARAMS ? m_params[itype][jtype] : params(itype,jtype);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT reff = Kokkos::cbrt(rsq * r + p.shieldcu_inv);
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) / r;
  const KK_FLOAT erfcd = Kokkos::exp(-alpha_kk * alpha_kk * rsq);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P) * alpha_kk * r);
  const KK_FLOAT erfcc = t * (static_cast<KK_FLOAT>(A1) + t * (static_cast<KK_FLOAT>(A2) +
    t * (static_cast<KK_FLOAT>(A3) + t * (static_cast<KK_FLOAT>(A4) + t * static_cast<KK_FLOAT>(A5))))) * erfcd;
  const KK_FLOAT dr = r - cut_coul_kk;

  KK_FLOAT ecoul = prefactor *
    (erfcc + r / reff - static_cast<KK_FLOAT>(1.0) - r * erfcc_cut_kk / cut_coul_kk - r / p.reffc +
     r / cut_coul_kk + r * p.f_shift * dr - r * p.s2d_shift * dr * dr * static_cast<KK_FLOAT>(0.5));
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) ecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor;

  return ecoul;
}

/* ----------------------------------------------------------------------
   self energy of local atom i from its charge, as PairCoulCTIP::self(),
   plus the self term of the shifted and damped Coulomb interaction
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairCoulCTIPKokkos<DeviceType>::operator()(TagPairCoulCTIPSelf, const int &ii, KK_ACC_FLOAT &ecoul) const {
  const int i = d_ilist[ii];
  const KK_FLOAT qi = q(i);
  const params_self &p = params_self_d(type(i));

  KK_FLOAT selfion;
  if (qi < p.qmin) {
    selfion = qi * ((p.chi - static_cast<KK_FLOAT>(2.0) * p.qmin * p.omega) +
                    qi * (static_cast<KK_FLOAT>(0.5) * p.eta + p.omega)) + p.qmin * p.qmin * p.omega;
  } else if (qi < p.qmax) {
    selfion = qi * (p.chi + qi * (static_cast<KK_FLOAT>(0.5) * p.eta));
  } else {
    selfion = qi * ((p.chi - static_cast<KK_FLOAT>(2.0) * p.qmax * p.omega) +
                    qi * (static_cast<KK_FLOAT>(0.5) * p.eta + p.omega)) + p.qmax * p.qmax * p.omega;
  }
  const KK_FLOAT eself = selfion + p.self_factor * qi * qi;

  if (eflag_global) ecoul += static_cast<KK_ACC_FLOAT>(eself);
  if (eflag_atom) d_eatom[i] += static_cast<KK_ACC_FLOAT>(eself);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulCTIPKokkos<DeviceType>::allocate()
{
  PairCoulCTIP::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_cut_ljsq = DAT::tdual_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  k_cut_coulsq = DAT::tdual_kkfloat_2d("pair:cut_coulsq",n+1,n+1);
  d_cut_coulsq = k_cut_coulsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulCTIP::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
  k_params_self = Kokkos::DualView<params_self*,Kokkos::LayoutRight,DeviceType>("PairCoulCTIP::params_self",n+1);
  params_self_d = k_params_self.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulCTIPKokkos<DeviceType>::init_style()
{
  PairCoulCTIP::init_style();

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
double PairCoulCTIPKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulCTIP::init_one(i,j);

  // the base class keeps its parameters per element; translate them to atom types.
  // types that are not mapped to an element (NULL) get zero parameters

  const int ielem = map[i];
  const int jelem = map[j];
  params_coul pc;
  pc.cutsq = static_cast<KK_FLOAT>(cutone*cutone);
  if ((ielem >= 0) && (jelem >= 0)) {
    const int ip = elem1param[ielem];
    const int jp = elem1param[jelem];
    pc.shieldcu_inv = static_cast<KK_FLOAT>(1.0 / shieldcu[ip][jp]);
    pc.reffc = static_cast<KK_FLOAT>(reffc[ip][jp]);
    pc.f_shift = static_cast<KK_FLOAT>(f_shift[ip][jp]);
    pc.s2d_shift = static_cast<KK_FLOAT>(s2d_shift[ip][jp]);
  }
  k_params.view_host()(i,j) = pc;
  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  // the device view "params" shadows the per-element parameters of the base class

  const Param *elem_params = PairCoulCTIP::params;
  for (const int itype : {i, j}) {
    params_self ps;
    if (map[itype] >= 0) {
      const int ip = elem1param[map[itype]];
      ps.chi = static_cast<KK_FLOAT>(elem_params[ip].chi);
      ps.eta = static_cast<KK_FLOAT>(elem_params[ip].eta);
      ps.qmin = static_cast<KK_FLOAT>(elem_params[ip].qmin);
      ps.qmax = static_cast<KK_FLOAT>(elem_params[ip].qmax);
      ps.omega = static_cast<KK_FLOAT>(elem_params[ip].omega);
      ps.self_factor = static_cast<KK_FLOAT>(self_factor[ip][ip]);
    }
    k_params_self.view_host()(itype) = ps;
  }
  k_params_self.modify_host();

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  }
  k_cutsq.view_host()(i,j) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = static_cast<KK_FLOAT>(cutone*cutone);
  k_cut_ljsq.modify_host();
  k_cut_coulsq.view_host()(i,j) = k_cut_coulsq.view_host()(j,i) = static_cast<KK_FLOAT>(cutone*cutone);
  k_cut_coulsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairCoulCTIPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulCTIPKokkos<LMPHostType>;
#endif
}

