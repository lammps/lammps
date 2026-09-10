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

#include "pair_born_coul_dsf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PIS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBornCoulDSFKokkos<DeviceType>::PairBornCoulDSFKokkos(LAMMPS *lmp)
  : PairBornCoulDSF(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBornCoulDSFKokkos<DeviceType>::~PairBornCoulDSFKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_ljsq,cut_ljsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulDSFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cut_ljsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  // the base style already computed e_shift and f_shift in init_style()

  alpha_kk = static_cast<KK_FLOAT>(alpha);
  e_shift_kk = static_cast<KK_FLOAT>(e_shift);
  f_shift_kk = static_cast<KK_FLOAT>(f_shift);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
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

  // damped-shifted-force self-energy per atom
  for (int i = 0; i < nlocal; i++) {
    double qisq = atom->q[i]*atom->q[i];
    eng_coul += -(e_shift/2.0 + alpha/MY_PIS) * qisq * force->qqrd2e;
  }

  EV_FLOAT ev;

  copymode = 1;

  ev = pair_compute<PairBornCoulDSFKokkos<DeviceType>,void>
    (this,(NeighListKokkos<DeviceType>*)list);

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
    // add the self-energy to the per-atom energy after the device sync
    for (int i = 0; i < nlocal; i++) {
      double qisq = atom->q[i]*atom->q[i];
      eatom[i] += -(e_shift/2.0 + alpha/MY_PIS) * qisq * force->qqrd2e;
    }
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   Born VdW force
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulDSFKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rhoinv = STACKPARAMS ? m_params[itype][jtype].rhoinv : params(itype,jtype).rhoinv;
  const KK_FLOAT sigma  = STACKPARAMS ? m_params[itype][jtype].sigma  : params(itype,jtype).sigma;
  const KK_FLOAT born1  = STACKPARAMS ? m_params[itype][jtype].born1  : params(itype,jtype).born1;
  const KK_FLOAT born2  = STACKPARAMS ? m_params[itype][jtype].born2  : params(itype,jtype).born2;
  const KK_FLOAT born3  = STACKPARAMS ? m_params[itype][jtype].born3  : params(itype,jtype).born3;
  const KK_FLOAT rexp = Kokkos::exp((sigma - r) * rhoinv);
  const KK_FLOAT forceborn = born1*r*rexp - born2*r6inv + born3*r2inv*r6inv;
  return forceborn*r2inv;
}

/* ----------------------------------------------------------------------
   Wolf Coulomb force
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulDSFKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT prefactor = qqrd2e*qtmp*q(j)/r;
  const KK_FLOAT erfcd = Kokkos::exp(-alpha_kk*alpha_kk*rsq);
  const KK_FLOAT erfcc = Kokkos::erfc(alpha_kk*r);

  KK_FLOAT forcecoul = prefactor * (erfcc/r +
    static_cast<KK_FLOAT>(2.0)*alpha_kk/static_cast<KK_FLOAT>(MY_PIS) * erfcd +
    r*f_shift_kk) * r;
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return forcecoul*r2inv;
}

/* ----------------------------------------------------------------------
   Born VdW energy
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulDSFKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
               const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rhoinv = STACKPARAMS ? m_params[itype][jtype].rhoinv : params(itype,jtype).rhoinv;
  const KK_FLOAT sigma  = STACKPARAMS ? m_params[itype][jtype].sigma  : params(itype,jtype).sigma;
  const KK_FLOAT a      = STACKPARAMS ? m_params[itype][jtype].a      : params(itype,jtype).a;
  const KK_FLOAT c      = STACKPARAMS ? m_params[itype][jtype].c      : params(itype,jtype).c;
  const KK_FLOAT d      = STACKPARAMS ? m_params[itype][jtype].d      : params(itype,jtype).d;
  const KK_FLOAT offset = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;
  const KK_FLOAT rexp = Kokkos::exp((sigma - r) * rhoinv);
  return a*rexp - c*r6inv + d*r6inv*r2inv - offset;
}

/* ----------------------------------------------------------------------
   Wolf Coulomb energy
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulDSFKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
               const int& /*itype*/, const int& /*jtype*/,
               const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT prefactor = qqrd2e*qtmp*q(j)/r;
  const KK_FLOAT erfcc = Kokkos::erfc(alpha_kk*r);

  KK_FLOAT ecoul = prefactor * (erfcc - r*e_shift_kk - rsq*f_shift_kk);
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulDSFKokkos<DeviceType>::allocate()
{
  PairBornCoulDSF::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_born_wolf**,Kokkos::LayoutRight,DeviceType>(
    "PairBornCoulDSF::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulDSFKokkos<DeviceType>::init_style()
{
  PairBornCoulDSF::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

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
double PairBornCoulDSFKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBornCoulDSF::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i,j).a       = static_cast<KK_FLOAT>(a[i][j]);
  k_params.view_host()(i,j).c       = static_cast<KK_FLOAT>(c[i][j]);
  k_params.view_host()(i,j).d       = static_cast<KK_FLOAT>(d[i][j]);
  k_params.view_host()(i,j).sigma   = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).rhoinv  = static_cast<KK_FLOAT>(rhoinv[i][j]);
  k_params.view_host()(i,j).born1   = static_cast<KK_FLOAT>(born1[i][j]);
  k_params.view_host()(i,j).born2   = static_cast<KK_FLOAT>(born2[i][j]);
  k_params.view_host()(i,j).born3   = static_cast<KK_FLOAT>(born3[i][j]);
  k_params.view_host()(i,j).offset  = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cut_ljsq  = static_cast<KK_FLOAT>(cut_ljsqm);
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

namespace LAMMPS_NS {
template class PairBornCoulDSFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBornCoulDSFKokkos<LMPHostType>;
#endif
}
