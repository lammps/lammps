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

#include "pair_born_coul_wolf_kokkos.h"

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
PairBornCoulWolfKokkos<DeviceType>::PairBornCoulWolfKokkos(LAMMPS *lmp)
  : PairBornCoulWolf(lmp)
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
PairBornCoulWolfKokkos<DeviceType>::~PairBornCoulWolfKokkos()
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
void PairBornCoulWolfKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  // Wolf self-energy shift factors (computed on host, used as scalars in kernel)
  m_alf = static_cast<KK_FLOAT>(alf);
  e_shift = static_cast<KK_FLOAT>(erfc(alf*cut_coul)/cut_coul);
  f_shift = static_cast<KK_FLOAT>(-(e_shift + 2.0*alf/MY_PIS *
            exp(-alf*alf*cut_coul*cut_coul)) / cut_coul);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;

  // Wolf self-energy per atom
  for (int i = 0; i < nlocal; i++) {
    double qisq = atom->q[i]*atom->q[i];
    eng_coul += -(e_shift/2.0 + m_alf/MY_PIS) * qisq * qqrd2e;
  }

  EV_FLOAT ev;

  copymode = 1;

  ev = pair_compute<PairBornCoulWolfKokkos<DeviceType>,void>
    (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
    // Add Wolf self-energy to per-atom energy after device sync
    for (int i = 0; i < nlocal; i++) {
      double qisq = atom->q[i]*atom->q[i];
      eatom[i] += -(e_shift/2.0 + m_alf/MY_PIS) * qisq * qqrd2e;
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
KK_FLOAT PairBornCoulWolfKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = 1.0/rsq;
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
KK_FLOAT PairBornCoulWolfKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r2inv = 1.0/rsq;
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) / r;
  const KK_FLOAT erfcc = erfc(m_alf*r);
  const KK_FLOAT erfcd = Kokkos::exp(-m_alf*m_alf*rsq);
  const KK_FLOAT dvdrr = (erfcc*r2inv + 2.0*m_alf/MY_PIS * erfcd/r) + f_shift;
  KK_FLOAT forcecoul = dvdrr * rsq * prefactor;
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
  return forcecoul * r2inv;
}

/* ----------------------------------------------------------------------
   Born VdW energy
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulWolfKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
               const int& itype, const int& jtype) const
{
  const KK_FLOAT r2inv = 1.0/rsq;
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
KK_FLOAT PairBornCoulWolfKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
               const int& /*itype*/, const int& /*jtype*/,
               const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT prefactor = qqrd2e * qtmp * q(j) / r;
  const KK_FLOAT erfcc = erfc(m_alf*r);
  KK_FLOAT ecoul = (erfcc - e_shift*r) * prefactor;
  if (factor_coul < static_cast<KK_FLOAT>(1.0)) ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulWolfKokkos<DeviceType>::allocate()
{
  PairBornCoulWolf::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_born_wolf**,Kokkos::LayoutRight,DeviceType>(
    "PairBornCoulWolf::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulWolfKokkos<DeviceType>::init_style()
{
  PairBornCoulWolf::init_style();

  Kokkos::deep_copy(d_cut_coulsq,cut_coulsq);

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
double PairBornCoulWolfKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBornCoulWolf::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i,j).a       = a[i][j];
  k_params.view_host()(i,j).c       = c[i][j];
  k_params.view_host()(i,j).d       = d[i][j];
  k_params.view_host()(i,j).sigma   = sigma[i][j];
  k_params.view_host()(i,j).rhoinv  = rhoinv[i][j];
  k_params.view_host()(i,j).born1   = born1[i][j];
  k_params.view_host()(i,j).born2   = born2[i][j];
  k_params.view_host()(i,j).born3   = born3[i][j];
  k_params.view_host()(i,j).offset  = offset[i][j];
  k_params.view_host()(i,j).cut_ljsq  = cut_ljsqm;
  k_params.view_host()(i,j).cut_coulsq = cut_coulsq;

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_ljsqm;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsq;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cut_ljsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairBornCoulWolfKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBornCoulWolfKokkos<LMPHostType>;
#endif
}
