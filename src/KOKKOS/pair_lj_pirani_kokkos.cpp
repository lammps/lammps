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

#include "pair_lj_pirani_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJPiraniKokkos<DeviceType>::PairLJPiraniKokkos(LAMMPS *lmp) : PairLJPirani(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJPiraniKokkos<DeviceType>::~PairLJPiraniKokkos()
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
void PairLJPiraniKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev = pair_compute<PairLJPiraniKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
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

  copymode = 0;
}

/* ----------------------------------------------------------------------
   ILJ (Improved Lennard-Jones / Pirani) force and energy
   V = eps * [ gamma/(n-gamma) * (1/rx)^n - n/(n-gamma) * (1/rx)^gamma ]
   where n(x) = alpha*rx^2 + beta, rx = r/rm

   F = -dV/dr, fpair = F/r used by pair_compute template
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJPiraniKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT alpha_val   = (STACKPARAMS?m_params[itype][jtype].alpha:params(itype,jtype).alpha);
  const KK_FLOAT beta_val    = (STACKPARAMS?m_params[itype][jtype].beta:params(itype,jtype).beta);
  const KK_FLOAT gamma_val   = (STACKPARAMS?m_params[itype][jtype].gamma:params(itype,jtype).gamma);
  const KK_FLOAT rm_val      = (STACKPARAMS?m_params[itype][jtype].rm:params(itype,jtype).rm);
  const KK_FLOAT epsilon_val = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);
  const KK_FLOAT rx = r/rm_val;
  const KK_FLOAT n_x = alpha_val*rx*rx + beta_val;
  const KK_FLOAT inv_rx = static_cast<KK_FLOAT>(1.0)/rx;
  const KK_FLOAT pow_rx_n_x   = Kokkos::pow(inv_rx, n_x);
  const KK_FLOAT pow_rx_gamma = Kokkos::pow(inv_rx, gamma_val);
  const KK_FLOAT denom = n_x - gamma_val;
  const KK_FLOAT denom_sq = denom*denom;

  // filj1..6 are contributions to dV/dr; forceilj = -dV/dr
  const KK_FLOAT filj1 = -static_cast<KK_FLOAT>(2.0)*alpha_val*gamma_val*rx*pow_rx_n_x
                          /(denom_sq*rm_val);
  const KK_FLOAT filj2 = +static_cast<KK_FLOAT>(2.0)*alpha_val*rx*n_x*pow_rx_gamma
                          /(denom_sq*rm_val);
  const KK_FLOAT filj3 = -static_cast<KK_FLOAT>(2.0)*alpha_val*rx*pow_rx_gamma
                          /(rm_val*denom);
  const KK_FLOAT filj4 = +static_cast<KK_FLOAT>(2.0)*alpha_val*gamma_val*(rx/rm_val)
                          *Kokkos::log(inv_rx)*pow_rx_n_x/denom;
  const KK_FLOAT filj5 = -gamma_val*n_x*pow_rx_n_x/(r*denom);
  const KK_FLOAT filj6 = +gamma_val*n_x*pow_rx_gamma/(r*denom);

  const KK_FLOAT forceilj = -epsilon_val*(filj1+filj2+filj3+filj4+filj5+filj6);
  // fpair = F_x / delx = F/r
  return forceilj/r;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJPiraniKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT alpha_val   = (STACKPARAMS?m_params[itype][jtype].alpha:params(itype,jtype).alpha);
  const KK_FLOAT beta_val    = (STACKPARAMS?m_params[itype][jtype].beta:params(itype,jtype).beta);
  const KK_FLOAT gamma_val   = (STACKPARAMS?m_params[itype][jtype].gamma:params(itype,jtype).gamma);
  const KK_FLOAT rm_val      = (STACKPARAMS?m_params[itype][jtype].rm:params(itype,jtype).rm);
  const KK_FLOAT epsilon_val = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);
  const KK_FLOAT offset_val  = (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);
  const KK_FLOAT rx = r/rm_val;
  const KK_FLOAT n_x = alpha_val*rx*rx + beta_val;
  const KK_FLOAT inv_rx = static_cast<KK_FLOAT>(1.0)/rx;
  const KK_FLOAT denom = n_x - gamma_val;
  const KK_FLOAT ilj1 = epsilon_val*gamma_val*Kokkos::pow(inv_rx, n_x)/denom;
  const KK_FLOAT ilj2 = -epsilon_val*n_x*Kokkos::pow(inv_rx, gamma_val)/denom;
  return ilj1 + ilj2 - offset_val;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJPiraniKokkos<DeviceType>::allocate()
{
  PairLJPirani::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_pirani**,Kokkos::LayoutRight,DeviceType>("PairLJPirani::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJPiraniKokkos<DeviceType>::init_style()
{
  PairLJPirani::init_style();

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
double PairLJPiraniKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJPirani::init_one(i,j);

  k_params.view_host()(i,j).alpha   = static_cast<KK_FLOAT>(alpha[i][j]);
  k_params.view_host()(i,j).beta    = static_cast<KK_FLOAT>(beta[i][j]);
  k_params.view_host()(i,j).gamma   = static_cast<KK_FLOAT>(gamma[i][j]);
  k_params.view_host()(i,j).rm      = static_cast<KK_FLOAT>(rm[i][j]);
  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).offset  = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cutsq   = static_cast<KK_FLOAT>(cutone*cutone);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJPiraniKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJPiraniKokkos<LMPHostType>;
#endif
}
