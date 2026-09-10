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

#include "pair_wf_cut_kokkos.h"

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
PairWFCutKokkos<DeviceType>::PairWFCutKokkos(LAMMPS *lmp) : PairWFCut(lmp)
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
PairWFCutKokkos<DeviceType>::~PairWFCutKokkos()
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
void PairWFCutKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  EV_FLOAT ev = pair_compute<PairWFCutKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairWFCutKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const int nu_i = (STACKPARAMS?m_params[itype][jtype].nu:params(itype,jtype).nu);
  const int mu_i = (STACKPARAMS?m_params[itype][jtype].mu:params(itype,jtype).mu);
  const KK_FLOAT nm_val = (STACKPARAMS?m_params[itype][jtype].nm:params(itype,jtype).nm);
  const KK_FLOAT e0nm_val = (STACKPARAMS?m_params[itype][jtype].e0nm:params(itype,jtype).e0nm);
  const KK_FLOAT rcmu_val = (STACKPARAMS?m_params[itype][jtype].rcmu:params(itype,jtype).rcmu);
  const KK_FLOAT sigma_mu_val = (STACKPARAMS?m_params[itype][jtype].sigma_mu:params(itype,jtype).sigma_mu);
  const KK_FLOAT rminv = Kokkos::pow(r2inv, static_cast<KK_FLOAT>(mu_i));
  const KK_FLOAT rm = sigma_mu_val * rminv - static_cast<KK_FLOAT>(1.0);
  const KK_FLOAT rn = rcmu_val * rminv - static_cast<KK_FLOAT>(1.0);
  const KK_FLOAT rn2nu = Kokkos::pow(rn, static_cast<KK_FLOAT>(2*nu_i));
  const KK_FLOAT rn2nu1 = Kokkos::pow(rn, static_cast<KK_FLOAT>(2*nu_i - 1));
  const KK_FLOAT forcenm = static_cast<KK_FLOAT>(2*mu_i) * sigma_mu_val * rn2nu
    + static_cast<KK_FLOAT>(4) * nm_val * rcmu_val * rm * rn2nu1;

  return e0nm_val * forcenm * Kokkos::pow(r2inv, static_cast<KK_FLOAT>(mu_i + 1));
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairWFCutKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const int nu_i = (STACKPARAMS?m_params[itype][jtype].nu:params(itype,jtype).nu);
  const int mu_i = (STACKPARAMS?m_params[itype][jtype].mu:params(itype,jtype).mu);
  const KK_FLOAT e0nm_val = (STACKPARAMS?m_params[itype][jtype].e0nm:params(itype,jtype).e0nm);
  const KK_FLOAT rcmu_val = (STACKPARAMS?m_params[itype][jtype].rcmu:params(itype,jtype).rcmu);
  const KK_FLOAT sigma_mu_val = (STACKPARAMS?m_params[itype][jtype].sigma_mu:params(itype,jtype).sigma_mu);
  const KK_FLOAT rminv = Kokkos::pow(r2inv, static_cast<KK_FLOAT>(mu_i));
  const KK_FLOAT rm = sigma_mu_val * rminv - static_cast<KK_FLOAT>(1.0);
  const KK_FLOAT rn = rcmu_val * rminv - static_cast<KK_FLOAT>(1.0);

  return e0nm_val * rm * Kokkos::pow(rn, static_cast<KK_FLOAT>(2*nu_i));
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairWFCutKokkos<DeviceType>::allocate()
{
  PairWFCut::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_wf**,Kokkos::LayoutRight,DeviceType>("PairWFCut::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairWFCutKokkos<DeviceType>::init_style()
{
  PairWFCut::init_style();

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
double PairWFCutKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairWFCut::init_one(i,j);

  k_params.view_host()(i,j).nm = static_cast<KK_FLOAT>(nm[i][j]);
  k_params.view_host()(i,j).e0nm = static_cast<KK_FLOAT>(e0nm[i][j]);
  k_params.view_host()(i,j).rcmu = static_cast<KK_FLOAT>(rcmu[i][j]);
  k_params.view_host()(i,j).sigma_mu = static_cast<KK_FLOAT>(sigma_mu[i][j]);
  k_params.view_host()(i,j).nu = nu[i][j];
  k_params.view_host()(i,j).mu = mu[i][j];
  k_params.view_host()(i,j).cutsq = static_cast<KK_FLOAT>(cutone*cutone);
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
template class PairWFCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairWFCutKokkos<LMPHostType>;
#endif
}
