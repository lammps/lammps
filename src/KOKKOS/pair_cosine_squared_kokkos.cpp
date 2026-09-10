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

#include "pair_cosine_squared_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCosineSquaredKokkos<DeviceType>::PairCosineSquaredKokkos(LAMMPS *lmp) : PairCosineSquared(lmp)
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
PairCosineSquaredKokkos<DeviceType>::~PairCosineSquaredKokkos()
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
void PairCosineSquaredKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  EV_FLOAT ev = pair_compute<PairCosineSquaredKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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
KK_FLOAT PairCosineSquaredKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT sigma_val = (STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma);

  if (r <= sigma_val) {
    if (STACKPARAMS?m_params[itype][jtype].wcaflag:params(itype,jtype).wcaflag) {
      const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      return r6inv*((STACKPARAMS?m_params[itype][jtype].lj12_f:params(itype,jtype).lj12_f)*r6inv -
                    (STACKPARAMS?m_params[itype][jtype].lj6_f:params(itype,jtype).lj6_f))*r2inv;
    } else {
      return static_cast<KK_FLOAT>(0.0);
    }
  } else {
    const KK_FLOAT dr = r - sigma_val;
    const KK_FLOAT pi_over_w =
      (STACKPARAMS?m_params[itype][jtype].pi_over_w:params(itype,jtype).pi_over_w);
    const KK_FLOAT pi_over_2w =
      (STACKPARAMS?m_params[itype][jtype].pi_over_2w:params(itype,jtype).pi_over_2w);
    return -(STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon) *
             pi_over_2w * Kokkos::sin(pi_over_w * dr) / r;
  }
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCosineSquaredKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT sigma_val = (STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma);
  const KK_FLOAT epsilon_val = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);

  if (r <= sigma_val) {
    if (STACKPARAMS?m_params[itype][jtype].wcaflag:params(itype,jtype).wcaflag) {
      const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      return r6inv*((STACKPARAMS?m_params[itype][jtype].lj12_e:params(itype,jtype).lj12_e)*r6inv -
                    (STACKPARAMS?m_params[itype][jtype].lj6_e:params(itype,jtype).lj6_e)) +
             (STACKPARAMS?m_params[itype][jtype].lj_shift:params(itype,jtype).lj_shift);
    } else {
      return -epsilon_val;
    }
  } else {
    const KK_FLOAT dr = r - sigma_val;
    const KK_FLOAT pi_over_2w =
      (STACKPARAMS?m_params[itype][jtype].pi_over_2w:params(itype,jtype).pi_over_2w);
    const KK_FLOAT cosone = Kokkos::cos(pi_over_2w * dr);
    return -epsilon_val * cosone * cosone;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCosineSquaredKokkos<DeviceType>::allocate()
{
  PairCosineSquared::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_cos_sq**,Kokkos::LayoutRight,DeviceType>("PairCosineSquared::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCosineSquaredKokkos<DeviceType>::init_style()
{
  PairCosineSquared::init_style();

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
double PairCosineSquaredKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCosineSquared::init_one(i,j);

  const int wca = wcaflag[i][j];
  const double lj_shift_val = (wca && (sigma[i][j] == cut[i][j])) ? epsilon[i][j] : 0.0;
  const double w_val = w[i][j];
  const double pi_over_2w_val = (w_val > 0.0) ? MY_PI / (2.0 * w_val) : 0.0;
  const double pi_over_w_val  = (w_val > 0.0) ? MY_PI / w_val : 0.0;

  k_params.view_host()(i,j).sigma = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).w = static_cast<KK_FLOAT>(w_val);
  k_params.view_host()(i,j).wcaflag = wca;
  k_params.view_host()(i,j).lj12_e = wca ? static_cast<KK_FLOAT>(lj12_e[i][j]) : static_cast<KK_FLOAT>(0.0);
  k_params.view_host()(i,j).lj6_e  = wca ? static_cast<KK_FLOAT>(lj6_e[i][j])  : static_cast<KK_FLOAT>(0.0);
  k_params.view_host()(i,j).lj12_f = wca ? static_cast<KK_FLOAT>(lj12_f[i][j]) : static_cast<KK_FLOAT>(0.0);
  k_params.view_host()(i,j).lj6_f  = wca ? static_cast<KK_FLOAT>(lj6_f[i][j])  : static_cast<KK_FLOAT>(0.0);
  k_params.view_host()(i,j).lj_shift = static_cast<KK_FLOAT>(lj_shift_val);
  k_params.view_host()(i,j).pi_over_2w = static_cast<KK_FLOAT>(pi_over_2w_val);
  k_params.view_host()(i,j).pi_over_w  = static_cast<KK_FLOAT>(pi_over_w_val);
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
template class PairCosineSquaredKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCosineSquaredKokkos<LMPHostType>;
#endif
}
