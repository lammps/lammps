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
------------------------------------------------------------------------- */

#include "pair_lj_class2_coul_cut_soft_gapsys_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::PairLJClass2CoulCutSoftGapsysKokkos(LAMMPS *lmp):PairLJClass2CoulCutSoftGapsys(lmp)
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
PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::~PairLJClass2CoulCutSoftGapsysKokkos()
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
void PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  sigmaq_kk = static_cast<KK_FLOAT>(sigmaq);
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>,void >
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
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

}

/* ----------------------------------------------------------------------
   compute soft-core LJ pair force between atoms i and j,
   already divided by r and without the special bond factor, as in the base style
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  if (rsq < (STACKPARAMS?m_params[itype][jtype].cut_lj_innersq:params(itype,jtype).cut_lj_innersq)) {
    const KK_FLOAT b1 = (STACKPARAMS?m_params[itype][jtype].b1:params(itype,jtype).b1);
    const KK_FLOAT b2 = (STACKPARAMS?m_params[itype][jtype].b2:params(itype,jtype).b2);
    return b1 + b2 / Kokkos::sqrt(rsq);
  }
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_ljsq:params(itype,jtype).cut_ljsq)) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj2 = (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT rinv = Kokkos::sqrt(r2inv);
  const KK_FLOAT r3inv = r2inv * rinv;
  const KK_FLOAT r6inv = r3inv * r3inv;
  return r6inv * (lj1 * r3inv - lj2) * r2inv;
}

/* ----------------------------------------------------------------------
   compute soft-core coulomb pair force between atoms i and j,
   already divided by r; the soft-core radius depends on the charge product
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT qiqj = qtmp * q(j);
  const KK_FLOAT cut_inner = (static_cast<KK_FLOAT>(1.0) + sigmaq_kk * Kokkos::fabs(qiqj)) *
    (STACKPARAMS?m_params[itype][jtype].cut_factor:params(itype,jtype).cut_factor);

  if (rsq < cut_inner * cut_inner)
    return factor_coul * qqrd2e * qiqj *
      (static_cast<KK_FLOAT>(-2.0) / cut_inner + static_cast<KK_FLOAT>(3.0) / Kokkos::sqrt(rsq)) /
      (cut_inner * cut_inner);
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq)) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  return factor_coul * qqrd2e * qiqj * Kokkos::sqrt(r2inv) * r2inv;
}

/* ----------------------------------------------------------------------
   compute soft-core LJ pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const KK_FLOAT offset = (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);

  if (rsq < (STACKPARAMS?m_params[itype][jtype].cut_lj_innersq:params(itype,jtype).cut_lj_innersq)) {
    const KK_FLOAT a3 = (STACKPARAMS?m_params[itype][jtype].a3:params(itype,jtype).a3);
    const KK_FLOAT b1 = (STACKPARAMS?m_params[itype][jtype].b1:params(itype,jtype).b1);
    const KK_FLOAT b2 = (STACKPARAMS?m_params[itype][jtype].b2:params(itype,jtype).b2);
    return static_cast<KK_FLOAT>(-0.5) * b1 * rsq - b2 * Kokkos::sqrt(rsq) + a3 - offset;
  }
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_ljsq:params(itype,jtype).cut_ljsq)) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj3 = (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3);
  const KK_FLOAT lj4 = (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT rinv = Kokkos::sqrt(r2inv);
  const KK_FLOAT r3inv = r2inv * rinv;
  const KK_FLOAT r6inv = r3inv * r3inv;
  return r6inv * (lj3 * r3inv - lj4) - offset;
}

/* ----------------------------------------------------------------------
   compute soft-core coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT qiqj = qtmp * q(j);
  const KK_FLOAT cut_inner = (static_cast<KK_FLOAT>(1.0) + sigmaq_kk * Kokkos::fabs(qiqj)) *
    (STACKPARAMS?m_params[itype][jtype].cut_factor:params(itype,jtype).cut_factor);

  if (rsq < cut_inner * cut_inner)
    return factor_coul * qqrd2e * qiqj *
      ((rsq / cut_inner - static_cast<KK_FLOAT>(3.0) * Kokkos::sqrt(rsq)) / cut_inner +
       static_cast<KK_FLOAT>(3.0)) / cut_inner;
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq)) return static_cast<KK_FLOAT>(0.0);

  return factor_coul * qqrd2e * qiqj / Kokkos::sqrt(rsq);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::allocate()
{
  PairLJClass2CoulCutSoftGapsys::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_cut_ljsq = DAT::tdual_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  k_cut_coulsq = DAT::tdual_kkfloat_2d("pair:cut_coulsq",n+1,n+1);
  d_cut_coulsq = k_cut_coulsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_coul_gapsys**,Kokkos::LayoutRight,DeviceType>("PairLJClass2CoulCutSoftGapsys::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::init_style()
{
  PairLJClass2CoulCutSoftGapsys::init_style();

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

template<class DeviceType>
double PairLJClass2CoulCutSoftGapsysKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJClass2CoulCutSoftGapsys::init_one(i,j);

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).offset = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsq[i][j]);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq[i][j]);
  k_params.view_host()(i,j).cut_lj_innersq =
    static_cast<KK_FLOAT>(cut_lj_inner[i][j]*cut_lj_inner[i][j]);
  k_params.view_host()(i,j).cut_factor = static_cast<KK_FLOAT>(cut_factor[i][j]);
  k_params.view_host()(i,j).a3 = static_cast<KK_FLOAT>(coeff_a3[i][j]);
  k_params.view_host()(i,j).b1 = static_cast<KK_FLOAT>(coeff_b1[i][j]);
  k_params.view_host()(i,j).b2 = static_cast<KK_FLOAT>(coeff_b2[i][j]);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  // the kernels apply the LJ and coulomb cutoffs, pair_kokkos.h only the overall cutoff

  const double cutonesq = cutone*cutone;
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutonesq);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cutonesq);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cutonesq);
  }
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutonesq;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = static_cast<KK_FLOAT>(cutonesq);
  k_cut_ljsq.modify_host();
  k_cut_coulsq.view_host()(i,j) = k_cut_coulsq.view_host()(j,i) = static_cast<KK_FLOAT>(cutonesq);
  k_cut_coulsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJClass2CoulCutSoftGapsysKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJClass2CoulCutSoftGapsysKokkos<LMPHostType>;
#endif
}

