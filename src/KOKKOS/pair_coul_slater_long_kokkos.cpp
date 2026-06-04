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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_slater_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulSlaterLongKokkos<DeviceType>::PairCoulSlaterLongKokkos(LAMMPS *lmp) : PairCoulSlaterLong(lmp)
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
PairCoulSlaterLongKokkos<DeviceType>::~PairCoulSlaterLongKokkos()
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
void PairCoulSlaterLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);

  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  lamda_kk   = static_cast<KK_FLOAT>(lamda);

  copymode = 1;
  // coul/slater/long does not use force tables
  EV_FLOAT ev = pair_compute<PairCoulSlaterLongKokkos<DeviceType>,CoulLongTable<0>>
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

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute Coulomb force with Slater correction between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulSlaterLongKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT &rsq, const int & /*i*/, const int &j,
              const int &itype, const int &jtype,
              const KK_FLOAT &factor_coul, const KK_FLOAT &qtmp) const {
  const KK_FLOAT r      = sqrt(rsq);
  const KK_FLOAT r2inv  = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT grij   = g_ewald_kk * r;
  const KK_FLOAT expm2  = exp(-grij*grij);
  const KK_FLOAT t      = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc   = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT slater_term = exp(-static_cast<KK_FLOAT>(2.0)*r/lamda_kk)
    * (static_cast<KK_FLOAT>(1.0) + (static_cast<KK_FLOAT>(2.0)*r/lamda_kk*(static_cast<KK_FLOAT>(1.0)+r/lamda_kk)));
  const KK_FLOAT scale = STACKPARAMS ? m_params[itype][jtype].scale : params(itype,jtype).scale;
  const KK_FLOAT prefactor = qqrd2e * scale * qtmp * q[j] / r;
  KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2 - slater_term);
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    forcecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor * (static_cast<KK_FLOAT>(1.0) - slater_term);
  return forcecoul * r2inv;
}

/* ----------------------------------------------------------------------
   compute Coulomb energy with Slater correction between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulSlaterLongKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT &rsq, const int & /*i*/, const int &j,
              const int &itype, const int &jtype,
              const KK_FLOAT &factor_coul, const KK_FLOAT &qtmp) const {
  const KK_FLOAT r      = sqrt(rsq);
  const KK_FLOAT grij   = g_ewald_kk * r;
  const KK_FLOAT expm2  = exp(-grij*grij);
  const KK_FLOAT t      = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc   = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;
  const KK_FLOAT scale = STACKPARAMS ? m_params[itype][jtype].scale : params(itype,jtype).scale;
  const KK_FLOAT prefactor = qqrd2e * scale * qtmp * q[j] / r;
  KK_FLOAT ecoul = prefactor * (erfc - (static_cast<KK_FLOAT>(1.0) + r/lamda_kk) * exp(-static_cast<KK_FLOAT>(2.0)*r/lamda_kk));
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    ecoul -= (static_cast<KK_FLOAT>(1.0) - factor_coul) * prefactor
             * (static_cast<KK_FLOAT>(1.0) - (static_cast<KK_FLOAT>(1.0) + r/lamda_kk) * exp(-static_cast<KK_FLOAT>(2.0)*r/lamda_kk));
  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulSlaterLongKokkos<DeviceType>::allocate()
{
  PairCoulSlaterLong::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq   = typename AT::t_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_coul_slater**,Kokkos::LayoutRight,DeviceType>("PairCoulSlaterLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulSlaterLongKokkos<DeviceType>::init_style()
{
  PairCoulSlaterLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));
  Kokkos::deep_copy(d_cut_ljsq,static_cast<KK_FLOAT>(cut_coulsq));

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa) error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
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
double PairCoulSlaterLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulSlaterLong::init_one(i,j);

  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);
  k_params.view_host()(i,j).scale      = static_cast<KK_FLOAT>(scale[i][j]);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i]       = m_cutsq[i][j]       = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_coulsq[j][i]  = m_cut_coulsq[i][j]  = static_cast<KK_FLOAT>(cut_coulsq);
    m_cut_ljsq[j][i]    = m_cut_ljsq[i][j]    = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairCoulSlaterLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulSlaterLongKokkos<LMPHostType>;
#endif
}
