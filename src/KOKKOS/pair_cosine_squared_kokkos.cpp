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
   Contributing authors: Pietro Sillano (University of Turin)
------------------------------------------------------------------------- */

#include "pair_cosine_squared_kokkos.h"

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
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCosineSquaredKokkos<DeviceType>::PairCosineSquaredKokkos(LAMMPS *lmp) : PairCosineSquared(lmp)
{
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
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairCosineSquaredKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space, datamask_modify);
  else atomKK->modified(execution_space, F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  EV_FLOAT ev = pair_compute<PairCosineSquaredKokkos<DeviceType>,void>(
      this, (NeighListKokkos<DeviceType>*) list);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
F_FLOAT PairCosineSquaredKokkos<DeviceType>::
compute_fpair(const F_FLOAT &rsq, const int &, const int &,
              const int &itype, const int &jtype) const
{
  const F_FLOAT sigma_   = STACKPARAMS ? m_params[itype][jtype].sigma   : params(itype,jtype).sigma;
  const F_FLOAT w_       = STACKPARAMS ? m_params[itype][jtype].w       : params(itype,jtype).w;
  const F_FLOAT epsilon_ = STACKPARAMS ? m_params[itype][jtype].epsilon : params(itype,jtype).epsilon;
  const F_FLOAT wcaflag_ = STACKPARAMS ? m_params[itype][jtype].wcaflag : params(itype,jtype).wcaflag;

  const F_FLOAT r = sqrt(rsq);

  if (r <= sigma_) {
    if (wcaflag_ > static_cast<F_FLOAT>(0.5)) {
      const F_FLOAT lj12_f_ = STACKPARAMS ? m_params[itype][jtype].lj12_f : params(itype,jtype).lj12_f;
      const F_FLOAT lj6_f_  = STACKPARAMS ? m_params[itype][jtype].lj6_f  : params(itype,jtype).lj6_f;
      const F_FLOAT r2inv = static_cast<F_FLOAT>(1.0) / rsq;
      const F_FLOAT r6inv = r2inv*r2inv*r2inv;
      return r6inv*(lj12_f_*r6inv - lj6_f_)*r2inv;
    }
    return static_cast<F_FLOAT>(0.0);
  }

  // cosine-squared region: sigma < r < cut
  const F_FLOAT arg = static_cast<F_FLOAT>(MY_PI)*(r - sigma_) / w_;
  return -(static_cast<F_FLOAT>(MY_PI)*epsilon_ / (static_cast<F_FLOAT>(2.0)*w_)) * sin(arg) / r;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
F_FLOAT PairCosineSquaredKokkos<DeviceType>::
compute_evdwl(const F_FLOAT &rsq, const int &, const int &,
              const int &itype, const int &jtype) const
{
  const F_FLOAT sigma_   = STACKPARAMS ? m_params[itype][jtype].sigma   : params(itype,jtype).sigma;
  const F_FLOAT w_       = STACKPARAMS ? m_params[itype][jtype].w       : params(itype,jtype).w;
  const F_FLOAT epsilon_ = STACKPARAMS ? m_params[itype][jtype].epsilon : params(itype,jtype).epsilon;
  const F_FLOAT wcaflag_ = STACKPARAMS ? m_params[itype][jtype].wcaflag : params(itype,jtype).wcaflag;

  const F_FLOAT r = sqrt(rsq);

  if (r <= sigma_) {
    if (wcaflag_ > static_cast<F_FLOAT>(0.5)) {
      const F_FLOAT lj12_e_   = STACKPARAMS ? m_params[itype][jtype].lj12_e   : params(itype,jtype).lj12_e;
      const F_FLOAT lj6_e_    = STACKPARAMS ? m_params[itype][jtype].lj6_e    : params(itype,jtype).lj6_e;
      const F_FLOAT wca_shift_= STACKPARAMS ? m_params[itype][jtype].wca_shift : params(itype,jtype).wca_shift;
      const F_FLOAT r2inv = static_cast<F_FLOAT>(1.0) / rsq;
      const F_FLOAT r6inv = r2inv*r2inv*r2inv;
      return r6inv*(lj12_e_*r6inv - lj6_e_) + wca_shift_;
    }
    return -epsilon_;
  }

  // cosine-squared region: sigma < r < cut
  const F_FLOAT cosone = cos(static_cast<F_FLOAT>(MY_PI)*(r - sigma_) /
                              (static_cast<F_FLOAT>(2.0)*w_));
  return -epsilon_*cosone*cosone;
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
  memoryKK->create_kokkos(k_cutsq, cutsq, n+1, n+1, "pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_cos_sq**,Kokkos::LayoutRight,DeviceType>(
      "PairCosineSquared::params", n+1, n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCosineSquaredKokkos<DeviceType>::init_style()
{
  PairCosineSquared::init_style();

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR, "Cannot use Kokkos pair style with rRESPA inner/middle");
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
double PairCosineSquaredKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCosineSquared::init_one(i, j);

  // wca_shift: precomputed energy shift for the WCA-only case (sigma == cut)
  const F_FLOAT wca_shift =
      (wcaflag[i][j] && sigma[i][j] == cut[i][j]) ?
      static_cast<F_FLOAT>(epsilon[i][j]) : static_cast<F_FLOAT>(0.0);

  k_params.view_host()(i,j).sigma     = static_cast<F_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).w         = static_cast<F_FLOAT>(w[i][j]);
  k_params.view_host()(i,j).epsilon   = static_cast<F_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).wcaflag   = static_cast<F_FLOAT>(wcaflag[i][j]);
  k_params.view_host()(i,j).lj12_e    = static_cast<F_FLOAT>(lj12_e[i][j]);
  k_params.view_host()(i,j).lj6_e     = static_cast<F_FLOAT>(lj6_e[i][j]);
  k_params.view_host()(i,j).lj12_f    = static_cast<F_FLOAT>(lj12_f[i][j]);
  k_params.view_host()(i,j).lj6_f     = static_cast<F_FLOAT>(lj6_f[i][j]);
  k_params.view_host()(i,j).wca_shift = wca_shift;
  k_params.view_host()(i,j).cutsq     = static_cast<F_FLOAT>(cutone*cutone);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  if (i < MAX_TYPES_STACKPARAMS+1 && j < MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<F_FLOAT>(cutone*cutone);
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
