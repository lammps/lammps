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
   Contributing authors: Tim Linke & Niels Gronbech-Jensen (UC Davis)
                         Kokkos port: LAMMPS developers
------------------------------------------------------------------------- */

#include "fix_gjf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixGJFKokkos<DeviceType>::FixGJFKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixGJF(lmp, narg, arg),
#ifndef LMP_KOKKOS_DEBUG_RNG
  rand_pool(seed + comm->me)
#else
  rand_pool(seed + comm->me, lmp)
#endif
{
  kokkosable = 1;
  exchange_comm_device = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(random, seed + comm->me);
#endif

  // Replace base class lv (double **) with a Kokkos DualView.
  // The base constructor already allocated and zero-initialized lv.
  // Save the pointer, set base to nullptr, create DualView, copy data.
  auto lv_tmp = lv;
  lv = nullptr;

  grow_arrays(atom->nmax);   // creates k_lv; sets lv to h_view data ptr

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    k_lv.view_host()(i,0) = lv_tmp[i][0];
    k_lv.view_host()(i,1) = lv_tmp[i][1];
    k_lv.view_host()(i,2) = lv_tmp[i][2];
  }
  k_lv.modify_host();

  memory->destroy(lv_tmp);

  d_count = typename AT::t_int_scalar("gjf:count");
  h_count = Kokkos::create_mirror_view(d_count);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixGJFKokkos<DeviceType>::~FixGJFKokkos()
{
  if (copymode) return;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif

  memoryKK->destroy_kokkos(k_lv, lv);
  lv = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::init()
{
  FixGJF::init();

  if (tbiasflag)
    error->all(FLERR, "Fix gjf/kk does not yet support temperature bias removal");

  if (tstyle == 2)   // 2 == ATOM style (ATOM enum is private to fix_gjf.cpp)
    error->all(FLERR, "Fix gjf/kk does not yet support per-atom temperature");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::initial_integrate(int /* vflag */)
{
  const int nlocal = atom->nlocal;
  l_rmass_flag = (atom->rmass != nullptr) ? 1 : 0;

  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | MASK_MASK |
               RMASS_MASK | TYPE_MASK);

  x    = atomKK->k_x.view<DeviceType>();
  v    = atomKK->k_v.view<DeviceType>();
  f    = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  if (l_rmass_flag) rmass = atomKK->k_rmass.view<DeviceType>();
  else              mass  = atomKK->k_mass.view<DeviceType>();

  // vhalf mode: copy lv → v at start of step (all but the first step)
  if (!osflag && lv_allocated) {
    k_lv.sync<DeviceType>();
    d_lv = k_lv.view<DeviceType>();
    auto l_v = v;
    auto l_d_lv = d_lv;
    auto l_mask = mask;
    auto l_groupbit = groupbit;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nlocal),
      LAMMPS_LAMBDA(const int i) {
        if (l_mask[i] & l_groupbit) {
          l_v(i,0) = l_d_lv(i,0);
          l_v(i,1) = l_d_lv(i,1);
          l_v(i,2) = l_d_lv(i,2);
        }
      });
    atomKK->modified(execution_space, V_MASK);
  }

  // compute_target() handles CONSTANT/EQUAL temperature on host
  compute_target();

  // cache integration constants into member variables for use in device functor
  const double dt = update->dt;
  l_gjfc2    = static_cast<KK_FLOAT>(gjfc2);
  l_c1sqrt   = static_cast<KK_FLOAT>(sqrt(gjfc1));
  l_c3sqrt   = static_cast<KK_FLOAT>(sqrt(gjfc3));
  l_csq      = static_cast<KK_FLOAT>(sqrt(gjfc3 / gjfc1));
  l_dtf      = static_cast<KK_FLOAT>(0.5 * dt * force->ftm2v);
  l_dt       = static_cast<KK_FLOAT>(dt);
  l_tsqrt    = static_cast<KK_FLOAT>(tsqrt);
  l_t_period = static_cast<KK_FLOAT>(t_period);
  l_boltz    = static_cast<KK_FLOAT>(force->boltz);
  l_mvv2e    = static_cast<KK_FLOAT>(force->mvv2e);
  l_ftm2v    = static_cast<KK_FLOAT>(force->ftm2v);

  k_lv.sync<DeviceType>();
  d_lv = k_lv.view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixGJFInitialIntegrate>(0, nlocal), *this);
  copymode = 0;

  k_lv.modify<DeviceType>();
  atomKK->modified(execution_space, X_MASK | V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixGJFKokkos<DeviceType>::operator()(TagFixGJFInitialIntegrate, const int &i) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT m;
    if (l_rmass_flag) m = rmass[i];
    else              m = mass[type[i]];

    // random force magnitude: beta = tsqrt * sqrt(2*dt*m*boltz/t_period/mvv2e) / ftm2v
    const KK_FLOAT beta = l_tsqrt
        * Kokkos::sqrt(2.0 * l_dt * m * l_boltz / l_t_period / l_mvv2e) / l_ftm2v;

    rand_type rand_gen = rand_pool.get_state();
    const KK_FLOAT fran0 = beta * rand_gen.normal();
    const KK_FLOAT fran1 = beta * rand_gen.normal();
    const KK_FLOAT fran2 = beta * rand_gen.normal();
    rand_pool.free_state(rand_gen);

    const KK_FLOAT dtfm = l_dtf / m;

    // Eq. 24a: v += csq * dtfm * f
    v(i,0) += l_csq * dtfm * f(i,0);
    v(i,1) += l_csq * dtfm * f(i,1);
    v(i,2) += l_csq * dtfm * f(i,2);

    // Eq. 24b: x += 0.5 * csq * dt * v
    x(i,0) += 0.5 * l_csq * l_dt * v(i,0);
    x(i,1) += 0.5 * l_csq * l_dt * v(i,1);
    x(i,2) += 0.5 * l_csq * l_dt * v(i,2);

    // Eq. 24c: lv = c1sqrt * v + ftm2v * (c3sqrt / (2*m)) * fran
    const KK_FLOAT lv_coeff = l_ftm2v * l_c3sqrt / (2.0 * m);
    d_lv(i,0) = l_c1sqrt * v(i,0) + lv_coeff * fran0;
    d_lv(i,1) = l_c1sqrt * v(i,1) + lv_coeff * fran1;
    d_lv(i,2) = l_c1sqrt * v(i,2) + lv_coeff * fran2;

    // Eq. 24d: v = (gjfc2/c1sqrt) * lv + ftm2v * csq * (0.5/m) * fran
    const KK_FLOAT v_coeff1 = l_gjfc2 / l_c1sqrt;
    const KK_FLOAT v_coeff2 = l_ftm2v * l_csq * 0.5 / m;
    v(i,0) = v_coeff1 * d_lv(i,0) + v_coeff2 * fran0;
    v(i,1) = v_coeff1 * d_lv(i,1) + v_coeff2 * fran1;
    v(i,2) = v_coeff1 * d_lv(i,2) + v_coeff2 * fran2;

    // Eq. 24e: x += 0.5 * csq * dt * v
    x(i,0) += 0.5 * l_csq * l_dt * v(i,0);
    x(i,1) += 0.5 * l_csq * l_dt * v(i,1);
    x(i,2) += 0.5 * l_csq * l_dt * v(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::final_integrate()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  l_rmass_flag = (atom->rmass != nullptr) ? 1 : 0;

  atomKK->sync(execution_space, V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);

  v    = atomKK->k_v.view<DeviceType>();
  f    = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  if (l_rmass_flag) rmass = atomKK->k_rmass.view<DeviceType>();
  else              mass  = atomKK->k_mass.view<DeviceType>();

  // these member vars were set in initial_integrate; they have the same values
  l_dtf = static_cast<KK_FLOAT>(0.5 * update->dt * force->ftm2v);
  l_csq = static_cast<KK_FLOAT>(sqrt(gjfc3 / gjfc1));

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixGJFFinalIntegrate>(0, nlocal), *this);
  copymode = 0;
  atomKK->modified(execution_space, V_MASK);

  lv_allocated = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixGJFKokkos<DeviceType>::operator()(TagFixGJFFinalIntegrate, const int &i) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT m;
    if (l_rmass_flag) m = rmass[i];
    else              m = mass[type[i]];

    // Eq. 24f: v += csq * dtfm * f
    const KK_FLOAT dtfm = l_dtf / m;
    v(i,0) += l_csq * dtfm * f(i,0);
    v(i,1) += l_csq * dtfm * f(i,1);
    v(i,2) += l_csq * dtfm * f(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::end_of_step()
{
  // vhalf mode: swap v (Eq. 24f from this step) and lv (Eq. 24c from this step)
  // after swap: v holds the half-step velocity; lv holds 24f for next step

  const int nlocal = atom->nlocal;

  atomKK->sync(execution_space, V_MASK | MASK_MASK);

  v    = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  k_lv.sync<DeviceType>();
  d_lv = k_lv.view<DeviceType>();

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixGJFEndOfStep>(0, nlocal), *this);
  copymode = 0;

  k_lv.modify<DeviceType>();
  atomKK->modified(execution_space, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixGJFKokkos<DeviceType>::operator()(TagFixGJFEndOfStep, const int &i) const
{
  if (mask[i] & groupbit) {
    // swap v and lv
    const KK_FLOAT tmp0 = v(i,0);
    const KK_FLOAT tmp1 = v(i,1);
    const KK_FLOAT tmp2 = v(i,2);
    v(i,0) = d_lv(i,0);
    v(i,1) = d_lv(i,1);
    v(i,2) = d_lv(i,2);
    d_lv(i,0) = tmp0;
    d_lv(i,1) = tmp1;
    d_lv(i,2) = tmp2;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::grow_arrays(int nmax)
{
  memoryKK->grow_kokkos(k_lv, lv, nmax, "fix_gjf:lv");
  d_lv = k_lv.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_lv.sync_host();
  FixGJF::copy_arrays(i, j, 0);
  k_lv.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixGJFKokkos<DeviceType>::pack_exchange_item(const int &mysend, int &offset,
                                                   const bool & /*final*/) const
{
  const int i = d_exchange_sendlist(mysend);

  int m = nsend + offset;
  d_buf[mysend] = m;
  d_buf[m++] = d_lv(i,0);
  d_buf[m++] = d_lv(i,1);
  d_buf[m++] = d_lv(i,2);
  if (mysend == nsend - 1) d_count() = m;
  offset = m - nsend;

  const int j = d_copylist(mysend);
  if (j > -1) {
    d_lv(i,0) = d_lv(j,0);
    d_lv(i,1) = d_lv(j,1);
    d_lv(i,2) = d_lv(j,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixGJFKokkos<DeviceType>::pack_exchange_kokkos(
  const int &nsend, DAT::tdual_double_2d_lr &k_buf,
  DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
  ExecutionSpace space)
{
  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_exchange_sendlist.sync<DeviceType>();

  d_buf = typename AT::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0) * k_buf.extent(1));
  d_copylist = k_copylist.view<DeviceType>();
  d_exchange_sendlist = k_exchange_sendlist.view<DeviceType>();
  this->nsend = nsend;

  k_lv.template sync<DeviceType>();
  d_lv = k_lv.view<DeviceType>();
  Kokkos::deep_copy(d_count, 0);

  copymode = 1;
  FixGJFKokkosPackExchangeFunctor<DeviceType> pack_exchange_functor(this);
  Kokkos::parallel_scan(nsend, pack_exchange_functor);
  copymode = 0;

  k_buf.modify<DeviceType>();
  if (space == Host) k_buf.sync_host();
  else k_buf.sync_device();

  k_lv.template modify<DeviceType>();

  Kokkos::deep_copy(h_count, d_count);
  return h_count();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixGJFKokkos<DeviceType>::operator()(TagFixGJFUnpackExchange, const int &i) const
{
  const int index = d_indices(i);
  if (index > -1) {
    int m = d_buf[i];
    if (i >= nrecv1)
      m = nextrarecv1 + d_buf[nextrarecv1 + i - nrecv1];
    d_lv(index,0) = static_cast<KK_FLOAT>(d_buf[m++]);
    d_lv(index,1) = static_cast<KK_FLOAT>(d_buf[m++]);
    d_lv(index,2) = static_cast<KK_FLOAT>(d_buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGJFKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_double_2d_lr &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int nrecv1, int nextrarecv1, ExecutionSpace /*space*/)
{
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  d_buf = typename AT::t_double_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0) * k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  this->nrecv1 = nrecv1;
  this->nextrarecv1 = nextrarecv1;

  k_lv.template sync<DeviceType>();
  d_lv = k_lv.view<DeviceType>();
  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType, TagFixGJFUnpackExchange>(0, nrecv), *this);
  copymode = 0;
  k_lv.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixGJFKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_lv.sync_host();
  int m = FixGJF::pack_exchange(i, buf);
  k_lv.modify_host();
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixGJFKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_lv.sync_host();
  int m = FixGJF::unpack_exchange(nlocal, buf);
  k_lv.modify_host();
  return m;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixGJFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixGJFKokkos<LMPHostType>;
#endif
}
