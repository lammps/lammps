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

#include "fix_dt_reset_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "output.h"
#include "pair.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDtResetKokkos<DeviceType>::FixDtResetKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDtReset(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDtResetKokkos<DeviceType>::init()
{
  FixDtReset::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDtResetKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(execution_space, V_MASK | F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();

  int nlocal = atom->nlocal;

  double dt;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixDtResetRMass>(0,nlocal), *this, Kokkos::Min<double>(dt));
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixDtResetMass>(0,nlocal), *this, Kokkos::Min<double>(dt));
  copymode = 0;

  MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, world);

  if (minbound) dt = MAX(dt, tmin);
  if (maxbound) dt = MIN(dt, tmax);

  // if timestep didn't change, just return
  // else reset update->dt and other classes that depend on it
  // rRESPA, pair style, fixes

  if (dt == update->dt) return;

  laststep = update->ntimestep;

  // calls to other classes that need to know timestep size changed
  // similar logic is in Input::timestep()

  update->update_time();
  update->dt = dt;
  update->dt_default = 0;
  if (force->pair) force->pair->reset_dt();
  for (auto &ifix : modify->get_fix_list()) ifix->reset_dt();
  output->reset_dt();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixDtResetKokkos<DeviceType>::operator()(TagFixDtResetMass, const int &i, double &dt_min) const {

  KK_FLOAT dt, dtv, dtf, dte, dtsq;
  KK_FLOAT vsq, fsq, massinv;
  KK_FLOAT delx, dely, delz, delr;

  if (mask[i] & groupbit) {

    const KK_FLOAT xmax_kk = static_cast<KK_FLOAT>(xmax);
    const KK_FLOAT ftm2v_kk = static_cast<KK_FLOAT>(ftm2v);
    massinv = static_cast<KK_FLOAT>(1.0) / mass[type[i]];
    vsq = v(i,0) * v(i,0) + v(i,1) * v(i,1) + v(i,2) * v(i,2);
    fsq = static_cast<KK_FLOAT>(f(i,0) * f(i,0) + f(i,1) * f(i,1) + f(i,2) * f(i,2));
    dtv = dtf = dte = static_cast<KK_FLOAT>(BIG);
    if (vsq > static_cast<KK_FLOAT>(0.0)) dtv = xmax_kk / Kokkos::sqrt(vsq);
    if (fsq > static_cast<KK_FLOAT>(0.0)) dtf = Kokkos::sqrt(static_cast<KK_FLOAT>(2.0) * xmax_kk / (ftm2v_kk * Kokkos::sqrt(fsq) * massinv));
    dt = MIN(dtv, dtf);
    if ((emax > 0.0) && (fsq * vsq > static_cast<KK_FLOAT>(0.0))) {
      dte = static_cast<KK_FLOAT>(emax) / Kokkos::sqrt(fsq * vsq) / Kokkos::sqrt(ftm2v_kk * static_cast<KK_FLOAT>(mvv2e));
      dt = MIN(dt, dte);
    }
    dtsq = dt * dt;
    delx = dt * v(i,0) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,0)) * ftm2v_kk;
    dely = dt * v(i,1) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,1)) * ftm2v_kk;
    delz = dt * v(i,2) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,2)) * ftm2v_kk;
    delr = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
    if (delr > xmax_kk) dt *= xmax_kk / delr;
    dt_min = MIN(dt_min,static_cast<double>(dt));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixDtResetKokkos<DeviceType>::operator()(TagFixDtResetRMass, const int &i, double &dt_min) const {

  KK_FLOAT dt, dtv, dtf, dte, dtsq;
  KK_FLOAT vsq, fsq, massinv;
  KK_FLOAT delx, dely, delz, delr;

  if (mask[i] & groupbit) {

    const KK_FLOAT xmax_kk = static_cast<KK_FLOAT>(xmax);
    const KK_FLOAT ftm2v_kk = static_cast<KK_FLOAT>(ftm2v);
    massinv = static_cast<KK_FLOAT>(1.0) / rmass[i];
    vsq = v(i,0) * v(i,0) + v(i,1) * v(i,1) + v(i,2) * v(i,2);
    fsq = static_cast<KK_FLOAT>(f(i,0) * f(i,0) + f(i,1) * f(i,1) + f(i,2) * f(i,2));
    dtv = dtf = dte = static_cast<KK_FLOAT>(BIG);
    if (vsq > static_cast<KK_FLOAT>(0.0)) dtv = xmax_kk / Kokkos::sqrt(vsq);
    if (fsq > static_cast<KK_FLOAT>(0.0)) dtf = Kokkos::sqrt(static_cast<KK_FLOAT>(2.0) * xmax_kk / (ftm2v_kk * Kokkos::sqrt(fsq) * massinv));
    dt = MIN(dtv, dtf);
    if ((emax > 0.0) && (fsq * vsq > static_cast<KK_FLOAT>(0.0))) {
      dte = static_cast<KK_FLOAT>(emax) / Kokkos::sqrt(fsq * vsq) / Kokkos::sqrt(ftm2v_kk * static_cast<KK_FLOAT>(mvv2e));
      dt = MIN(dt, dte);
    }
    dtsq = dt * dt;
    delx = dt * v(i,0) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,0)) * ftm2v_kk;
    dely = dt * v(i,1) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,1)) * ftm2v_kk;
    delz = dt * v(i,2) + static_cast<KK_FLOAT>(0.5) * dtsq * massinv * static_cast<KK_FLOAT>(f(i,2)) * ftm2v_kk;
    delr = Kokkos::sqrt(delx * delx + dely * dely + delz * delz);
    if (delr > xmax_kk) dt *= xmax_kk / delr;
    dt_min = MIN(dt_min,static_cast<double>(dt));
  }
}

namespace LAMMPS_NS {
template class FixDtResetKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixDtResetKokkos<LMPHostType>;
#endif
}

