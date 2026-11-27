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

#include "fix_damping_cundall_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "input.h"
#include "memory_kokkos.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

enum { NONE, TYPE, VARIABLE };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDampingCundallKokkos<DeviceType>::FixDampingCundallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDampingCundall(lmp, narg, arg), maxatom_scaleval(0)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | MASK_MASK | TYPE_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDampingCundallKokkos<DeviceType>::~FixDampingCundallKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDampingCundallKokkos<DeviceType>::init()
{
  FixDampingCundall::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix damping/cundall/kk");

  // per-type scale factors are fixed once the style is defined

  if (scalestyle == TYPE) {
    const int ntypes = atom->ntypes;
    k_scalegamma = DAT::tdual_kkfloat_1d("damping/cundall:scalegamma",ntypes+1);
    for (int i = 1; i <= ntypes; i++)
      k_scalegamma.view_host()(i) = static_cast<KK_FLOAT>(scalegamma[i]);
    k_scalegamma.template modify<LMPHostType>();
    k_scalegamma.template sync<DeviceType>();
    d_scalegamma = k_scalegamma.template view<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDampingCundallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // an atom-style variable can only be evaluated on the host; push the
  // result into a device view for the kernel to read

  if (scalestyle == VARIABLE) {
    if (atom->nmax > maxatom_scaleval) {
      maxatom_scaleval = atom->nmax;
      memory->destroy(scaleval);
      memory->create(scaleval,maxatom_scaleval,"fix_damping/cundall:scaleval");
      k_scaleval = DAT::tdual_kkfloat_1d("damping/cundall:scaleval",maxatom_scaleval);
      d_scaleval = k_scaleval.template view<DeviceType>();
    }
    atomKK->sync(Host,ALL_MASK);
    input->variable->compute_atom(scalevar,igroup,scaleval,1,0);
    const int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      k_scaleval.view_host()(i) = static_cast<KK_FLOAT>(scaleval[i]);
    k_scaleval.template modify<LMPHostType>();
    k_scaleval.template sync<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  d_f = atomKK->k_f.view<DeviceType>();
  d_torque = atomKK->k_torque.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_omega = atomKK->k_omega.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();

  const int nlocal = atom->nlocal;

  copymode = 1;
  if (scalestyle == TYPE)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixDampingCundall<TYPE>>(0,nlocal),*this);
  else if (scalestyle == VARIABLE)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixDampingCundall<VARIABLE>>(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixDampingCundall<NONE>>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
}

/* ----------------------------------------------------------------------
   apply damping force/torque to finite-size atoms in group
   add a fraction of the current force/torque if work is negative
   subtract a fraction of the current force/torque if work is positive
   applied over each component independently (non-objective)
   magnitude depends on atom type

   see, e.g. Yade-DEM implementation of NewtonIntegrator::cundallDamp1st()
   gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/NewtonIntegrator.cpp
------------------------------------------------------------------------- */

template<class DeviceType>
template<int SCALE>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixDampingCundallKokkos<DeviceType>::operator()(TagFixDampingCundall<SCALE>,
                                                     const int &i) const
{
  if (d_mask[i] & groupbit) {
    const KK_FLOAT gamma_lin_kk = static_cast<KK_FLOAT>(gamma_lin);
    const KK_FLOAT gamma_ang_kk = static_cast<KK_FLOAT>(gamma_ang);

    KK_FLOAT gamma_l, gamma_a;
    if (SCALE == TYPE) {
      gamma_l = gamma_lin_kk * d_scalegamma(d_type(i));
      gamma_a = gamma_ang_kk * d_scalegamma(d_type(i));
    } else if (SCALE == VARIABLE) {
      gamma_l = gamma_lin_kk * d_scaleval(i);
      gamma_a = gamma_ang_kk * d_scaleval(i);
    } else {    // scalestyle NONE
      gamma_l = gamma_lin_kk;
      gamma_a = gamma_ang_kk;
    }

    const KK_FLOAT signf0 = (static_cast<KK_FLOAT>(d_f(i,0))*d_v(i,0) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    const KK_FLOAT signf1 = (static_cast<KK_FLOAT>(d_f(i,1))*d_v(i,1) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    const KK_FLOAT signf2 = (static_cast<KK_FLOAT>(d_f(i,2))*d_v(i,2) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    d_f(i,0) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_l*signf0);
    d_f(i,1) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_l*signf1);
    d_f(i,2) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_l*signf2);

    const KK_FLOAT signt0 = (static_cast<KK_FLOAT>(d_torque(i,0))*d_omega(i,0) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    const KK_FLOAT signt1 = (static_cast<KK_FLOAT>(d_torque(i,1))*d_omega(i,1) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    const KK_FLOAT signt2 = (static_cast<KK_FLOAT>(d_torque(i,2))*d_omega(i,2) >= static_cast<KK_FLOAT>(0.0)) ? 1.0 : -1.0;
    d_torque(i,0) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_a*signt0);
    d_torque(i,1) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_a*signt1);
    d_torque(i,2) *= static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(1.0) - gamma_a*signt2);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDampingCundallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixDampingCundallKokkos<LMPHostType>;
#endif
}
