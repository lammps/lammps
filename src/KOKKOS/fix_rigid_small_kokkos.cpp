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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_rigid_small_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "rigid_const.h"

using namespace LAMMPS_NS;
using namespace RigidConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg), FixRigidSmallBaseKokkos<DeviceType>(atom, domain)
{
  kokkosable = 1;
  forward_comm_device = 1;
  reverse_comm_device = 1;
  exchange_comm_device = 1;
  sort_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  const int nmax = atom->nmax;
  const int nlocal = atom->nlocal;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init()
{
  FixRigidSmall::init();
  atomKK->k_mass.modify_host();
  atomKK->k_mass.template sync<DeviceType>();
#ifdef LMP_KOKKOS_DEBUG_RNG
  this->rand_pool.init(random,seed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag)
{
  this->setup_base(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag)
{
  copymode = 1;
  this->template initial_integrate_base<false,false,false>(vflag);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/)
{
  if (langflag) this->apply_langevin_thermostat_base();
  if (earlyflag) this->compute_forces_and_torques_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate()
{
  copymode = 1;
  this->template final_integrate_base<false,false,false>();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor()
{
  atomKK->sync(Host, ALL_MASK);
  this->sync_host_base();
  FixRigidSmall::setup_pre_neighbor();
  atomKK->modified(Host, X_MASK | IMAGE_MASK);
  this->modify_host_base();
  atomKK->sync(Device, X_MASK | IMAGE_MASK);
  this->sync_device_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax)
{
  this->grow_arrays_base(nmax);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_body()
{
  nmax_body += DELTA_BODY;
  this->grow_body_base(nmax_body);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  this->sync_host_base();
  FixRigidSmall::copy_arrays(i, j, delflag);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_arrays(int i)
{
  this->sync_host_base();
  FixRigidSmall::set_arrays(i);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::set_molecule(int nlocalprev, tagint tagprev,
                                                    int imol, double *xgeom,
                                                    double *vcm, double *quat)
{
  this->sync_host_base();
  FixRigidSmall::set_molecule(nlocalprev, tagprev, imol, xgeom, vcm, quat);
  this->modify_host_base();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  this->sync_host_base();
  return FixRigidSmall::pack_exchange(i, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int result = FixRigidSmall::unpack_exchange(nlocal, buf);
  this->modify_host_base();
  return result;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm(int n, int *list,
                                                       double *buf, int pbc_flag, int *pbc)
{
  this->k_body.sync_host();
  this->k_bodyown.sync_host();
  return FixRigidSmall::pack_forward_comm(n, list, buf, pbc_flag, pbc);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  FixRigidSmall::unpack_forward_comm(n, first, buf);
  this->k_body.modify_host();
  this->k_bodyown.modify_host();
  this->k_body.sync_device();
  this->k_bodyown.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  this->k_body.sync_host();
  this->k_bodyown.sync_host();
  return FixRigidSmall::pack_reverse_comm(n, first, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  FixRigidSmall::unpack_reverse_comm(n, list, buf);
  this->k_body.modify_host();
  this->k_bodyown.modify_host();
  this->k_body.sync_device();
  this->k_bodyown.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::compute_scalar()
{
  return this->template compute_scalar_base<false,false,false>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum()
{
  copymode = 1;
  this->zero_momentum_base();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation()
{
  copymode = 1;
  this->zero_rotation_base();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRigidSmallKokkos<LMPHostType>;
#endif
}
