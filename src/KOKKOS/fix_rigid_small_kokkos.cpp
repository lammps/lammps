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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::FixRigidSmallKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRigidSmall(lmp, narg, arg),
  FixRigidBaseKokkos<DeviceType, FixRigidSmall>(atom, domain)
{
  kokkosable = 1;
  forward_comm_device = 0;
  reverse_comm_device = 0;
  exchange_comm_device = 0;
  sort_device = 0;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRigidSmallKokkos<DeviceType>::~FixRigidSmallKokkos()
{
  if (copymode) return;
}

/* ----------------------------------------------------------------------
   FIX METHODS
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::init() {
  this->init_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_constructor() {
  this->post_constructor_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup_pre_neighbor() {
  return this->setup_pre_neighbor_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::setup(int vflag) {
  this->setup_base(vflag);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::pre_neighbor() {
  this->pre_neighbor_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::initial_integrate(int vflag) { this->initial_integrate_base(vflag);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::post_force(int /*vflag*/) {
  this->post_force_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::final_integrate() {
  this->final_integrate_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_arrays(int nmax) {
  this->grow_arrays_base(nmax);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_momentum() {
  this->zero_momentum_base();
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::zero_rotation() {
  this->zero_rotation_base();
}

template<class DeviceType>
double FixRigidSmallKokkos<DeviceType>::compute_scalar() {
  return this->compute_scalar_base();
}

/* ----------------------------------------------------------------------
   FixRigidSmall PROTECTED METHODS
------------------------------------------------------------------------- */

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::grow_body() {
  this->grow_body_base();
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

/* ----------------------------------------------------------------------
   HOST COMM METHODS
------------------------------------------------------------------------- */

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_exchange(int i, double *buf) {
  this->sync_host_base();
  return FixRigidSmall::pack_exchange(i, buf);
}

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf) {
  int result = FixRigidSmall::unpack_exchange(nlocal, buf);
  this->modify_host_base();
  return result;
}

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_forward_comm(int n, int *list,
                                                          double *buf, int pbc_flag, int *pbc) {
  this->k_body.sync_host();
  this->k_bodyown.sync_host();
  return FixRigidSmall::pack_forward_comm(n, list, buf, pbc_flag, pbc);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf) {
  FixRigidSmall::unpack_forward_comm(n, first, buf);
  this->k_body.modify_host();
  this->k_bodyown.modify_host();
}

template<class DeviceType>
int FixRigidSmallKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf) {
  this->k_body.sync_host();
  this->k_bodyown.sync_host();
  return FixRigidSmall::pack_reverse_comm(n, first, buf);
}

template<class DeviceType>
void FixRigidSmallKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf) {
  FixRigidSmall::unpack_reverse_comm(n, list, buf);
  this->k_body.modify_host();
  this->k_bodyown.modify_host();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRigidSmallKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
//template class FixRigidSmallKokkos<LMPHostType>;
#endif
}
