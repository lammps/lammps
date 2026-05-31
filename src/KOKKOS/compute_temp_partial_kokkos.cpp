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

#include "compute_temp_partial_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempPartialKokkos<DeviceType>::ComputeTempPartialKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempPartial(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbias = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempPartialKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  invoked_scalar = update->ntimestep;

  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double t = 0.0;
  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  t = t_kk.t0;

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempPartialKokkos<DeviceType>::operator()(TagComputeTempPartialScalar<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    KK_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += (xflag*v(i,0)*v(i,0) + yflag*v(i,1)*v(i,1) +
                zflag*v(i,2)*v(i,2)) * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  int i;

  invoked_vector = update->ntimestep;

  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;
  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialVector<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  t[0] = t_kk.t0;
  t[1] = t_kk.t1;
  t[2] = t_kk.t2;
  t[3] = t_kk.t3;
  t[4] = t_kk.t4;
  t[5] = t_kk.t5;

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempPartialKokkos<DeviceType>::operator()(TagComputeTempPartialVector<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    KK_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += massone * xflag*v(i,0)*v(i,0);
    t_kk.t1 += massone * yflag*v(i,1)*v(i,1);
    t_kk.t2 += massone * zflag*v(i,2)*v(i,2);
    t_kk.t3 += massone * xflag*yflag*v(i,0)*v(i,1);
    t_kk.t4 += massone * xflag*zflag*v(i,0)*v(i,2);
    t_kk.t5 += massone * yflag*zflag*v(i,1)*v(i,2);
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   the bias is the velocity in the non-thermostatted dimensions
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::remove_bias_all()
{
  remove_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::remove_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    maxbias = atom->nmax;
    vbiasall = typename AT::t_kkfloat_1d_3("temp/partial/kk:vbiasall", maxbias);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialRemoveBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempPartialKokkos<DeviceType>::operator()(TagComputeTempPartialRemoveBias, const int &i) const {
  if (mask[i] & groupbit) {
    if (!xflag) { vbiasall(i,0) = v(i,0); v(i,0) = 0.0; }
    if (!yflag) { vbiasall(i,1) = v(i,1); v(i,1) = 0.0; }
    if (!zflag) { vbiasall(i,2) = v(i,2); v(i,2) = 0.0; }
  }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::restore_bias_all()
{
  restore_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::restore_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialRestoreBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempPartialKokkos<DeviceType>::operator()(TagComputeTempPartialRestoreBias, const int &i) const {
  if (mask[i] & groupbit) {
    if (!xflag) v(i,0) += vbiasall(i,0);
    if (!yflag) v(i,1) += vbiasall(i,1);
    if (!zflag) v(i,2) += vbiasall(i,2);
  }
}

/* ----------------------------------------------------------------------
   reset thermal velocity of all atoms to be consistent with bias
   called from velocity command after it creates thermal velocities
   this re-zeroes components that should stay zero
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempPartialKokkos<DeviceType>::reapply_bias_all()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempPartialReapplyBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
  atomKK->sync(Host,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempPartialKokkos<DeviceType>::operator()(TagComputeTempPartialReapplyBias, const int &i) const {
  if (mask[i] & groupbit) {
    if (!xflag) v(i,0) = 0.0;
    if (!yflag) v(i,1) = 0.0;
    if (!zflag) v(i,2) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempPartialKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempPartialKokkos<LMPHostType>;
#endif
}
