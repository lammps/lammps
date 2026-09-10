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

#include "compute_temp_ramp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempRampKokkos<DeviceType>::ComputeTempRampKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempRamp(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbias = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempRampKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  invoked_scalar = update->ntimestep;

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRampScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRampScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  double t = t_kk.t0;

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
void ComputeTempRampKokkos<DeviceType>::operator()(TagComputeTempRampScalar<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    double vthermal[3];
    vthermal[0] = static_cast<double>(v(i,0));
    vthermal[1] = static_cast<double>(v(i,1));
    vthermal[2] = static_cast<double>(v(i,2));
    vthermal[v_dim] -= ramp_bias(i);
    double massone;
    if (RMASS) massone = static_cast<double>(rmass[i]);
    else massone = static_cast<double>(mass[type[i]]);
    t_kk.t0 += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
                vthermal[2]*vthermal[2]) * massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRampKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  int i;

  invoked_vector = update->ntimestep;

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRampVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRampVector<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  double t[6];
  t[0] = t_kk.t0; t[1] = t_kk.t1; t[2] = t_kk.t2;
  t[3] = t_kk.t3; t[4] = t_kk.t4; t[5] = t_kk.t5;

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRampKokkos<DeviceType>::operator()(TagComputeTempRampVector<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    double vthermal[3];
    vthermal[0] = static_cast<double>(v(i,0));
    vthermal[1] = static_cast<double>(v(i,1));
    vthermal[2] = static_cast<double>(v(i,2));
    vthermal[v_dim] -= ramp_bias(i);
    double massone;
    if (RMASS) massone = static_cast<double>(rmass[i]);
    else massone = static_cast<double>(mass[type[i]]);
    t_kk.t0 += massone * vthermal[0]*vthermal[0];
    t_kk.t1 += massone * vthermal[1]*vthermal[1];
    t_kk.t2 += massone * vthermal[2]*vthermal[2];
    t_kk.t3 += massone * vthermal[0]*vthermal[1];
    t_kk.t4 += massone * vthermal[0]*vthermal[2];
    t_kk.t5 += massone * vthermal[1]*vthermal[2];
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   the bias is the ramped velocity profile at the atom's coordinate
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRampKokkos<DeviceType>::remove_bias_all()
{
  remove_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRampKokkos<DeviceType>::remove_bias_all_kk()
{
  atomKK->sync(execution_space,X_MASK|V_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    maxbias = atom->nmax;
    vbiasall = typename AT::t_kkfloat_1d_3("temp/ramp/kk:vbiasall", maxbias);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempRampRemoveBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRampKokkos<DeviceType>::operator()(TagComputeTempRampRemoveBias, const int &i) const {
  if (mask[i] & groupbit) {
    const double vramp = ramp_bias(i);
    vbiasall(i,v_dim) = static_cast<KK_FLOAT>(vramp);
    v(i,v_dim) -= static_cast<KK_FLOAT>(vramp);
  }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRampKokkos<DeviceType>::restore_bias_all()
{
  restore_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRampKokkos<DeviceType>::restore_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempRampRestoreBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRampKokkos<DeviceType>::operator()(TagComputeTempRampRestoreBias, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,v_dim) += vbiasall(i,v_dim);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempRampKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempRampKokkos<LMPHostType>;
#endif
}
