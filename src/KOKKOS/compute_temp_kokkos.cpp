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

#include "compute_temp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempKokkos<DeviceType>::ComputeTempKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTemp(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempKokkos<DeviceType>::compute_scalar()
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
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  t = t_kk.t0; // could make this more efficient

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;

  return scalar;
}

template<class DeviceType>
template<int RMASS>
KOKKOS_INLINE_FUNCTION
void ComputeTempKokkos<DeviceType>::operator()(TagComputeTempScalar<RMASS>, const int &i, CTEMP& t_kk) const {
  if (RMASS) {
    if (mask[i] & groupbit)
      t_kk.t0 += (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2)) * rmass[i];
  } else {
    if (mask[i] & groupbit)
      t_kk.t0 += (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2)) *
        mass[type[i]];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempKokkos<DeviceType>::compute_vector()
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
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempVector<0> >(0,nlocal),*this,t_kk);
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
KOKKOS_INLINE_FUNCTION
void ComputeTempKokkos<DeviceType>::operator()(TagComputeTempVector<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    F_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += massone * v(i,0)*v(i,0);
    t_kk.t1 += massone * v(i,1)*v(i,1);
    t_kk.t2 += massone * v(i,2)*v(i,2);
    t_kk.t3 += massone * v(i,0)*v(i,1);
    t_kk.t4 += massone * v(i,0)*v(i,2);
    t_kk.t5 += massone * v(i,1)*v(i,2);
  }
}

namespace LAMMPS_NS {
template class ComputeTempKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempKokkos<LMPHostType>;
#endif
}
