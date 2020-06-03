/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_temp_kokkos.h"
#include <mpi.h>
#include <cstring>
#include "atom_kokkos.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
ComputeTempKokkos<Space>::ComputeTempKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTemp(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;

  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
double ComputeTempKokkos<Space>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);
  DualViewHelper<Space>::sync(atomKK->k_mass);

  invoked_scalar = update->ntimestep;

  v = DualViewHelper<Space>::view(atomKK->k_v);
  if (atomKK->rmass)
    rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  else
    mass = DualViewHelper<Space>::view(atomKK->k_mass);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
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

template<ExecutionSpace Space>
template<int RMASS>
KOKKOS_INLINE_FUNCTION
void ComputeTempKokkos<Space>::operator()(TagComputeTempScalar<RMASS>, const int &i, CTEMP& t_kk) const {
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

template<ExecutionSpace Space>
void ComputeTempKokkos<Space>::compute_vector()
{
  atomKK->sync(execution_space,datamask_read);

  int i;

  invoked_vector = update->ntimestep;

  v = DualViewHelper<Space>::view(atomKK->k_v);
  if (atomKK->rmass)
    rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  else
    mass = DualViewHelper<Space>::view(atomKK->k_mass);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
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

template<ExecutionSpace Space>
template<int RMASS>
KOKKOS_INLINE_FUNCTION
void ComputeTempKokkos<Space>::operator()(TagComputeTempVector<RMASS>, const int &i, CTEMP& t_kk) const {
  if (mask[i] & groupbit) {
    KK_FLOAT massone = 0.0;
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
template class ComputeTempKokkos<Device>;
template class ComputeTempKokkos<Host>;
}

