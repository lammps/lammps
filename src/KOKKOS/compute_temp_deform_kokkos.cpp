// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Emily Kahl (Uni. of QLD, e.kahl@uq.edu.au)
------------------------------------------------------------------------- */

#include "compute_temp_deform_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempDeformKokkos<DeviceType>::ComputeTempDeformKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempDeform(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbias = 0;
}

template<class DeviceType>
ComputeTempDeformKokkos<DeviceType>::~ComputeTempDeformKokkos()
{


}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempDeformKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  invoked_scalar = update->ntimestep;

  v = atomKK->k_v.view<DeviceType>();
  x = atomKK->k_x.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double t = 0.0;
  CTEMP t_kk;

  domainKK->x2lamda(nlocal);
  h_rate = domainKK->h_rate;
  h_ratelo = domainKK->h_ratelo;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  domainKK->lamda2x(nlocal);

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
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformScalar<RMASS>, const int &i, CTEMP& t_kk) const {

  double vstream[3],vthermal[3];

  vstream[0] = h_rate[0]*x(i,0) + h_rate[5]*x(i,1) + h_rate[4]*x(i,2) + h_ratelo[0];
  vstream[1] = h_rate[1]*x(i,1) + h_rate[3]*x(i,2) + h_ratelo[1];
  vstream[2] = h_rate[2]*x(i,2) + h_ratelo[2];
  vthermal[0] = v(i,0) - vstream[0];
  vthermal[1] = v(i,1) - vstream[1];
  vthermal[2] = v(i,2) - vstream[2];
  if (RMASS) {
    if (mask[i] & groupbit)
      t_kk.t0 += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + vthermal[2]*vthermal[2]) * rmass[i];
  } else {
    if (mask[i] & groupbit)
      t_kk.t0 += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + vthermal[2]*vthermal[2]) * mass[type[i]];
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space,datamask_read);

  int i;

  invoked_vector = update->ntimestep;

  v = atomKK->k_v.view<DeviceType>();
  x = atomKK->k_x.view<DeviceType>();
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

  domainKK->x2lamda(nlocal);
  h_rate = domainKK->h_rate;
  h_ratelo = domainKK->h_ratelo;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformVector<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  domainKK->lamda2x(nlocal);

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
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformVector<RMASS>, const int &i, CTEMP& t_kk) const {

  double vstream[3],vthermal[3];

  vstream[0] = h_rate[0]*x(i,0) + h_rate[5]*x(i,1) + h_rate[4]*x(i,2) + h_ratelo[0];
  vstream[1] = h_rate[1]*x(i,1) + h_rate[3]*x(i,2) + h_ratelo[1];
  vstream[2] = h_rate[2]*x(i,2) + h_ratelo[2];
  vthermal[0] = v(i,0) - vstream[0];
  vthermal[1] = v(i,1) - vstream[1];
  vthermal[2] = v(i,2) - vstream[2];

  if (mask[i] & groupbit) {
    F_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += massone * vthermal[0]*vthermal[0];
    t_kk.t1 += massone * vthermal[1]*vthermal[1];
    t_kk.t2 += massone * vthermal[2]*vthermal[2];
    t_kk.t3 += massone * vthermal[0]*vthermal[1];
    t_kk.t4 += massone * vthermal[0]*vthermal[2];
    t_kk.t5 += massone * vthermal[1]*vthermal[2];
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::remove_bias_all()
{
  atomKK->sync(execution_space,datamask_read);
  v = atomKK->k_v.view<DeviceType>();
  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    maxbias = atom->nmax;
    vbiasall = typename ArrayTypes<DeviceType>::t_v_array("temp/deform/kk:vbiasall", maxbias);
  }

  domainKK->x2lamda(nlocal);

  h_rate = domain->h_rate;
  h_ratelo = domain->h_ratelo;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformRemoveBias >(0,nlocal),*this);
  copymode = 0;

  domainKK->lamda2x(nlocal);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformRemoveBias, const int &i) const {
  if (mask[i] & groupbit) {
    vbiasall(i,0) = h_rate[0]*x(i,0) + h_rate[5]*x(i,1) + h_rate[4]*x(i,2) + h_ratelo[0];
    vbiasall(i,1) = h_rate[1]*x(i,1) + h_rate[3]*x(i,2) + h_ratelo[1];
    vbiasall(i,2) = h_rate[2]*x(i,2) + h_ratelo[2];
    v(i,0) -= vbiasall(i,0);
    v(i,1) -= vbiasall(i,1);
    v(i,2) -= vbiasall(i,2);
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeTempDeformKokkos<DeviceType>::restore_bias_all()
{
  atomKK->sync(execution_space,datamask_read);
  v = atomKK->k_v.view<DeviceType>();
  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempDeformRestoreBias >(0,nlocal),*this);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeTempDeformKokkos<DeviceType>::operator()(TagComputeTempDeformRestoreBias, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) += vbiasall(i,0);
    v(i,1) += vbiasall(i,1);
    v(i,2) += vbiasall(i,2);
  }
}

namespace LAMMPS_NS {
template class ComputeTempDeformKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempDeformKokkos<LMPHostType>;
#endif
}
