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

#include <stdlib.h>
#include <string.h>
#include "fix_momentum_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contributing author: Dan Ibanez (SNL)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMomentumKokkos<DeviceType>::FixMomentumKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMomentum(lmp, narg, arg)
{
#ifdef KOKKOS_HAVE_CUDA
  if (angular) error->all(FLERR, "Kokkos+CUDA fix momentum doesn't support angular");
#endif
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMomentumKokkos<DeviceType>::init()
{
//this is done in FixMomentum::init(), but will need to change
//once Group is ported to Kokkos
//masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixMomentumKokkos<DeviceType>::get_kinetic_energy(
    int nlocal,
    typename AT::t_v_array_randomread v,
    typename AT::t_int_1d_randomread mask)
{
  double ke=0.0;
  auto groupbit2 = groupbit;
  // D.I. : does this atom->rmass check make sense in Kokkos mode ?
  if (atom->rmass) {
    atomKK->sync(execution_space, RMASS_MASK);
    typename AT::t_float_1d_randomread rmass = atomKK->k_rmass.view<DeviceType>();
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, double& update) {
      if (mask(i) & groupbit2)
        update += rmass(i) *
          (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    }, ke);
  } else {
    // D.I. : why is there no MASS_MASK ?
    atomKK->sync(execution_space, TYPE_MASK);
    typename AT::t_int_1d_randomread type = atomKK->k_type.view<DeviceType>();
    typename AT::t_float_1d_randomread mass = atomKK->k_mass.view<DeviceType>();
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, double& update) {
      if (mask(i) & groupbit2)
        update += mass(type(i)) *
          (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    }, ke);
  }
  double ke_total;
  MPI_Allreduce(&ke,&ke_total,1,MPI_DOUBLE,MPI_SUM,world);
  return ke_total;
}

template<class DeviceType>
void FixMomentumKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(execution_space, V_MASK | MASK_MASK);

  typename AT::t_v_array v = atomKK->k_v.view<DeviceType>();
  typename AT::t_int_1d_randomread mask = atomKK->k_mask.view<DeviceType>();

  const int nlocal = atom->nlocal;
  double ekin_old,ekin_new;
  ekin_old = ekin_new = 0.0;

  if (dynamic)
    masstotal = group->mass(igroup); // change once Group is ported to Kokkos

  // do nothing if group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;

  // compute kinetic energy before momentum removal, if needed

  if (rescale) ekin_old = get_kinetic_energy(nlocal, v, mask);

  auto groupbit2 = groupbit;
  if (linear) {
    double vcm[3];
    group->vcm(igroup,masstotal,vcm); // will need to change when Group has Kokkos

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set
    
    auto xflag2 = xflag;
    auto yflag2 = yflag;
    auto zflag2 = zflag;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask(i) & groupbit2) {
        if (xflag2) v(i,0) -= vcm[0];
        if (yflag2) v(i,1) -= vcm[1];
        if (zflag2) v(i,2) -= vcm[2];
      }
    });
  }

#ifndef KOKKOS_HAVE_CUDA
  if (angular) {
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    group->xcm(igroup,masstotal,xcm); // change when Group has Kokkos ?
    group->angmom(igroup,xcm,angmom);
    group->inertia(igroup,xcm,inertia);
    group->omega(angmom,inertia,omega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    atomKK->sync(execution_space, X_MASK | IMAGE_MASK);
    typename AT::t_x_array_randomread x = atomKK->k_x.view<DeviceType>();
    typename AT::t_imageint_1d_randomread image = atomKK->k_image.view<DeviceType>();
    int nlocal = atom->nlocal;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask[i] & groupbit2) {
        double dx,dy,dz;
        double unwrap[3];
        domain->unmap(&x(i,0),image(i),unwrap); // this will not work in CUDA
        dx = unwrap[0] - xcm[0];
        dy = unwrap[1] - xcm[1];
        dz = unwrap[2] - xcm[2];
        v(i,0) -= omega[1]*dz - omega[2]*dy;
        v(i,1) -= omega[2]*dx - omega[0]*dz;
        v(i,2) -= omega[0]*dy - omega[1]*dx;
      }
    });
  }
#endif

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    ekin_new = get_kinetic_energy(nlocal, v, mask);

    double factor = 1.0;
    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask(i) & groupbit2) {
        v(i,0) *= factor;
        v(i,1) *= factor;
        v(i,2) *= factor;
      }
    });
  }
}

namespace LAMMPS_NS {
template class FixMomentumKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixMomentumKokkos<LMPHostType>;
#endif
}

