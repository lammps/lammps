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

#include "fix_momentum_kokkos.h"
#include <cstdlib>
#include <cstring>
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "domain_kokkos.h"
#include "group.h"
#include "error.h"
#include "force.h"
#include "kokkos_few.h"

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
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
static double get_kinetic_energy(
    AtomKokkos* atomKK,
    MPI_Comm world,
    int groupbit,
    int nlocal,
    typename ArrayTypes<DeviceType>::t_v_array_randomread v,
    typename ArrayTypes<DeviceType>::t_int_1d_randomread mask)
{
  using AT = ArrayTypes<DeviceType>;
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  double ke=0.0;
  if (atomKK->rmass) {
    atomKK->sync(execution_space, RMASS_MASK);
    typename AT::t_float_1d_randomread rmass = atomKK->k_rmass.view<DeviceType>();
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, double& update) {
      if (mask(i) & groupbit)
        update += rmass(i) *
          (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    }, ke);
  } else {
    // D.I. : why is there no MASS_MASK ?
    atomKK->sync(execution_space, TYPE_MASK);
    typename AT::t_int_1d_randomread type = atomKK->k_type.view<DeviceType>();
    typename AT::t_float_1d_randomread mass = atomKK->k_mass.view<DeviceType>();
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, double& update) {
      if (mask(i) & groupbit)
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

  if (rescale) {
    ekin_old = get_kinetic_energy<DeviceType>(atomKK, world, groupbit, nlocal, v, mask);
  }

  auto groupbit2 = groupbit;
  if (linear) {
    /* this is needed because Group is not Kokkos-aware ! */
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space,
        V_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);
    Few<double, 3> tmpvcm;
    group->vcm(igroup,masstotal,&tmpvcm[0]);
    const Few<double, 3> vcm(tmpvcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    auto xflag2 = xflag;
    auto yflag2 = yflag;
    auto zflag2 = zflag;
    const Few<double,3> &vcm_ref = vcm;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask(i) & groupbit2) {
        if (xflag2) v(i,0) -= vcm_ref[0];
        if (yflag2) v(i,1) -= vcm_ref[1];
        if (zflag2) v(i,2) -= vcm_ref[2];
      }
    });
    atomKK->modified(execution_space, V_MASK);
  }

  if (angular) {
    Few<double, 3> tmpxcm, tmpangmom, tmpomega;
    double inertia[3][3];
    /* syncs for each Kokkos-unaware Group method */
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space,
        X_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->xcm(igroup,masstotal,&tmpxcm[0]);
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space,
        X_MASK | V_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->angmom(igroup,&tmpxcm[0],&tmpangmom[0]);
    atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space,
        X_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->inertia(igroup,&tmpxcm[0],inertia);
    group->omega(&tmpangmom[0],inertia,&tmpomega[0]);
    const Few<double, 3> xcm(tmpxcm), angmom(tmpangmom), omega(tmpomega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    atomKK->sync(execution_space, X_MASK | IMAGE_MASK);
    typename AT::t_x_array_randomread x = atomKK->k_x.view<DeviceType>();
    typename AT::t_imageint_1d_randomread image = atomKK->k_image.view<DeviceType>();
    int nlocal = atom->nlocal;

    auto prd = Few<double,3>(domain->prd);
    auto h = Few<double,6>(domain->h);
    auto triclinic = domain->triclinic;

    const Few<double,3> &xcm_ref = xcm;
    const Few<double,3> &omega_ref = omega;

    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask[i] & groupbit2) {
        Few<double,3> x_i;
        x_i[0] = x(i,0);
        x_i[1] = x(i,1);
        x_i[2] = x(i,2);
        auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
        auto dx = unwrap[0] - xcm_ref[0];
        auto dy = unwrap[1] - xcm_ref[1];
        auto dz = unwrap[2] - xcm_ref[2];
        v(i,0) -= omega_ref[1]*dz - omega_ref[2]*dy;
        v(i,1) -= omega_ref[2]*dx - omega_ref[0]*dz;
        v(i,2) -= omega_ref[0]*dy - omega_ref[1]*dx;
      }
    });
    atomKK->modified(execution_space, V_MASK);
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    ekin_new = get_kinetic_energy<DeviceType>(atomKK, world, groupbit, nlocal, v, mask);

    double factor = 1.0;
    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask(i) & groupbit2) {
        v(i,0) *= factor;
        v(i,1) *= factor;
        v(i,2) *= factor;
      }
    });
    atomKK->modified(execution_space, V_MASK);
  }
}

namespace LAMMPS_NS {
template class FixMomentumKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixMomentumKokkos<LMPHostType>;
#endif
}
