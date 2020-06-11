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

template<ExecutionSpace Space>
FixMomentumKokkos<Space>::FixMomentumKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMomentum(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
static KK_FLOAT get_kinetic_energy(
    AtomKokkos* atomKK,
    MPI_Comm world,
    int groupbit,
    int nlocal,
    typename ArrayTypes<Space>::t_float_1d_3_randomread v,
    typename ArrayTypes<Space>::t_int_1d_randomread mask)
{
  typedef typename GetDeviceType<Space>::value DeviceType;
  using AT = ArrayTypes<Space>;
  auto execution_space = Space;

  KK_FLOAT ke=0.0;
  if (atomKK->rmass) {
    atomKK->sync(Space, RMASS_MASK);
    typename AT::t_float_1d_randomread rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, KK_FLOAT& update) {
      if (mask(i) & groupbit)
        update += rmass(i) *
          (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    }, (KK_FLOAT)ke);
  } else {
    // D.I. : why is there no MASS_MASK ?
    atomKK->sync(Space, TYPE_MASK);
    typename AT::t_int_1d_randomread type = DualViewHelper<Space>::view(atomKK->k_type);
    typename AT::t_float_1d_randomread mass = DualViewHelper<Space>::view(atomKK->k_mass);
    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(int i, KK_FLOAT& update) {
      if (mask(i) & groupbit)
        update += mass(type(i)) *
          (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    }, (KK_FLOAT)ke);
  }
  KK_FLOAT ke_total;
  MPI_Allreduce(&ke,&ke_total,1,MPI_DOUBLE,MPI_SUM,world);
  return ke_total;
}

template<ExecutionSpace Space>
void FixMomentumKokkos<Space>::end_of_step()
{
  atomKK->sync(Space, V_MASK | MASK_MASK);

  typename AT::t_float_1d_3 v = DualViewHelper<Space>::view(atomKK->k_v);
  typename AT::t_int_1d_randomread mask = DualViewHelper<Space>::view(atomKK->k_mask);

  const int nlocal = atom->nlocal;
  KK_FLOAT ekin_old,ekin_new;
  ekin_old = ekin_new = 0.0;

  if (dynamic)
    masstotal = group->mass(igroup); // change once Group is ported to Kokkos

  // do nothing if group is empty, i.e. mass is zero;

  if (masstotal == 0.0) return;

  // compute kinetic energy before momentum removal, if needed

  if (rescale) {
    ekin_old = get_kinetic_energy<Space>(atomKK, world, groupbit, nlocal, v, mask);
  }

  auto groupbit2 = groupbit;
  if (linear) {
    /* this is needed because Group is not Kokkos-aware ! */
    atomKK->sync(Host,
        V_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);
    Few<double, 3> tmpvcm;
    group->vcm(igroup,masstotal,&tmpvcm[0]);
    const Few<double, 3> vcm(tmpvcm);

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
    atomKK->modified(Space, V_MASK);
  }

  if (angular) {
    Few<double, 3> tmpxcm, tmpangmom, tmpomega;
    double inertia[3][3];
    /* syncs for each Kokkos-unaware Group method */
    atomKK->sync(Host,
        X_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->xcm(igroup,masstotal,&tmpxcm[0]);
    atomKK->sync(Host,
        X_MASK | V_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->angmom(igroup,&tmpxcm[0],&tmpangmom[0]);
    atomKK->sync(Host,
        X_MASK | MASK_MASK | TYPE_MASK | IMAGE_MASK | RMASS_MASK);
    group->inertia(igroup,&tmpxcm[0],inertia);
    group->omega(&tmpangmom[0],inertia,&tmpomega[0]);
    const Few<double, 3> xcm(tmpxcm), angmom(tmpangmom), omega(tmpomega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    atomKK->sync(Space, X_MASK | IMAGE_MASK);
    typename AT::t_float_1d_3_lr_randomread x = DualViewHelper<Space>::view(atomKK->k_x);
    typename AT::t_imageint_1d_randomread image = DualViewHelper<Space>::view(atomKK->k_image);
    int nlocal = atom->nlocal;

    auto prd = Few<double,3>(domain->prd);
    auto h = Few<double,6>(domain->h);
    auto triclinic = domain->triclinic;
    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask[i] & groupbit2) {
        Few<double,3> x_i;
        x_i[0] = x(i,0);
        x_i[1] = x(i,1);
        x_i[2] = x(i,2);
        auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
        auto dx = unwrap[0] - xcm[0];
        auto dy = unwrap[1] - xcm[1];
        auto dz = unwrap[2] - xcm[2];
        v(i,0) -= omega[1]*dz - omega[2]*dy;
        v(i,1) -= omega[2]*dx - omega[0]*dz;
        v(i,2) -= omega[0]*dy - omega[1]*dx;
      }
    });
    atomKK->modified(Space, V_MASK);
  }

  // compute kinetic energy after momentum removal, if needed

  if (rescale) {

    ekin_new = get_kinetic_energy<Space>(atomKK, world, groupbit, nlocal, v, mask);

    KK_FLOAT factor = 1.0;
    if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
    Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
      if (mask(i) & groupbit2) {
        v(i,0) *= factor;
        v(i,1) *= factor;
        v(i,2) *= factor;
      }
    });
    atomKK->modified(Space, V_MASK);
  }
}

namespace LAMMPS_NS {
template class FixMomentumKokkos<Device>;
template class FixMomentumKokkos<Host>;
}

