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

#include "fix_nve_asphere_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "math_extra_kokkos.h"

using namespace LAMMPS_NS;

static constexpr double INERTIA = 0.2;       // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVEAsphereKokkos<DeviceType>::FixNVEAsphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVEAsphere(lmp, narg, arg)
{
  kokkosable = 1;
  fuse_integrate_flag = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
  vatom = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereKokkos<DeviceType>::init()
{
  FixNVEAsphere::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | ANGMOM_MASK | TORQUE_MASK |
                                RMASS_MASK | ELLIPSOID_MASK | BONUS_MASK | MASK_MASK);

  bonus = avecEllipKK->k_bonus.view<DeviceType>();
  ellipsoid = atomKK->k_ellipsoid.view<DeviceType>();

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  angmom = atomKK->k_angmom.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVEAsphereKokkosInitialIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

  atomKK->modified(execution_space, X_MASK | V_MASK | ANGMOM_MASK |
                                    ELLIPSOID_MASK | BONUS_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEAsphereKokkos<DeviceType>::initial_integrate_item(const int i) const
{
  // set timestep here since dt may have changed or come via rRESPA

  const KK_FLOAT dtq = 0.5 * dtv;
  KK_FLOAT inertia[3], omega[3];
  double *shape, *quat;
  KK_FLOAT angm[3];

  if (mask(i) & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);

    // update angular momentum by 1/2 step into a local array
    angm[0] = angmom(i,0) + dtf * torque(i,0);
    angm[1] = angmom(i,1) + dtf * torque(i,1);
    angm[2] = angmom(i,2) + dtf * torque(i,2);

    // principal moments of inertia
    quat = bonus(ellipsoid(i)).quat;
    shape = bonus(ellipsoid(i)).shape;

    inertia[0] = INERTIA*rmass(i) *
                 (shape[1]*shape[1] + shape[2]*shape[2]);
    inertia[1] = INERTIA*rmass(i) *
                 (shape[0]*shape[0] + shape[2]*shape[2]);
    inertia[2] = INERTIA*rmass(i) *
                 (shape[0]*shape[0] + shape[1]*shape[1]);

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion
    MathExtraKokkos::mq_to_omega(angm, quat, inertia, omega);
    MathExtraKokkos::richardson(quat, angm, omega, inertia, dtq);

    // write back updated angular momentum
    angmom(i,0) = angm[0];
    angmom(i,1) = angm[1];
    angmom(i,2) = angm[2];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space, V_MASK | F_MASK | ANGMOM_MASK | TORQUE_MASK |
                                RMASS_MASK | MASK_MASK);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  angmom = atomKK->k_angmom.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVEAsphereKokkosFinalIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

  atomKK->modified(execution_space, V_MASK | ANGMOM_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEAsphereKokkos<DeviceType>::final_integrate_item(const int i) const
{
  if (mask(i) & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    angmom(i,0) += dtf * torque(i,0);
    angmom(i,1) += dtf * torque(i,1);
    angmom(i,2) += dtf * torque(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereKokkos<DeviceType>::fused_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | V_MASK | F_MASK | ANGMOM_MASK | TORQUE_MASK |
                                RMASS_MASK | ELLIPSOID_MASK | BONUS_MASK | MASK_MASK);

  bonus = avecEllipKK->k_bonus.view<DeviceType>();
  ellipsoid = atomKK->k_ellipsoid.view<DeviceType>();

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  angmom = atomKK->k_angmom.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVEAsphereKokkosFusedIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

  atomKK->modified(execution_space, X_MASK | V_MASK | ANGMOM_MASK |
                                    ELLIPSOID_MASK | BONUS_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEAsphereKokkos<DeviceType>::fused_integrate_item(const int i) const
{
  const KK_FLOAT dtq = 0.5 * dtv;
  KK_FLOAT inertia[3], omega[3];
  double *shape, *quat;
  KK_FLOAT angm[3];

  if (mask(i) & groupbit) {
    const KK_FLOAT dtfm = 2.0 * dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    angmom(i,0) += dtf * torque(i,0);
    angmom(i,1) += dtf * torque(i,1);
    angmom(i,2) += dtf * torque(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);

    // update angular momentum by 1/2 step into a local array
    angm[0] = angmom(i,0) + dtf * torque(i,0);
    angm[1] = angmom(i,1) + dtf * torque(i,1);
    angm[2] = angmom(i,2) + dtf * torque(i,2);

    // principal moments of inertia

    quat = bonus(ellipsoid(i)).quat;
    shape = bonus(ellipsoid(i)).shape;

    inertia[0] = INERTIA*rmass(i) *
                 (shape[1]*shape[1] + shape[2]*shape[2]);
    inertia[1] = INERTIA*rmass(i) *
                 (shape[0]*shape[0] + shape[2]*shape[2]);
    inertia[2] = INERTIA*rmass(i) *
                 (shape[0]*shape[0] + shape[1]*shape[1]);

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion

    MathExtraKokkos::mq_to_omega(angm, quat, inertia, omega);
    MathExtraKokkos::richardson(quat, angm, omega, inertia, dtq);

    // write back updated angular momentum
    angmom(i,0) = angm[0];
    angmom(i,1) = angm[1];
    angmom(i,2) = angm[2];
  }
}

namespace LAMMPS_NS {
template class FixNVEAsphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVEAsphereKokkos<LMPHostType>;
#endif
}
