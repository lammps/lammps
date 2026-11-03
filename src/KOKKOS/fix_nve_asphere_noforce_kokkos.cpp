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

#include "fix_nve_asphere_noforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "math_extra_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr double INERTIA = 0.2;    // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVEAsphereNoforceKokkos<DeviceType>::FixNVEAsphereNoforceKokkos(LAMMPS *lmp, int narg,
                                                                  char **arg) :
  FixNVEAsphereNoforce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereNoforceKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
  vatom = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereNoforceKokkos<DeviceType>::init()
{
  FixNVEAsphereNoforce::init();

  // AtomVecEllipsoidKokkos has no bonus_super array, so the superellipsoid
  // branch of the base style has no device counterpart

  if (atom->superellipsoid_flag)
    error->all(FLERR,"Cannot (yet) use superellipsoid particles with "
               "fix nve/asphere/noforce/kk");

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix nve/asphere/noforce/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEAsphereNoforceKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | V_MASK | ANGMOM_MASK | RMASS_MASK |
                                ELLIPSOID_MASK | BONUS_MASK | MASK_MASK);

  bonus = avecEllipKK->k_bonus.view<DeviceType>();
  ellipsoid = atomKK->k_ellipsoid.view<DeviceType>();

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  angmom = atomKK->k_angmom.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVEAsphereNoforceKokkosInitialIntegrateFunctor<DeviceType> functor(this);
  Kokkos::parallel_for(nlocal,functor);

  // richardson() reads angmom and writes only the bonus quaternion

  atomKK->modified(execution_space, X_MASK | ELLIPSOID_MASK | BONUS_MASK);
}

/* ----------------------------------------------------------------------
   update positions and quaternions, but leave v and angmom alone: this fix
   is fix nve/asphere without the force and torque half-kicks
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVEAsphereNoforceKokkos<DeviceType>::initial_integrate_item(const int i) const
{
  const KK_FLOAT dtv_kk = static_cast<KK_FLOAT>(dtv);
  const KK_FLOAT dtq_kk = static_cast<KK_FLOAT>(dtq);
  KK_FLOAT inertia[3], omega[3], angm[3];
  double *shape, *quat;

  if (mask(i) & groupbit) {
    x(i,0) += dtv_kk * v(i,0);
    x(i,1) += dtv_kk * v(i,1);
    x(i,2) += dtv_kk * v(i,2);

    // principal moments of inertia

    quat = bonus(ellipsoid(i)).quat;
    shape = bonus(ellipsoid(i)).shape;

    inertia[0] = static_cast<KK_FLOAT>(INERTIA*static_cast<double>(rmass(i)) *
                 (shape[1]*shape[1] + shape[2]*shape[2]));
    inertia[1] = static_cast<KK_FLOAT>(INERTIA*static_cast<double>(rmass(i)) *
                 (shape[0]*shape[0] + shape[2]*shape[2]));
    inertia[2] = static_cast<KK_FLOAT>(INERTIA*static_cast<double>(rmass(i)) *
                 (shape[0]*shape[0] + shape[1]*shape[1]));

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion

    angm[0] = angmom(i,0);
    angm[1] = angmom(i,1);
    angm[2] = angmom(i,2);

    MathExtraKokkos::mq_to_omega(angm, quat, inertia, omega);
    MathExtraKokkos::richardson(quat, angm, omega, inertia, dtq_kk);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNVEAsphereNoforceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVEAsphereNoforceKokkos<LMPHostType>;
#endif
}
