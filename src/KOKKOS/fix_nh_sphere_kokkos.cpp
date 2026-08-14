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

#include "fix_nh_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute.h"
#include "domain.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixNHSphereKokkos<DeviceType>::FixNHSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  if (!this->atom->omega_flag)
    this->error->all(FLERR, "Fix {} requires atom attribute omega", this->style);
  if (!this->atom->radius_flag)
    this->error->all(FLERR, "Fix {} requires atom attribute radius", this->style);

  // inertia = moment of inertia prefactor for sphere or disc

  inertia = 0.4;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "disc") == 0) {
      inertia = 0.5;
      if (this->domain->dimension != 2)
        this->error->all(FLERR, "Fix {} disc option requires 2d simulation", this->style);
    }
    iarg++;
  }
  fprintf(stderr, "flags dipole %d  dlm %d\n", this->dipole_flag, this->dlm_flag);

}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::init()
{
  // check that all group particles are finite-size

  double *radius = this->atom->radius;
  int *mask = this->atom->mask;
  int nlocal = this->atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & this->groupbit)
      if (radius[i] == 0.0)
        this->error->one(FLERR, "Fix {} requires extended particles", this->style);

  FixNHKokkos<DeviceType>::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of translational and rotational velocities
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nve_v()
{
  // translational velocity update (Kokkos kernel in FixNHKokkos)
  // also sets this->rmass and this->mask views

  FixNHKokkos<DeviceType>::nve_v();

  // rotational velocity (omega) update

  this->atomKK->sync(this->execution_space, OMEGA_MASK | TORQUE_MASK | RADIUS_MASK);

  omega_kk  = this->atomKK->k_omega.template view<DeviceType>();
  torque_kk = this->atomKK->k_torque.template view<DeviceType>();
  radius_kk = this->atomKK->k_radius.template view<DeviceType>();

  int nlocal = this->atomKK->nlocal;
  if (this->igroup == this->atomKK->firstgroup) nlocal = this->atomKK->nfirst;

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nve_v_omega>(0, nlocal), *this);
  this->copymode = 0;

  this->atomKK->modified(this->execution_space, OMEGA_MASK);
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nve_v_omega, const int &i) const
{
  if (this->mask(i) & this->groupbit) {
    const KK_FLOAT dtfrotate  = this->dtf / inertia;
    const KK_FLOAT dtirotate  = dtfrotate / (radius_kk(i) * radius_kk(i) * this->rmass(i));
    omega_kk(i, 0) += dtirotate * torque_kk(i, 0);
    omega_kk(i, 1) += dtirotate * torque_kk(i, 1);
    omega_kk(i, 2) += dtirotate * torque_kk(i, 2);
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions, and dipole (if requested)
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nve_x()
{
  // position update (Kokkos kernel in FixNHKokkos)
  // also sets this->x, this->v, this->mask views

  FixNHKokkos<DeviceType>::nve_x();

  // dipole orientation update (simple cross-product integrator only)

  if (this->dipole_flag) {
    if (this->dlm_flag)
      this->error->all(FLERR, "Fix {}: DLM dipole integrator not supported in /kk mode", this->style);

    this->atomKK->sync(this->execution_space, MU_MASK | OMEGA_MASK | MASK_MASK);
    mu_kk    = this->atomKK->k_mu.template view<DeviceType>();
    omega_kk = this->atomKK->k_omega.template view<DeviceType>();

    int nlocal = this->atomKK->nlocal;
    if (this->igroup == this->atomKK->firstgroup) nlocal = this->atomKK->nfirst;

    this->copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nve_x_dipole>(0, nlocal),
                         *this);
    this->copymode = 0;

    this->atomKK->modified(this->execution_space, MU_MASK);
  }
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nve_x_dipole, const int &i) const
{
  if (this->mask(i) & this->groupbit && mu_kk(i, 3) > 0.0) {
    const KK_FLOAT g0  = mu_kk(i, 0) + this->dtv * (omega_kk(i, 1) * mu_kk(i, 2) - omega_kk(i, 2) * mu_kk(i, 1));
    const KK_FLOAT g1  = mu_kk(i, 1) + this->dtv * (omega_kk(i, 2) * mu_kk(i, 0) - omega_kk(i, 0) * mu_kk(i, 2));
    const KK_FLOAT g2  = mu_kk(i, 2) + this->dtv * (omega_kk(i, 0) * mu_kk(i, 1) - omega_kk(i, 1) * mu_kk(i, 0));
    const KK_FLOAT msq = g0 * g0 + g1 * g1 + g2 * g2;
    const KK_FLOAT scale = mu_kk(i, 3) / Kokkos::sqrt(msq);
    mu_kk(i, 0) = g0 * scale;
    mu_kk(i, 1) = g1 * scale;
    mu_kk(i, 2) = g2 * scale;
  }
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of translational and rotational velocities
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nh_v_temp()
{
  // translational velocity scaling (handles bias removal/restore)

  FixNHKokkos<DeviceType>::nh_v_temp();

  // rotational velocity scaling by same factor_eta

  this->atomKK->sync(this->execution_space, OMEGA_MASK | MASK_MASK);
  omega_kk = this->atomKK->k_omega.template view<DeviceType>();

  int nlocal = this->atomKK->nlocal;
  if (this->igroup == this->atomKK->firstgroup) nlocal = this->atomKK->nfirst;

  this->copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nh_v_temp_omega>(0, nlocal),
                       *this);
  this->copymode = 0;

  this->atomKK->modified(this->execution_space, OMEGA_MASK);
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nh_v_temp_omega, const int &i) const
{
  if (this->mask(i) & this->groupbit) {
    omega_kk(i, 0) *= this->factor_eta;
    omega_kk(i, 1) *= this->factor_eta;
    omega_kk(i, 2) *= this->factor_eta;
  }
}

namespace LAMMPS_NS {
template class FixNHSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNHSphereKokkos<LMPHostType>;
#endif
}
