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

#include "fix_nve_sphere_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
FixNVESphereKokkos<Space>::FixNVESphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;

  datamask_read = F_MASK | TORQUE_MASK | RMASS_MASK | RADIUS_MASK | MASK_MASK;
  datamask_modify = X_MASK | V_MASK | OMEGA_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVESphereKokkos<Space>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVESphereKokkos<Space>::init()
{
  FixNVESphere::init();

  if (extra == DIPOLE) {
    error->all(FLERR,"Fix nve/sphere/kk doesn't yet support dipole");
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVESphereKokkos<Space>::initial_integrate(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  v = DualViewHelper<Space>::view(atomKK->k_v);
  omega = DualViewHelper<Space>::view(atomKK->k_omega);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  torque = DualViewHelper<Space>::view(atomKK->k_torque);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  radius = DualViewHelper<Space>::view(atomKK->k_radius);

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereKokkosInitialIntegrateFunctor<Space> f(this);
  Kokkos::parallel_for(nlocal,f);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVESphereKokkos<Space>::initial_integrate_item(const int i) const
{
  const KK_FLOAT dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);

    const KK_FLOAT dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,0) += dtirotate * torque(i,0);
    omega(i,1) += dtirotate * torque(i,1);
    omega(i,2) += dtirotate * torque(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVESphereKokkos<Space>::final_integrate()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  v = DualViewHelper<Space>::view(atomKK->k_v);
  omega = DualViewHelper<Space>::view(atomKK->k_omega);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  torque = DualViewHelper<Space>::view(atomKK->k_torque);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  radius = DualViewHelper<Space>::view(atomKK->k_radius);

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereKokkosFinalIntegrateFunctor<Space> f(this);
  Kokkos::parallel_for(nlocal,f);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVESphereKokkos<Space>::final_integrate_item(const int i) const
{
  const KK_FLOAT dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    const KK_FLOAT dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,0) += dtirotate * torque(i,0);
    omega(i,1) += dtirotate * torque(i,1);
    omega(i,2) += dtirotate * torque(i,2);
  }
}

namespace LAMMPS_NS {
template class FixNVESphereKokkos<Device>;
template class FixNVESphereKokkos<Host>;
}
