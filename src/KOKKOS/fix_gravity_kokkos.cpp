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

#include "fix_gravity_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixGravityKokkos<DeviceType>::FixGravityKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixGravity(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | F_MASK | RMASS_MASK | MASK_MASK | TYPE_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGravityKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // update gravity due to variables

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) vert = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) phi = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) theta = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) xdir = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) ydir = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zdir = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    set_acceleration();
  }

  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;

  eflag = 0;
  egrav = 0.0;

  if (atomKK->rmass) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixGravityRMass>(0,nlocal), *this, egrav);
  }
  else {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixGravityMass>(0,nlocal), *this, egrav);
  }
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixGravityKokkos<DeviceType>::operator()(TagFixGravityRMass, const int i, double &eg) const
{
  if (mask[i] & groupbit) {
    double massone = rmass[i];
    f(i,0) += massone*xacc;
    f(i,1) += massone*yacc;
    f(i,2) += massone*zacc;
    eg -= massone * (xacc*x(i,0) + yacc*x(i,1) + zacc*x(i,2));
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixGravityKokkos<DeviceType>::operator()(TagFixGravityMass, const int i, double &eg) const
{
  if (mask[i] & groupbit) {
    double massone = mass[type[i]];
    f(i,0) += massone*xacc;
    f(i,1) += massone*yacc;
    f(i,2) += massone*zacc;
    eg -= massone * (xacc*x(i,0) + yacc*x(i,1) + zacc*x(i,2));
  }
}

namespace LAMMPS_NS {
template class FixGravityKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixGravityKokkos<LMPHostType>;
#endif
}
