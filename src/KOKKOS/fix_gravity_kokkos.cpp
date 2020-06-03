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

#include "fix_gravity_kokkos.h"
#include "atom_masks.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "atom_kokkos.h"
#include "atom_vec.h"

using namespace LAMMPS_NS;

enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
FixGravityKokkos<Space>::FixGravityKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixGravity(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;

  datamask_read = X_MASK | F_MASK | RMASS_MASK | MASK_MASK | TYPE_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixGravityKokkos<Space>::post_force(int vflag)
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

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  if (atomKK->rmass)
    rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  else
    mass = DualViewHelper<Space>::view(atomKK->k_mass);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;

  eflag = 0;
  egrav = 0.0;

  if (atomKK->rmass) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixGravityRMass>(0,nlocal), *this, (KK_FLOAT)egrav);
  }
  else {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixGravityMass>(0,nlocal), *this, (KK_FLOAT)egrav);
  }
  copymode = 0;
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixGravityKokkos<Space>::operator()(TagFixGravityRMass, const int i, KK_FLOAT &eg) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT massone = rmass[i];
    f(i,0) += massone*xacc;
    f(i,1) += massone*yacc;
    f(i,2) += massone*zacc;
    eg -= massone * (xacc*x(i,0) + yacc*x(i,1) + zacc*x(i,2));
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixGravityKokkos<Space>::operator()(TagFixGravityMass, const int i, KK_FLOAT &eg) const
{
  if (mask[i] & groupbit) {
    KK_FLOAT massone = mass[type[i]];
    f(i,0) += massone*xacc;
    f(i,1) += massone*yacc;
    f(i,2) += massone*zacc;
    eg -= massone * (xacc*x(i,0) + yacc*x(i,1) + zacc*x(i,2));
  }
}

namespace LAMMPS_NS {
template class FixGravityKokkos<Device>;
template class FixGravityKokkos<Host>;
}
