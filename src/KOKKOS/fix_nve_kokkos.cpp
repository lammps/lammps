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

#include "fix_nve_kokkos.h"
#include <cstdio>
#include <cstring>
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
FixNVEKokkos<Space>::FixNVEKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = X_MASK | V_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVEKokkos<Space>::init()
{
  FixNVE::init();

  atomKK->k_mass.modify_host();
  DualViewHelper<Space>::sync(atomKK->k_mass);
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVEKokkos<Space>::initial_integrate(int vflag)
{
  atomKK->sync(Space,datamask_read);
  atomKK->modified(Space,datamask_modify);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  v = DualViewHelper<Space>::view(atomKK->k_v);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  mass = DualViewHelper<Space>::view(atomKK->k_mass);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (rmass.data()) {
    FixNVEKokkosInitialIntegrateFunctor<Space,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  } else {
    FixNVEKokkosInitialIntegrateFunctor<Space,0> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<Space>::initial_integrate_item(int i) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<Space>::initial_integrate_rmass_item(int i) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass[i];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVEKokkos<Space>::final_integrate()
{
  atomKK->sync(Space,V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  atomKK->modified(Space,V_MASK);

  v = DualViewHelper<Space>::view(atomKK->k_v);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  mass = DualViewHelper<Space>::view(atomKK->k_mass);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (rmass.data()) {
    FixNVEKokkosFinalIntegrateFunctor<Space,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  } else {
    FixNVEKokkosFinalIntegrateFunctor<Space,0> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  }

  // debug
  //atomKK->sync(Host,datamask_read);
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<Space>::final_integrate_item(int i) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<Space>::final_integrate_rmass_item(int i) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT dtfm = dtf / rmass[i];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixNVEKokkos<Space>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}

namespace LAMMPS_NS {
template class FixNVEKokkos<Device>;
template class FixNVEKokkos<Host>;
}

