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

#include "fix_dpd_energy_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
FixDPDenergyKokkos<Space>::FixDPDenergyKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDPDenergy(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  pairDPDEKK = dynamic_cast<decltype(pairDPDEKK)>(pairDPDE);
  if (!pairDPDEKK)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy/kk with fix dpd/energy/kk");
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixDPDenergyKokkos<Space>::take_half_step()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  using AT = ArrayTypes<Space>;

  atomKK->sync(Space, UCOND_MASK);
  typename AT::t_float_1d uCond = DualViewHelper<Space>::view(atomKK->k_uCond);
  atomKK->sync(Space, UMECH_MASK);
  typename AT::t_float_1d uMech = DualViewHelper<Space>::view(atomKK->k_uMech);

  DualViewHelper<Space>::sync(pairDPDEKK->k_duCond);
  typename AT::t_float_1d_const duCond = DualViewHelper<Space>::view(pairDPDEKK->k_duCond);
  DualViewHelper<Space>::sync(pairDPDEKK->k_duMech);
  typename AT::t_float_1d_const duMech = DualViewHelper<Space>::view(pairDPDEKK->k_duMech);

  auto dt = update->dt;

  Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
    uCond(i) += 0.5*dt*duCond(i);
    uMech(i) += 0.5*dt*duMech(i);
  });

  atomKK->modified(Space, UCOND_MASK);
  atomKK->modified(Space, UMECH_MASK);
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixDPDenergyKokkos<Space>::initial_integrate(int)
{
  take_half_step();
}

/* ---------------------------------------------------------------------- */

template <ExecutionSpace Space>
void FixDPDenergyKokkos<Space>::final_integrate()
{
  take_half_step();
}

namespace LAMMPS_NS {
template class FixDPDenergyKokkos<Device>;
template class FixDPDenergyKokkos<Host>;
}
