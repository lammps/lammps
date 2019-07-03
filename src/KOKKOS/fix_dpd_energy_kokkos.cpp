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
#include <cstdio>
#include <cstring>
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
FixDPDenergyKokkos<DeviceType>::FixDPDenergyKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDPDenergy(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  pairDPDEKK = dynamic_cast<decltype(pairDPDEKK)>(pairDPDE);
  if (!pairDPDEKK)
    error->all(FLERR,"Must use pair_style dpd/fdt/energy/kk with fix dpd/energy/kk");
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixDPDenergyKokkos<DeviceType>::take_half_step()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  using AT = ArrayTypes<DeviceType>;

  atomKK->sync(execution_space, UCOND_MASK);
  typename AT::t_efloat_1d uCond = atomKK->k_uCond.view<DeviceType>();
  atomKK->sync(execution_space, UMECH_MASK);
  typename AT::t_efloat_1d uMech = atomKK->k_uMech.view<DeviceType>();

  pairDPDEKK->k_duCond.template sync<DeviceType>();
  typename AT::t_efloat_1d_const duCond = pairDPDEKK->k_duCond.template view<DeviceType>();
  pairDPDEKK->k_duMech.template sync<DeviceType>();
  typename AT::t_efloat_1d_const duMech = pairDPDEKK->k_duMech.template view<DeviceType>();

  auto dt = update->dt;

  Kokkos::parallel_for(nlocal, LAMMPS_LAMBDA(int i) {
    uCond(i) += 0.5*dt*duCond(i);
    uMech(i) += 0.5*dt*duMech(i);
  });

  atomKK->modified(execution_space, UCOND_MASK);
  atomKK->modified(execution_space, UMECH_MASK);
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixDPDenergyKokkos<DeviceType>::initial_integrate(int)
{
  take_half_step();
}

/* ---------------------------------------------------------------------- */

template <typename DeviceType>
void FixDPDenergyKokkos<DeviceType>::final_integrate()
{
  take_half_step();
}

namespace LAMMPS_NS {
template class FixDPDenergyKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDPDenergyKokkos<LMPHostType>;
#endif
}
