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

#include "fix_press_langevin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "irregular.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double DELTAFLIP = 0.1;
static constexpr double TILTMAX = 1.5;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPressLangevinKokkos<DeviceType>::FixPressLangevinKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPressLangevin(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  // this fix owns no per-atom kernel: the GJF barostat integrates 6 piston
  // degrees of freedom on the host (the random forces in couple_beta() are
  // 6 scalars drawn on rank 0 and broadcast), the pressure is virial-only
  // (compute pressure NULL virial; EMPTY_MASK datamask), and the only per-atom
  // work -- the box remap -- runs on device via DomainKokkos (see remap() and
  // pre_exchange() below).

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass

   mirrors FixPressLangevin::remap() exactly (including the triclinic tilt
   handling and TILTMAX check) except the per-atom coordinate conversions run
   on the device: domain is a DomainKokkos instance, so x2lamda(int)/
   x2lamda(int,int) and the lamda2x variants dispatch to device kernels that
   sync X themselves.  The dilate-partial path uses the device group variant
   instead of a host per-atom loop.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::remap()
{
  int nlocal = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap)
    domain->x2lamda(nlocal);
  else
    domain->x2lamda(nlocal, groupbit);

  for (auto &ifix : rfix) ifix->deform(0);

  // reset global and local box to new size/shape

  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      double oldlo = domain->boxlo[i];
      double oldhi = domain->boxhi[i];
      double ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo - ctr) * dilation[i] + ctr;
      domain->boxhi[i] = (oldhi - ctr) * dilation[i] + ctr;
    }
  }

  if (p_flag[3]) domain->xy += dilation[3];
  if (p_flag[4]) domain->xz += dilation[4];
  if (p_flag[5]) domain->yz += dilation[5];

  if (domain->yz < -TILTMAX * domain->yprd || domain->yz > TILTMAX * domain->yprd ||
      domain->xz < -TILTMAX * domain->xprd || domain->xz > TILTMAX * domain->xprd ||
      domain->xy < -TILTMAX * domain->xprd || domain->xy > TILTMAX * domain->xprd)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix {} has tilted box too far in one step - periodic cell is too far from "
               "equilibrium state",
               style);

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap)
    domain->lamda2x(nlocal);
  else
    domain->lamda2x(nlocal, groupbit);

  for (auto &ifix : rfix) ifix->deform(1);
}

/* ----------------------------------------------------------------------
   flip box if tilt exceeds limits, then migrate atoms to new procs

   mirrors FixPressLangevin::pre_exchange() but runs the per-atom remap on the
   device via DomainKokkos::remap_all()/image_flip()/x2lamda()/lamda2x().  Only
   the irregular atom migration is intrinsically host-side, so X is synced to
   the host around that call (and back).  This hook is only active when tilt
   flips can occur (triclinic / pre-tilted boxes) and fires rarely.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::pre_exchange()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;

  // flip is only triggered when tilt exceeds 0.5 by DELTAFLIP
  // this avoids immediate re-flipping due to tilt oscillations

  double xtiltmax = (0.5 + DELTAFLIP) * xprd;
  double ytiltmax = (0.5 + DELTAFLIP) * yprd;

  int flipxy, flipxz, flipyz;
  flipxy = flipxz = flipyz = 0;

  if (domain->yperiodic) {
    if (domain->yz < -ytiltmax) {
      domain->yz += yprd;
      domain->xz += domain->xy;
      flipyz = 1;
    } else if (domain->yz >= ytiltmax) {
      domain->yz -= yprd;
      domain->xz -= domain->xy;
      flipyz = -1;
    }
  }

  if (domain->xperiodic) {
    if (domain->xz < -xtiltmax) {
      domain->xz += xprd;
      flipxz = 1;
    } else if (domain->xz >= xtiltmax) {
      domain->xz -= xprd;
      flipxz = -1;
    }
    if (domain->xy < -xtiltmax) {
      domain->xy += xprd;
      flipxy = 1;
    } else if (domain->xy >= xtiltmax) {
      domain->xy -= xprd;
      flipxy = -1;
    }
  }

  int flip = 0;
  if (flipxy || flipxz || flipyz) flip = 1;

  if (flip) {
    domain->set_global_box();
    domain->set_local_box();

    domainKK->image_flip(flipxy, flipxz, flipyz);

    domainKK->remap_all();

    domainKK->x2lamda(atom->nlocal);
    atomKK->sync(Host, ALL_MASK);
    irregular->migrate_atoms();
    atomKK->modified(Host, ALL_MASK);
    domainKK->lamda2x(atom->nlocal);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPressLangevinKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPressLangevinKokkos<LMPHostType>;
#endif
}
