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

#include "fix_press_berendsen_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute.h"
#include "domain.h"
#include "force.h"
#include "kokkos.h"
#include "kspace.h"
#include "modify.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { ISO, ANISO };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPressBerendsenKokkos<DeviceType>::FixPressBerendsenKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPressBerendsen(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  // this fix performs no per-atom kernel of its own: the kinetic energy is
  // reduced on-device by the temperature compute (temp/kk), and the box remap
  // runs on device via DomainKokkos::x2lamda / lamda2x (see remap() below).

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  // the base class created the internal "temp" compute relying on the -sf kk
  // suffix.  When this style is requested with an explicit /kk suffix but
  // without -sf kk, that compute is not a KOKKOS style and would force a
  // host/device sync every step.  Recreate it as temp/kk in that case.

  if (tflag) {
    Compute *c = modify->get_compute_by_id(id_temp);
    if (c && !c->kokkosable) {
      modify->delete_compute(id_temp);
      modify->add_compute(fmt::format("{} all temp/kk", id_temp));
    }
  }
}

/* ----------------------------------------------------------------------
   warn if the temperature compute is not a KOKKOS style (e.g. set via
   fix_modify temp to a non-kk compute): correct but forces per-step syncs.
   (the pressure compute has no /kk variant, so it is not checked)
------------------------------------------------------------------------- */

template<class DeviceType>
void FixPressBerendsenKokkos<DeviceType>::init()
{
  FixPressBerendsen::init();
  KokkosLMP::warn_nonkokkos_compute(lmp, style, temperature, "temperature");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressBerendsenKokkos<DeviceType>::end_of_step()
{
  // compute new T,P
  // a kokkosable temperature (temp/kk) syncs the velocities to the device
  // itself; a non-kokkosable temperature (e.g. via fix_modify temp) needs an
  // explicit host sync so it does not read a stale velocity snapshot.

  if (pstyle == ISO) {
    if (temperature->kokkosable)
      temperature->compute_scalar();
    else {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->compute_scalar();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
    pressure->compute_scalar();
  } else {
    if (temperature->kokkosable)
      temperature->compute_vector();
    else {
      atomKK->sync(temperature->execution_space, temperature->datamask_read);
      temperature->compute_vector();
      atomKK->modified(temperature->execution_space, temperature->datamask_modify);
      atomKK->sync(execution_space, temperature->datamask_modify);
    }
    pressure->compute_vector();
  }
  couple();

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i] - p_start[i]);
      dilation[i] = pow(1.0 - update->dt / p_period[i] * (p_target[i] - p_current[i]) / bulkmodulus,
                        1.0 / 3.0);
    }
  }

  // remap simulation box and atoms (on device)
  // redo KSpace coeffs since volume has changed

  remap();
  if (kspace_flag) force->kspace->setup();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep + 1);
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass

   identical to FixPressBerendsen::remap() except the per-atom coordinate
   conversions run on the device: domain is a DomainKokkos instance under
   the KOKKOS package, so the virtual x2lamda(int)/x2lamda(int,int) and
   lamda2x variants dispatch to device kernels that sync X themselves.
------------------------------------------------------------------------- */

template<class DeviceType>
void FixPressBerendsenKokkos<DeviceType>::remap()
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

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap)
    domain->lamda2x(nlocal);
  else
    domain->lamda2x(nlocal, groupbit);

  for (auto &ifix : rfix) ifix->deform(1);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPressBerendsenKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPressBerendsenKokkos<LMPHostType>;
#endif
}
