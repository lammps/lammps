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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio@gmail.com)
------------------------------------------------------------------------- */

#include "fix_deform_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "force.h"
#include "input.h"
#include "irregular.h"
#include "kspace.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDeformKokkos<DeviceType>::FixDeformKokkos(LAMMPS *lmp, int narg, char **arg) : FixDeform(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memoryKK->create_kokkos(k_set,6,"fix_deform:set");
  d_set = k_set.template view<DeviceType>();
}

template<class DeviceType>
FixDeformKokkos<DeviceType>::~FixDeformKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_set,set);
}


template <class DeviceType>
void FixDeformKokkos<DeviceType>::init()
{
  FixDeform::init();
  memcpy((void *)k_set.h_view.data(), (void *)set, sizeof(Set)*6);
  k_set.template modify<LMPHostType>();
  k_set.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
  box flipped on previous step
  reset box tilts for flipped config and create new box in domain
  image_flip() adjusts image flags due to box shape change induced by flip
  remap() puts atoms outside the new box back into the new box
  perform irregular on atoms in lamda coords to migrate atoms to new procs
  important that image_flip comes before remap, since remap may change
    image flags to new values, making eqs in doc of Domain:image_flip incorrect
------------------------------------------------------------------------- */

template <class DeviceType>
void FixDeformKokkos<DeviceType>::pre_exchange()
{
  if (flip == 0) return;

  domainKK->yz = d_set[3].tilt_target = d_set[3].tilt_flip;
  domainKK->xz = d_set[4].tilt_target = d_set[4].tilt_flip;
  domainKK->xy = d_set[5].tilt_target = d_set[5].tilt_flip;
  domainKK->set_global_box();
  domainKK->set_local_box();

  domainKK->image_flip(flipxy, flipxz, flipyz);

  // FIXME: just replace with domainKK->remap_all(), is this correct ?
  //double **x = atom->x;
  //imageint *image = atom->image;
  //int nlocal = atom->nlocal;
  //for (int i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);
  //domain->x2lamda(atom->nlocal);
  //irregular->migrate_atoms();
  //domain->lamda2x(atom->nlocal);
  domainKK->remap_all();

  flip = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDeformKokkos<DeviceType>::end_of_step()
{
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  // set new box size for strain-based dims

  apply_strain();

  // set new box size for VOLUME dims that are linked to other dims
  // NOTE: still need to set h_rate for these dims

  apply_volume();

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);

  update_domain();

  // redo KSpace coeffs since box has changed

  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   apply strain controls
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDeformKokkos<DeviceType>::apply_strain()
{
  // for NONE, target is current box size
  // for TRATE, set target directly based on current time, also set h_rate
  // for WIGGLE, set target directly based on current time, also set h_rate
  // for VARIABLE, set target directly via variable eval, also set h_rate
  // for others except VOLUME, target is linear value between start and stop

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  for (int i = 0; i < 3; i++) {
    if (d_set[i].style == NONE) {
      d_set[i].lo_target = domainKK->boxlo[i];
      d_set[i].hi_target = domainKK->boxhi[i];
    } else if (d_set[i].style == TRATE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      double shift = 0.5 * ((d_set[i].hi_start - d_set[i].lo_start) * exp(d_set[i].rate * delt));
      d_set[i].lo_target = 0.5 * (d_set[i].lo_start + d_set[i].hi_start) - shift;
      d_set[i].hi_target = 0.5 * (d_set[i].lo_start + d_set[i].hi_start) + shift;
      h_rate[i] = d_set[i].rate * domainKK->h[i];
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (d_set[i].style == WIGGLE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      double shift = 0.5 * d_set[i].amplitude * sin(MY_2PI * delt / d_set[i].tperiod);
      d_set[i].lo_target = d_set[i].lo_start - shift;
      d_set[i].hi_target = d_set[i].hi_start + shift;
      h_rate[i] = MY_2PI / d_set[i].tperiod * d_set[i].amplitude *
        cos(MY_2PI * delt / d_set[i].tperiod);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (d_set[i].style == VARIABLE) {
      double del = input->variable->compute_equal(d_set[i].hvar);
      d_set[i].lo_target = d_set[i].lo_start - 0.5 * del;
      d_set[i].hi_target = d_set[i].hi_start + 0.5 * del;
      h_rate[i] = input->variable->compute_equal(d_set[i].hratevar);
      h_ratelo[i] = -0.5 * h_rate[i];
    } else if (d_set[i].style == FINAL || d_set[i].style == DELTA || d_set[i].style == SCALE ||
               d_set[i].style == VEL || d_set[i].style == ERATE) {
      d_set[i].lo_target = d_set[i].lo_start + delta * (d_set[i].lo_stop - d_set[i].lo_start);
      d_set[i].hi_target = d_set[i].hi_start + delta * (d_set[i].hi_stop - d_set[i].hi_start);
    }
  }

  // for triclinic, set new box shape
  // for NONE, target is current tilt
  // for TRATE, set target directly based on current time. also set h_rate
  // for WIGGLE, set target directly based on current time. also set h_rate
  // for VARIABLE, set target directly via variable eval. also set h_rate
  // for other styles, target is linear value between start and stop values

  if (triclinic) {
    for (int i = 3; i < 6; i++) {
      if (d_set[i].style == NONE) {
        if (i == 5) d_set[i].tilt_target = domainKK->xy;
        else if (i == 4) d_set[i].tilt_target = domainKK->xz;
        else if (i == 3) d_set[i].tilt_target = domainKK->yz;
      } else if (d_set[i].style == TRATE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        d_set[i].tilt_target = d_set[i].tilt_start * exp(d_set[i].rate * delt);
        h_rate[i] = d_set[i].rate * domainKK->h[i];
      } else if (d_set[i].style == WIGGLE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        d_set[i].tilt_target = d_set[i].tilt_start +
          d_set[i].amplitude * sin(MY_2PI * delt / d_set[i].tperiod);
        h_rate[i] = MY_2PI / d_set[i].tperiod * d_set[i].amplitude *
          cos(MY_2PI * delt / d_set[i].tperiod);
      } else if (d_set[i].style == VARIABLE) {
        double delta_tilt = input->variable->compute_equal(d_set[i].hvar);
        d_set[i].tilt_target = d_set[i].tilt_start + delta_tilt;
        h_rate[i] = input->variable->compute_equal(d_set[i].hratevar);
      } else {
        d_set[i].tilt_target = d_set[i].tilt_start + delta * (d_set[i].tilt_stop - d_set[i].tilt_start);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   apply volume controls
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDeformKokkos<DeviceType>::apply_volume()
{
  for (int i = 0; i < 3; i++) {
    if (d_set[i].style != VOLUME) continue;

    int dynamic1 = d_set[i].dynamic1;
    int dynamic2 = d_set[i].dynamic2;
    int fixed = d_set[i].fixed;
    double v0 = d_set[i].vol_start;
    double shift = 0.0;

    if (d_set[i].substyle == ONE_FROM_ONE) {
      shift = 0.5 * (v0 / (d_set[dynamic1].hi_target - d_set[dynamic1].lo_target) /
             (d_set[fixed].hi_start - d_set[fixed].lo_start));
    } else if (d_set[i].substyle == ONE_FROM_TWO) {
      shift = 0.5 * (v0 / (d_set[dynamic1].hi_target - d_set[dynamic1].lo_target) /
             (d_set[dynamic2].hi_target - d_set[dynamic2].lo_target));
    } else if (d_set[i].substyle == TWO_FROM_ONE) {
      shift = 0.5 * sqrt(v0 * (d_set[i].hi_start - d_set[i].lo_start) /
                 (d_set[dynamic1].hi_target - d_set[dynamic1].lo_target) /
                 (d_set[fixed].hi_start - d_set[fixed].lo_start));
    }

    h_rate[i] = (2.0 * shift / (domainKK->boxhi[i] - domainKK->boxlo[i]) - 1.0) / update->dt;
    h_ratelo[i] = -0.5 * h_rate[i];

    d_set[i].lo_target = 0.5 * (d_set[i].lo_start + d_set[i].hi_start) - shift;
    d_set[i].hi_target = 0.5 * (d_set[i].lo_start + d_set[i].hi_start) + shift;
  }
}

/* ----------------------------------------------------------------------
   Update box domain
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDeformKokkos<DeviceType>::update_domain()
{
  // tilt_target can be large positive or large negative value
  // add/subtract box lengths until tilt_target is closest to current value

  if (triclinic) {
    double *h = domainKK->h;
    for (int i = 3; i < 6; i++) {
      int idenom = 0;
      if (i == 5) idenom = 0;
      else if (i == 4) idenom = 0;
      else if (i == 3) idenom = 1;
      double denom = d_set[idenom].hi_target - d_set[idenom].lo_target;

      double current = h[i] / h[idenom];

      while (d_set[i].tilt_target / denom - current > 0.0)
        d_set[i].tilt_target -= denom;
      while (d_set[i].tilt_target / denom - current < 0.0)
        d_set[i].tilt_target += denom;
      if (fabs(d_set[i].tilt_target / denom - 1.0 - current) <
          fabs(d_set[i].tilt_target / denom - current))
        d_set[i].tilt_target -= denom;
    }
  }

  // if any tilt ratios exceed 0.5, set flip = 1 and compute new tilt values
  // do not flip in x or y if non-periodic (can tilt but not flip)
  //   this is b/c the box length would be changed (dramatically) by flip
  // if yz tilt exceeded, adjust C vector by one B vector
  // if xz tilt exceeded, adjust C vector by one A vector
  // if xy tilt exceeded, adjust B vector by one A vector
  // check yz first since it may change xz, then xz check comes after
  // flip is performed on next timestep, before reneighboring in pre-exchange()

  if (triclinic && flipflag) {
    double xprd = d_set[0].hi_target - d_set[0].lo_target;
    double yprd = d_set[1].hi_target - d_set[1].lo_target;
    double xprdinv = 1.0 / xprd;
    double yprdinv = 1.0 / yprd;
    if (d_set[3].tilt_target * yprdinv < -0.5 ||
        d_set[3].tilt_target * yprdinv > 0.5 ||
        d_set[4].tilt_target * xprdinv < -0.5 ||
        d_set[4].tilt_target * xprdinv > 0.5 ||
        d_set[5].tilt_target * xprdinv < -0.5 ||
        d_set[5].tilt_target * xprdinv > 0.5) {
      d_set[3].tilt_flip = d_set[3].tilt_target;
      d_set[4].tilt_flip = d_set[4].tilt_target;
      d_set[5].tilt_flip = d_set[5].tilt_target;

      flipxy = flipxz = flipyz = 0;

      if (domainKK->yperiodic) {
        if (d_set[3].tilt_flip * yprdinv < -0.5) {
          d_set[3].tilt_flip += yprd;
          d_set[4].tilt_flip += d_set[5].tilt_flip;
          flipyz = 1;
        } else if (d_set[3].tilt_flip * yprdinv > 0.5) {
          d_set[3].tilt_flip -= yprd;
          d_set[4].tilt_flip -= d_set[5].tilt_flip;
          flipyz = -1;
        }
      }
      if (domainKK->xperiodic) {
        if (d_set[4].tilt_flip * xprdinv < -0.5) {
          d_set[4].tilt_flip += xprd;
          flipxz = 1;
        }
        if (d_set[4].tilt_flip * xprdinv > 0.5) {
          d_set[4].tilt_flip -= xprd;
          flipxz = -1;
        }
        if (d_set[5].tilt_flip * xprdinv < -0.5) {
          d_set[5].tilt_flip += xprd;
          flipxy = 1;
        }
        if (d_set[5].tilt_flip * xprdinv > 0.5) {
          d_set[5].tilt_flip -= xprd;
          flipxy = -1;
        }
      }

      flip = 0;
      if (flipxy || flipxz || flipyz) flip = 1;
      if (flip) next_reneighbor = update->ntimestep + 1;
    }
  }

  // convert atoms and rigid bodies to lamda coords

  if (remapflag == Domain::X_REMAP) {
    atomKK->sync(execution_space, X_MASK | MASK_MASK );
    d_x = atomKK->k_x.template view<DeviceType>();
    d_mask = atomKK->k_mask.template view<DeviceType>();
    int nlocal = atomKK->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (d_mask(i) & groupbit)
        domainKK->x2lamda(&d_x(i,0), &d_x(i,0));

    for (auto &ifix : rfix)
      ifix->deform(0);
  }

  // reset global and local box to new size/shape
  // only if deform fix is controlling the dimension

  if (dimflag[0]) {
    domainKK->boxlo[0] = d_set[0].lo_target;
    domainKK->boxhi[0] = d_set[0].hi_target;
  }
  if (dimflag[1]) {
    domainKK->boxlo[1] = d_set[1].lo_target;
    domainKK->boxhi[1] = d_set[1].hi_target;
  }
  if (dimflag[2]) {
    domainKK->boxlo[2] = d_set[2].lo_target;
    domainKK->boxhi[2] = d_set[2].hi_target;
  }
  if (triclinic) {
    if (dimflag[3]) domainKK->yz = d_set[3].tilt_target;
    if (dimflag[4]) domainKK->xz = d_set[4].tilt_target;
    if (dimflag[5]) domainKK->xy = d_set[5].tilt_target;
  }

  domainKK->set_global_box();
  domainKK->set_local_box();

  // convert atoms and rigid bodies back to box coords

  if (remapflag == Domain::X_REMAP) {
    atomKK->sync(execution_space, X_MASK | MASK_MASK );
    d_x = atomKK->k_x.template view<DeviceType>();
    d_mask = atomKK->k_mask.template view<DeviceType>();
    int nlocal = atomKK->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (d_mask(i) & groupbit)
        domainKK->lamda2x(&d_x(i,0), &d_x(i,0));

    for (auto &ifix : rfix)
      ifix->deform(1);
  }
}


namespace LAMMPS_NS {
template class FixDeformKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixDeformKokkos<LMPHostType>;
#endif
}

