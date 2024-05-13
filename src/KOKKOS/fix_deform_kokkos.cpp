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
   Contributing author: Pieter in 't Veld (SNL)
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
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE=0,FINAL,DELTA,SCALE,VEL,ERATE,TRATE,VOLUME,WIGGLE,VARIABLE};
enum{ONE_FROM_ONE,ONE_FROM_TWO,TWO_FROM_ONE};

/* ---------------------------------------------------------------------- */

FixDeformKokkos::FixDeformKokkos(LAMMPS *lmp, int narg, char **arg) : FixDeform(lmp, narg, arg)
{
  kokkosable = 1;
  domainKK = (DomainKokkos *) domain;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
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

void FixDeformKokkos::pre_exchange()
{
  if (flip == 0) return;

  domain->yz = set[3].tilt_target = set[3].tilt_flip;
  domain->xz = set[4].tilt_target = set[4].tilt_flip;
  domain->xy = set[5].tilt_target = set[5].tilt_flip;
  domain->set_global_box();
  domain->set_local_box();

  domainKK->image_flip(flipxy,flipxz,flipyz);

  domainKK->remap_all();

  domainKK->x2lamda(atom->nlocal);
  atomKK->sync(Host,ALL_MASK);
  irregular->migrate_atoms();
  atomKK->modified(Host,ALL_MASK);
  domainKK->lamda2x(atom->nlocal);

  flip = 0;
}

/* ---------------------------------------------------------------------- */

void FixDeformKokkos::end_of_step()
{
  int i;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  // set new box size
  // for NONE, target is current box size
  // for TRATE, set target directly based on current time, also set h_rate
  // for WIGGLE, set target directly based on current time, also set h_rate
  // for VARIABLE, set target directly via variable eval, also set h_rate
  // for others except VOLUME, target is linear value between start and stop

  for (i = 0; i < 3; i++) {
    if (set[i].style == NONE) {
      set[i].lo_target = domain->boxlo[i];
      set[i].hi_target = domain->boxhi[i];
    } else if (set[i].style == TRATE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      set[i].lo_target = 0.5*(set[i].lo_start+set[i].hi_start) -
        0.5*((set[i].hi_start-set[i].lo_start) * exp(set[i].rate*delt));
      set[i].hi_target = 0.5*(set[i].lo_start+set[i].hi_start) +
        0.5*((set[i].hi_start-set[i].lo_start) * exp(set[i].rate*delt));
      h_rate[i] = set[i].rate * domain->h[i];
      h_ratelo[i] = -0.5*h_rate[i];
    } else if (set[i].style == WIGGLE) {
      double delt = (update->ntimestep - update->beginstep) * update->dt;
      set[i].lo_target = set[i].lo_start -
        0.5*set[i].amplitude * sin(MY_2PI*delt/set[i].tperiod);
      set[i].hi_target = set[i].hi_start +
        0.5*set[i].amplitude * sin(MY_2PI*delt/set[i].tperiod);
      h_rate[i] = MY_2PI/set[i].tperiod * set[i].amplitude *
        cos(MY_2PI*delt/set[i].tperiod);
      h_ratelo[i] = -0.5*h_rate[i];
    } else if (set[i].style == VARIABLE) {
      double del = input->variable->compute_equal(set[i].hvar);
      set[i].lo_target = set[i].lo_start - 0.5*del;
      set[i].hi_target = set[i].hi_start + 0.5*del;
      h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      h_ratelo[i] = -0.5*h_rate[i];
    } else if (set[i].style != VOLUME) {
      set[i].lo_target = set[i].lo_start +
        delta*(set[i].lo_stop - set[i].lo_start);
      set[i].hi_target = set[i].hi_start +
        delta*(set[i].hi_stop - set[i].hi_start);
    }
  }

  // set new box size for VOLUME dims that are linked to other dims
  // NOTE: still need to set h_rate for these dims

  for (i = 0; i < 3; i++) {
    if (set[i].style != VOLUME) continue;

    if (set[i].substyle == ONE_FROM_ONE) {
      set[i].lo_target = 0.5*(set[i].lo_start+set[i].hi_start) -
        0.5*(set[i].vol_start /
             (set[set[i].dynamic1].hi_target -
              set[set[i].dynamic1].lo_target) /
             (set[set[i].fixed].hi_start-set[set[i].fixed].lo_start));
      set[i].hi_target = 0.5*(set[i].lo_start+set[i].hi_start) +
        0.5*(set[i].vol_start /
             (set[set[i].dynamic1].hi_target -
              set[set[i].dynamic1].lo_target) /
             (set[set[i].fixed].hi_start-set[set[i].fixed].lo_start));

    } else if (set[i].substyle == ONE_FROM_TWO) {
      set[i].lo_target = 0.5*(set[i].lo_start+set[i].hi_start) -
        0.5*(set[i].vol_start /
             (set[set[i].dynamic1].hi_target -
              set[set[i].dynamic1].lo_target) /
             (set[set[i].dynamic2].hi_target -
              set[set[i].dynamic2].lo_target));
      set[i].hi_target = 0.5*(set[i].lo_start+set[i].hi_start) +
        0.5*(set[i].vol_start /
             (set[set[i].dynamic1].hi_target -
              set[set[i].dynamic1].lo_target) /
             (set[set[i].dynamic2].hi_target -
              set[set[i].dynamic2].lo_target));

    } else if (set[i].substyle == TWO_FROM_ONE) {
      set[i].lo_target = 0.5*(set[i].lo_start+set[i].hi_start) -
        0.5*sqrt(set[i].vol_start /
                 (set[set[i].dynamic1].hi_target -
                  set[set[i].dynamic1].lo_target) /
                 (set[set[i].fixed].hi_start -
                  set[set[i].fixed].lo_start) *
                 (set[i].hi_start - set[i].lo_start));
      set[i].hi_target = 0.5*(set[i].lo_start+set[i].hi_start) +
        0.5*sqrt(set[i].vol_start /
                 (set[set[i].dynamic1].hi_target -
                  set[set[i].dynamic1].lo_target) /
                 (set[set[i].fixed].hi_start -
                  set[set[i].fixed].lo_start) *
                 (set[i].hi_start - set[i].lo_start));
    }
  }

  // for triclinic, set new box shape
  // for NONE, target is current tilt
  // for TRATE, set target directly based on current time. also set h_rate
  // for WIGGLE, set target directly based on current time. also set h_rate
  // for VARIABLE, set target directly via variable eval. also set h_rate
  // for other styles, target is linear value between start and stop values

  if (triclinic) {
    double *h = domain->h;

    for (i = 3; i < 6; i++) {
      if (set[i].style == NONE) {
        if (i == 5) set[i].tilt_target = domain->xy;
        else if (i == 4) set[i].tilt_target = domain->xz;
        else if (i == 3) set[i].tilt_target = domain->yz;
      } else if (set[i].style == TRATE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        set[i].tilt_target = set[i].tilt_start * exp(set[i].rate*delt);
        h_rate[i] = set[i].rate * domain->h[i];
      } else if (set[i].style == WIGGLE) {
        double delt = (update->ntimestep - update->beginstep) * update->dt;
        set[i].tilt_target = set[i].tilt_start +
          set[i].amplitude * sin(MY_2PI*delt/set[i].tperiod);
        h_rate[i] = MY_2PI/set[i].tperiod * set[i].amplitude *
          cos(MY_2PI*delt/set[i].tperiod);
      } else if (set[i].style == VARIABLE) {
        double delta_tilt = input->variable->compute_equal(set[i].hvar);
        set[i].tilt_target = set[i].tilt_start + delta_tilt;
        h_rate[i] = input->variable->compute_equal(set[i].hratevar);
      } else {
        set[i].tilt_target = set[i].tilt_start +
          delta*(set[i].tilt_stop - set[i].tilt_start);
      }

      // tilt_target can be large positive or large negative value
      // add/subtract box lengths until tilt_target is closest to current value

      int idenom = 0;
      if (i == 5) idenom = 0;
      else if (i == 4) idenom = 0;
      else if (i == 3) idenom = 1;
      double denom = set[idenom].hi_target - set[idenom].lo_target;

      double current = h[i]/h[idenom];

      while (set[i].tilt_target/denom - current > 0.0)
        set[i].tilt_target -= denom;
      while (set[i].tilt_target/denom - current < 0.0)
        set[i].tilt_target += denom;
      if (fabs(set[i].tilt_target/denom - 1.0 - current) <
          fabs(set[i].tilt_target/denom - current))
        set[i].tilt_target -= denom;
    }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + nevery);

  // if any tilt ratios exceed 0.5, set flip = 1 and compute new tilt values
  // do not flip in x or y if non-periodic (can tilt but not flip)
  //   this is b/c the box length would be changed (dramatically) by flip
  // if yz tilt exceeded, adjust C vector by one B vector
  // if xz tilt exceeded, adjust C vector by one A vector
  // if xy tilt exceeded, adjust B vector by one A vector
  // check yz first since it may change xz, then xz check comes after
  // flip is performed on next timestep, before reneighboring in pre-exchange()

  if (triclinic && flipflag) {
    double xprd = set[0].hi_target - set[0].lo_target;
    double yprd = set[1].hi_target - set[1].lo_target;
    double xprdinv = 1.0 / xprd;
    double yprdinv = 1.0 / yprd;
    if (set[3].tilt_target*yprdinv < -0.5 ||
                                     set[3].tilt_target*yprdinv > 0.5 ||
        set[4].tilt_target*xprdinv < -0.5 ||
                                     set[4].tilt_target*xprdinv > 0.5 ||
        set[5].tilt_target*xprdinv < -0.5 ||
                                     set[5].tilt_target*xprdinv > 0.5) {
      set[3].tilt_flip = set[3].tilt_target;
      set[4].tilt_flip = set[4].tilt_target;
      set[5].tilt_flip = set[5].tilt_target;

      flipxy = flipxz = flipyz = 0;

      if (domain->yperiodic) {
        if (set[3].tilt_flip*yprdinv < -0.5) {
          set[3].tilt_flip += yprd;
          set[4].tilt_flip += set[5].tilt_flip;
          flipyz = 1;
        } else if (set[3].tilt_flip*yprdinv > 0.5) {
          set[3].tilt_flip -= yprd;
          set[4].tilt_flip -= set[5].tilt_flip;
          flipyz = -1;
        }
      }
      if (domain->xperiodic) {
        if (set[4].tilt_flip*xprdinv < -0.5) {
          set[4].tilt_flip += xprd;
          flipxz = 1;
        }
        if (set[4].tilt_flip*xprdinv > 0.5) {
          set[4].tilt_flip -= xprd;
          flipxz = -1;
        }
        if (set[5].tilt_flip*xprdinv < -0.5) {
          set[5].tilt_flip += xprd;
          flipxy = 1;
        }
        if (set[5].tilt_flip*xprdinv > 0.5) {
          set[5].tilt_flip -= xprd;
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
    int nlocal = atom->nlocal;

    domainKK->x2lamda(nlocal);

    if (rfix.size() > 0) {
      atomKK->sync(Host,ALL_MASK);
      for (auto &ifix : rfix)
        ifix->deform(0);
      atomKK->modified(Host,ALL_MASK);
    }
  }

  // reset global and local box to new size/shape
  // only if deform fix is controlling the dimension

  if (set[0].style) {
    domain->boxlo[0] = set[0].lo_target;
    domain->boxhi[0] = set[0].hi_target;
  }
  if (set[1].style) {
    domain->boxlo[1] = set[1].lo_target;
    domain->boxhi[1] = set[1].hi_target;
  }
  if (set[2].style) {
    domain->boxlo[2] = set[2].lo_target;
    domain->boxhi[2] = set[2].hi_target;
  }
  if (triclinic) {
    if (set[3].style) domain->yz = set[3].tilt_target;
    if (set[4].style) domain->xz = set[4].tilt_target;
    if (set[5].style) domain->xy = set[5].tilt_target;
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert atoms and rigid bodies back to box coords

  if (remapflag == Domain::X_REMAP) {
    int nlocal = atom->nlocal;

    domainKK->lamda2x(nlocal);

    if (rfix.size() > 0) {
      atomKK->sync(Host,ALL_MASK);
      for (auto &ifix : rfix)
        ifix->deform(1);
      atomKK->modified(Host,ALL_MASK);
    }
  }

  // redo KSpace coeffs since box has changed

  if (kspace_flag) force->kspace->setup();
}


