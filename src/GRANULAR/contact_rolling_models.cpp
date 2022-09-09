/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "contact_normal_models.h"
#include "contact_rolling_models.h"
#include "contact.h"
#include "error.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace Contact;
using namespace MathConst;
using namespace MathExtra;

/* ----------------------------------------------------------------------
   Default rolling friction model
------------------------------------------------------------------------- */

RollingModel::RollingModel(LAMMPS *lmp) : SubModel(lmp) {}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

RollingNone::RollingNone(LAMMPS *lmp) : RollingModel(lmp) {}

/* ----------------------------------------------------------------------
   SDS rolling friction model
------------------------------------------------------------------------- */

RollingSDS::RollingSDS(LAMMPS *lmp) : RollingModel(lmp)
{
  num_coeffs = 3;
  size_history = 3;
}

/* ---------------------------------------------------------------------- */

void RollingSDS::coeffs_to_local()
{
  k = coeffs[0];
  mu = coeffs[1];
  gamma = coeffs[2];

  if (k < 0.0 || mu < 0.0 || gamma < 0.0)
    error->all(FLERR, "Illegal SDS rolling model");
}

/* ---------------------------------------------------------------------- */

void RollingSDS::calculate_forces()
{
  int rhist0, rhist1, rhist2, frameupdate;
  double Frcrit, rolldotn, rollmag, prjmag, magfr, hist_temp[3], temp_dbl, temp_array[3];

  rhist0 = history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu * contact->normal_model->Fncrit;

  if (contact->history_update) {
    hist_temp[0] = contact->history[rhist0];
    hist_temp[1] = contact->history[rhist1];
    hist_temp[2] = contact->history[rhist2];
    rolldotn = dot3(hist_temp, contact->nx);

    frameupdate = fabs(rolldotn) * k  > EPSILON * Frcrit;
    if (frameupdate) { // rotate into tangential plane
      rollmag = len3(hist_temp);
      // projection
      temp_dbl = -rolldotn;
      scale3(temp_dbl, contact->nx, temp_array);
      sub3(hist_temp, temp_array, hist_temp);

      // also rescale to preserve magnitude
      prjmag = len3(hist_temp);
      if (prjmag > 0) temp_dbl = rollmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, hist_temp);
    }
    scale3(contact->dt, contact->vrl, temp_array);
    add3(hist_temp, temp_array, hist_temp);
  }

  scaleadd3(k, hist_temp, gamma, contact->vrl, contact->fr);
  negate3(contact->fr);

  // rescale frictional displacements and forces if needed

  magfr = len3(contact->fr);
  if (magfr > Frcrit) {
    rollmag = len3(hist_temp);
    if (rollmag != 0.0) {
      temp_dbl = -Frcrit / (magfr * k);
      scale3(temp_dbl, contact->fr, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = -gamma/k;
      scale3(temp_dbl, contact->vrl, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = Frcrit / magfr;
      scale3(temp_dbl, contact->fr);
    } else {
      zero3(contact->fr);
    }
  }

  contact->history[rhist0] = hist_temp[0];
  contact->history[rhist1] = hist_temp[1];
  contact->history[rhist2] = hist_temp[2];
}
