/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "gran_sub_mod_rolling.h"

#include "gran_sub_mod_rolling_kernel.h"

#include "error.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace Granular_NS::GranKernel;

/* ----------------------------------------------------------------------
   Default rolling friction model
------------------------------------------------------------------------- */

GranSubModRolling::GranSubModRolling(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  allow_synchronization = 0;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModRollingNone::GranSubModRollingNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModRolling(gm, lmp)
{
  allow_synchronization = 1;
}

/* ----------------------------------------------------------------------
   SDS rolling friction model
------------------------------------------------------------------------- */

GranSubModRollingSDS::GranSubModRollingSDS(GranularModel *gm, LAMMPS *lmp) :
    GranSubModRolling(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  allow_synchronization = 1;

  // the rolling displacement is the same seen from either particle:
  // vrl = Reff * (relrot x nx) is unchanged when both relrot and nx flip
  // sign, so this history must NOT be negated when it is transferred
  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = +1.0;
}

/* ---------------------------------------------------------------------- */

void GranSubModRollingSDS::coeffs_to_local()
{
  k = coeffs[0];
  gamma = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || mu < 0.0 || gamma < 0.0) error->all(FLERR, "Illegal SDS rolling model");
}

/* ---------------------------------------------------------------------- */

void GranSubModRollingSDS::calculate_forces()
{
  GranRollingParams<double> p;
  p.model = GRAN_ROLLING_NONE;    // host classes call the models directly
  p.k = k;
  p.gamma = gamma;
  p.mu = mu;

  GranRollingState<double> s;
  s.nx = gm->nx;
  s.nx_unrotated = gm->nx_unrotated;
  s.vrl = gm->vrl;
  s.dt = gm->dt;
  s.Fncrit = gm->normal_model->get_fncrit();
  s.synchronized_verlet = gm->synchronized_verlet;
  s.history_update = gm->history_update;

  gran_rolling_sds(p, s, &gm->history[history_index], gm->fr);
}
