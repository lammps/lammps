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

#include "gran_sub_mod_tangential.h"

#include "error.h"
#include "gran_sub_mod_damping.h"
#include "gran_sub_mod_normal.h"
#include "granular_model.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace Granular_NS::GranKernel;

/* ----------------------------------------------------------------------
   Default model
------------------------------------------------------------------------- */

GranSubModTangential::GranSubModTangential(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  allow_synchronization = 0;
  mindlin_force = 0;
  mindlin_rescale = 0;
}

/* ----------------------------------------------------------------------
   collect the per-contact inputs shared by all tangential kernels
------------------------------------------------------------------------- */

static GranKernel::GranTangentialState<double> tangential_state(GranularModel *gm,
                                                                double damp_prefactor)
{
  GranKernel::GranTangentialState<double> s;
  s.nx = gm->nx;
  s.nx_unrotated = gm->nx_unrotated;
  s.vtr = gm->vtr;
  s.vrel = gm->vrel;
  s.dt = gm->dt;
  s.contact_radius = gm->contact_radius;
  s.damp_prefactor = damp_prefactor;
  s.Fncrit = gm->normal_model->get_fncrit();
  s.synchronized_verlet = gm->synchronized_verlet;
  s.history_update = gm->history_update;
  return s;
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModTangentialNone::GranSubModTangentialNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  allow_synchronization = 1;
}

/* ----------------------------------------------------------------------
   Linear model with no history
------------------------------------------------------------------------- */

GranSubModTangentialLinearNoHistory::GranSubModTangentialLinearNoHistory(GranularModel *gm,
                                                                         LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 2;
  size_history = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::coeffs_to_local()
{
  k = 0.0;    // No tangential stiffness with no history
  xt = coeffs[0];
  mu = coeffs[1];

  if (k < 0.0 || xt < 0.0 || mu < 0.0)
    error->all(FLERR, "Illegal linear no history tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearNoHistory::calculate_forces()
{
  // classic pair gran/hooke (no history)
  GranTangentialParams<double> p;
  fill_kernel_params(p);
  const GranTangentialState<double> s = tangential_state(gm, gm->damping_model->get_damp_prefactor());

  damp = gran_tangential_damp(p, s.damp_prefactor);
  gran_tangential_linear_nohistory(p, s, gm->fs);
}

/* ----------------------------------------------------------------------
   Linear model with history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistory::GranSubModTangentialLinearHistory(GranularModel *gm,
                                                                     LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal linear tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistory::calculate_forces()
{
  // Note: this is the same as the base Mindlin calculation except k isn't scaled by contact radius
  GranTangentialParams<double> p;
  fill_kernel_params(p);
  const GranTangentialState<double> s = tangential_state(gm, gm->damping_model->get_damp_prefactor());

  damp = gran_tangential_damp(p, s.damp_prefactor);
  gran_tangential_linear_history(p, s, &gm->history[history_index], gm->fs);
}

/* ----------------------------------------------------------------------
   Linear model with history from pair gran/hooke/history
------------------------------------------------------------------------- */

GranSubModTangentialLinearHistoryClassic::GranSubModTangentialLinearHistoryClassic(
    GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialLinearHistory(gm, lmp)
{
  // Whether contact radius scales normal force (as in Hertz)
  contact_radius_flag = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialLinearHistoryClassic::calculate_forces()
{
  GranTangentialParams<double> p;
  fill_kernel_params(p);
  const GranTangentialState<double> s = tangential_state(gm, gm->damping_model->get_damp_prefactor());

  damp = gran_tangential_damp(p, s.damp_prefactor);
  gran_tangential_classic(p, s, &gm->history[history_index], gm->fs);
}

/* ----------------------------------------------------------------------
   Mindlin from pair gran/hertz/history
------------------------------------------------------------------------- */

GranSubModTangentialMindlinClassic::GranSubModTangentialMindlinClassic(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialLinearHistoryClassic(gm, lmp)
{
  contact_radius_flag = 1;    // Sets gran/hertz/history behavior
}

/* ----------------------------------------------------------------------
   Mindlin model
------------------------------------------------------------------------- */

GranSubModTangentialMindlin::GranSubModTangentialMindlin(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangential(gm, lmp)
{
  num_coeffs = 3;
  size_history = 3;
  contact_radius_flag = 1;
  allow_synchronization = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::coeffs_to_local()
{
  k = coeffs[0];
  xt = coeffs[1];
  mu = coeffs[2];

  if (k == -1) {
    if (!gm->normal_model->get_material_properties())
      error->all(FLERR,
                 "Must either specify tangential stiffness or material properties for normal model "
                 "for the Mindlin tangential style");

    double Emod = gm->normal_model->get_emod();
    double poiss = gm->normal_model->get_poiss();

    if (gm->contact_type == PAIR) {
      k = 8.0 * mix_stiffnessG(Emod, Emod, poiss, poiss);
    } else {
      k = 8.0 * mix_stiffnessG_wall(Emod, poiss);
    }
  }

  if (k < 0.0 || xt < 0.0 || mu < 0.0) error->all(FLERR, "Illegal Mindlin tangential model");
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  if (icoeffs[0] == -1 || jcoeffs[0] == -1)
    coeffs[0] = -1;
  else
    coeffs[0] = mix_geom(icoeffs[0], jcoeffs[0]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

void GranSubModTangentialMindlin::calculate_forces()
{
  GranTangentialParams<double> p;
  fill_kernel_params(p);
  const GranTangentialState<double> s = tangential_state(gm, gm->damping_model->get_damp_prefactor());

  damp = gran_tangential_damp(p, s.damp_prefactor);
  gran_tangential_mindlin(p, s, &gm->history[history_index], gm->fs);
}

/* ----------------------------------------------------------------------
   Mindlin force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinForce::GranSubModTangentialMindlinForce(GranularModel *gm, LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  mindlin_force = 1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescale::GranSubModTangentialMindlinRescale(GranularModel *gm,
                                                                       LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}

/* ----------------------------------------------------------------------
   Mindlin rescale force model
------------------------------------------------------------------------- */

GranSubModTangentialMindlinRescaleForce::GranSubModTangentialMindlinRescaleForce(GranularModel *gm,
                                                                                 LAMMPS *lmp) :
    GranSubModTangentialMindlin(gm, lmp)
{
  size_history = 4;
  mindlin_force = 1;
  mindlin_rescale = 1;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) transfer_history_factor[i] = -1.0;
  transfer_history_factor[3] = +1;
}
