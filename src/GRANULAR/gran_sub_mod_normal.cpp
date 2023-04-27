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

#include "gran_sub_mod_normal.h"
#include "error.h"
#include "granular_model.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;

using MathConst::MY_2PI;
using MathConst::MY_PI;

static constexpr double PI27SQ = 266.47931882941264802866;      // 27*PI**2
static constexpr double THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double FOURTHIRDS = (4.0 / 3.0);               // 4/3
static constexpr double JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI**2/16)

/* ----------------------------------------------------------------------
   Default normal model
------------------------------------------------------------------------- */

GranSubModNormal::GranSubModNormal(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  material_properties = 0;
  cohesive_flag = 0;
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormal::touch()
{
  bool touchflag = (gm->rsq < gm->radsum * gm->radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::pulloff_distance(double /*radi*/, double /*radj*/)
{
  // called outside of compute(), do not assume correct geometry defined in contact
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::calculate_contact_radius()
{
  return sqrt(gm->dR);
}

/* ---------------------------------------------------------------------- */

void GranSubModNormal::set_fncrit()
{
  Fncrit = fabs(gm->Fntot);
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModNormalNone::GranSubModNormalNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalNone::calculate_forces()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

GranSubModNormalHooke::GranSubModNormalHooke(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHooke::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hooke normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHooke::calculate_forces()
{
  return k * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

GranSubModNormalHertz::GranSubModNormalHertz(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertz::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHertz::calculate_forces()
{
  return k * gm->contact_radius * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

GranSubModNormalHertzMaterial::GranSubModNormalHertzMaterial(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormalHertz(gm, lmp)
{
  material_properties = 1;
  num_coeffs = 3;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  if (gm->contact_type == PAIR) {
    k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
  } else {
    k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz material normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

GranSubModNormalDMT::GranSubModNormalDMT(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];
  if (gm->contact_type == PAIR) {
    k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
  } else {
    k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal DMT normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalDMT::calculate_forces()
{
  Fne = k * gm->contact_radius * gm->delta;
  F_pulloff = 4.0 * MY_PI * cohesion * gm->Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

GranSubModNormalJKR::GranSubModNormalJKR(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (gm->contact_type == PAIR) {
    Emix = mix_stiffnessE(Emod, Emod, poiss, poiss);
  } else {
    Emix = mix_stiffnessE_wall(Emod, poiss);
  }

  k = FOURTHIRDS * Emix;

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal JKR normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormalJKR::touch()
{
  double delta_pulloff, dist_pulloff;
  bool touchflag;

  if (gm->touch) {
    // delta_pulloff defined as positive so center-to-center separation is > radsum
    delta_pulloff = JKRPREFIX * cbrt(gm->Reff * cohesion * cohesion / (Emix * Emix));
    dist_pulloff = gm->radsum + delta_pulloff;
    touchflag = gm->rsq < (dist_pulloff * dist_pulloff);
  } else {
    touchflag = gm->rsq < (gm->radsum * gm->radsum);
  }

  return touchflag;
}

/* ----------------------------------------------------------------------
  called outside of compute(), do not assume geometry defined in contact
------------------------------------------------------------------------- */

double GranSubModNormalJKR::pulloff_distance(double radi, double radj)
{
  double Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj);    // May not be defined
  if (Reff_tmp <= 0) return 0;
  // Defined as positive so center-to-center separation is > radsum
  return JKRPREFIX * cbrt(Reff_tmp * cohesion * cohesion / (Emix * Emix));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_contact_radius()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  R2 = gm->Reff * gm->Reff;
  dR2 = gm->dR * gm->dR;
  t0 = cohesion * cohesion * R2 * R2 * Emix;
  t1 = PI27SQ * t0;
  t2 = 8.0 * gm->dR * dR2 * Emix * Emix * Emix;
  t3 = 4.0 * dR2 * Emix;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2.0 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / Emix;
  sqrt2 = MAX(0, 2.0 * gm->dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4.0 * gm->dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (Emix * t6));

  return INVROOT6 * (t6 + sqrt(sqrt3));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_forces()
{
  double a2;
  a2 = gm->contact_radius * gm->contact_radius;
  Fne = k * gm->contact_radius * a2 / gm->Reff -
      MY_2PI * a2 * sqrt(4.0 * cohesion * Emix / (MY_PI * gm->contact_radius));
  F_pulloff = 3.0 * MY_PI * cohesion * gm->Reff;

  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}
