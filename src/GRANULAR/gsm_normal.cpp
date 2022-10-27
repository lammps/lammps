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

#include "gsm_normal.h"
#include "granular_model.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathConst;

#define PI27SQ 266.47931882941264802866    // 27*PI**2
#define THREEROOT3 5.19615242270663202362  // 3*sqrt(3)
#define SIXROOT6 14.69693845669906728801   // 6*sqrt(6)
#define INVROOT6 0.40824829046386307274    // 1/sqrt(6)
#define FOURTHIRDS (4.0/3.0)               // 4/3
#define ONETHIRD (1.0/3.0)                 // 1/3
#define THREEQUARTERS 0.75                 // 3/4

/* ----------------------------------------------------------------------
   Default normal model
------------------------------------------------------------------------- */

GSMNormal::GSMNormal(GranularModel *gm, LAMMPS *lmp) : GSM(gm, lmp)
{
  material_properties = 0;
}

/* ---------------------------------------------------------------------- */

bool GSMNormal::touch()
{
  bool touchflag = (gm->rsq < gm->radsum * gm->radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double GSMNormal::pulloff_distance(double radi, double radj)
{
  //called outside of compute(), do not assume correct geometry defined in contact
  return radi + radj;
}

/* ---------------------------------------------------------------------- */

double GSMNormal::calculate_area()
{
  return sqrt(gm->dR);
}

/* ---------------------------------------------------------------------- */

void GSMNormal::set_fncrit()
{
  Fncrit = fabs(gm->Fntot);
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GSMNormalNone::GSMNormalNone(GranularModel *gm, LAMMPS *lmp) : GSMNormal(gm, lmp) {}

/* ---------------------------------------------------------------------- */

double GSMNormalNone::calculate_forces()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

GSMNormalHooke::GSMNormalHooke(GranularModel *gm, LAMMPS *lmp) : GSMNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GSMNormalHooke::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hooke normal model");
}

/* ---------------------------------------------------------------------- */

double GSMNormalHooke::calculate_forces()
{
  Fne = knfac * gm->delta;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GSMNormalHooke::set_knfac()
{
  knfac = k;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

GSMNormalHertz::GSMNormalHertz(GranularModel *gm, LAMMPS *lmp) : GSMNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GSMNormalHertz::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz normal model");
}

/* ---------------------------------------------------------------------- */

double GSMNormalHertz::calculate_forces()
{
  Fne = knfac * gm->delta;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GSMNormalHertz::set_knfac()
{
  knfac = k * gm->area;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

GSMNormalHertzMaterial::GSMNormalHertzMaterial(GranularModel *gm, LAMMPS *lmp) : GSMNormalHertz(gm, lmp)
{
  material_properties = 1;
  num_coeffs = 3;
}

/* ---------------------------------------------------------------------- */

void GSMNormalHertzMaterial::coeffs_to_local()
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

void GSMNormalHertzMaterial::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0],icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

GSMNormalDMT::GSMNormalDMT(GranularModel *gm, LAMMPS *lmp) : GSMNormal(gm, lmp)
{
  allow_limit_damping = 0;
  material_properties = 1;
  num_coeffs = 4;
}

/* ---------------------------------------------------------------------- */

void GSMNormalDMT::coeffs_to_local()
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

void GSMNormalDMT::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0],icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double GSMNormalDMT::calculate_forces()
{
  Fne = knfac * gm->delta;
  F_pulloff = 4.0 * MathConst::MY_PI * cohesion * gm->Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GSMNormalDMT::set_knfac()
{
  knfac = k * gm->area;
}

/* ---------------------------------------------------------------------- */

void GSMNormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

GSMNormalJKR::GSMNormalJKR(GranularModel *gm, LAMMPS *lmp) : GSMNormal(gm, lmp)
{
  allow_limit_damping = 0;
  material_properties = 1;
  beyond_contact = 1;
  num_coeffs = 4;
}

/* ---------------------------------------------------------------------- */

void GSMNormalJKR::coeffs_to_local()
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

void GSMNormalJKR::mix_coeffs(double* icoeffs, double* jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0],icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool GSMNormalJKR::touch()
{
  double area_at_pulloff, R2, delta_pulloff, dist_pulloff;
  bool touchflag;

  if (gm->touch) {
    R2 = gm->Reff * gm->Reff;
    area_at_pulloff = cbrt(9.0 * MY_PI * cohesion * R2 / (4.0 * Emix));
    delta_pulloff = area_at_pulloff * area_at_pulloff / gm->Reff - 2.0 * sqrt(MY_PI * cohesion * area_at_pulloff / Emix);
    dist_pulloff = gm->radsum - delta_pulloff;
    touchflag = gm->rsq < (dist_pulloff * dist_pulloff);
  } else {
    touchflag = gm->rsq < (gm->radsum * gm->radsum);
  }

  return touchflag;
}

/* ----------------------------------------------------------------------
  called outside of compute(), do not assume geometry defined in contact
------------------------------------------------------------------------- */

double GSMNormalJKR::pulloff_distance(double radi, double radj)
{
  double area_at_pulloff, Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj); // May not be defined
  if (Reff_tmp <= 0) return 0;

  area_at_pulloff = cbrt(9.0 * MY_PI * cohesion * Reff_tmp * Reff_tmp / (4.0 * Emix));
  return area_at_pulloff * area_at_pulloff / Reff_tmp - 2.0 * sqrt(MY_PI * cohesion * area_at_pulloff / Emix);
}

/* ---------------------------------------------------------------------- */

double GSMNormalJKR::calculate_area()
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

double GSMNormalJKR::calculate_forces()
{
  double a2;
  a2 = gm->area * gm->area;
  Fne = knfac * a2 / gm->Reff - MY_2PI * a2 * sqrt(4.0 * cohesion * Emix / (MY_PI * gm->area));
  F_pulloff = 3.0 * MY_PI * cohesion * gm->Reff;

  return Fne;
}

/* ---------------------------------------------------------------------- */

void GSMNormalJKR::set_knfac()
{
  knfac = k * gm->area;
}

/* ---------------------------------------------------------------------- */

void GSMNormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}
