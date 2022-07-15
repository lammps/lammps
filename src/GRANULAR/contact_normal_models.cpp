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
#include "math_const.h"
#include "contact.h"

#include <cmath>

using namespace MathConst;

namespace Contact{

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

NormalModel::NormalModel()
{
  material_properties = 0;
}

/* ---------------------------------------------------------------------- */

void NormalModel::set_fncrit()
{
  Fncrit = fabs(contact.Fntot);
}

/* ---------------------------------------------------------------------- */

void NormalModel::pulloff_distance(double radi, double radj)
{
  return radi + radj;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

void NormalHooke::NormalHooke()
{
  num_coeffs = 2;
  allocate_coeffs();
}

/* ---------------------------------------------------------------------- */

void NormalHooke::coeffs_to_local()
{
  k_norm = coeffs[0];
  damp = coeffs[1];
}

/* ---------------------------------------------------------------------- */

void NormalHooke::mix_coeffs(NormalModel* imodel, NormalModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double NormalHooke::calculate_forces()
{
  contact.area = sqrt(contact.dR);
  knfac = k_norm * contact.area;
  Fne = knfac * contact.delta;
  return Fne;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

void NormalHertz::NormalHertz()
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void NormalHertz::coeffs_to_local()
{
  k_norm = coeffs[0];
  damp = coeffs[1];
}

/* ---------------------------------------------------------------------- */

void NormalHertz::mix_coeffs(NormalModel* imodel, NormalModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double NormalHertz::calculate_forces()
{
  contact.area = sqrt(contact.dR);
  knfac = contact.k_norm * contact.area;
  Fne = knfac * contact.delta;
  return Fne;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

void NormalHertzMaterial::NormalHertzMaterial()
{
  material_properties = 1;
  num_coeffs = 3;
}

/* ---------------------------------------------------------------------- */

void NormalHertzMaterial::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  k_norm = 4 / 3 * Emod;
}

/* ---------------------------------------------------------------------- */

void NormalHertzMaterial::mix_coeffs(NormalModel* imodel, NormalModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs[2] = mix_geom(imodel->coeffs[2], jmodel->coeffs[2]);
  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

void NormalDMT::NormalDMT(ContactModel &c)
{
  allow_limit_damping = 0;
  material_properties = 1;
  num_coeffs = 4;
}

/* ---------------------------------------------------------------------- */

void NormalDMT::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];
  k_norm = 4 / 3 * Emod;
}

/* ---------------------------------------------------------------------- */

void NormalDMT::mix_coeffs(NormalModel* imodel, NormalModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs[2] = mix_geom(imodel->coeffs[2], jmodel->coeffs[2]);
  coeffs[3] = mix_geom(imodel->coeffs[3], jmodel->coeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double NormalDMT::calculate_forces()
{
  contact.area = sqrt(contact.dR);
  knfac = k_norm * contact.area;
  Fne = knfac * contact.delta;
  F_pulloff = 4 * MathConst::MY_PI * cohesion * contact.Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void NormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2* F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

void NormalJKR::NormalJKR(ContactModel &c)
{
  allow_limit_damping = 0;
  material_properties = 1;
  contact.beyond_contact = beyond_contact = 1;
  num_coeffs = 4;
}

/* ---------------------------------------------------------------------- */

void NormalJKR::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];
  k_norm = 4/3*Emod;
  Escaled = k_norm * THREEQUARTERS;
}

/* ---------------------------------------------------------------------- */

void NormalJKR::mix_coeffs(NormalModel* imodel, NormalModel* jmodel)
{
  coeffs[0] = mix_geom(imodel->coeffs[0], jmodel->coeffs[0]);
  coeffs[1] = mix_geom(imodel->coeffs[1], jmodel->coeffs[1]);
  coeffs[2] = mix_geom(imodel->coeffs[2], jmodel->coeffs[2]);
  coeffs[3] = mix_geom(imodel->coeffs[3], jmodel->coeffs[3]);
  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool NormalJKR::touch(int touch)
{
  double R2, delta_pulloff, dist_pulloff;
  bool touchflag;

  if (touch) {
    R2 = contact.Reff * contact.Reff;
    a = cbrt(9.0 * MY_PI * cohesion * R2 / (4 * Escaled));
    delta_pulloff = a * a / contact.Reff - 2 * sqrt(MY_PI * cohesion * a / Escaled);
    dist_pulloff = contact.radsum - delta_pulloff;
    touchflag = (contact.rsq < dist_pulloff * dist_pulloff);
  } else {
    touchflag = (rcontact.sq < contact.radsum * contact.radsum);
  }
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double NormalJKR::calculate_forces()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3, a2;

  R2 = Reff * Reff;
  dR2 = dR * dR;
  t0 = cohesion * cohesion * R2 * R2 * Escaled;
  t1 = PI27SQ*t0;
  t2 = 8 * dR * dR2 * Escaled * Escaled * Escaled;
  t3 = 4 * dR2 * Escaled;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / Escaled;
  sqrt2 = MAX(0, 2 * dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4 * dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (Escaled * t6));
  contact.area = INVROOT6 * (t6 + sqrt(sqrt3));
  a2 = contact.area * contact.area;
  Fne = Escaled * contact.area * a2 / Reff - MY_2PI * a2 * sqrt(4 * cohesion * Escaled / (MY_PI * contact.area));
  F_pulloff = 3 * MY_PI * cohesion * Reff;

  knfac = Escaled * contact.area;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void NormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2 * F_pulloff);
}

/* ---------------------------------------------------------------------- */

void NormalJKR::pulloff_distance(double radi, double radj)
{
  double a_tmp, Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj);
  if (Reff_tmp <= 0) return 0;

  a_tmp = cbrt(9 * MY_PI * cohesion * Reff_tmp * Reff_tmp / (4 * Ecaled));
  return a_tmp * a_tmp / Reff_tmp - 2 * sqrt(MY_PI * cohesion * a_tmp / Ecaled);
}

