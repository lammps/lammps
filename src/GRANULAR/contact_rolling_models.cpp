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

#include "normal_contact_models.h"
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

// ************************
// Default behaviors where needed
// ************************
void NormalModel::set_fncrit(){
  contact->Fncrit = fabs(contact->Fntot);
}

void NormalModel::mix_coeffs(NormalModel* imodel, NormalModel* jmodel){
  for (int i = 0; i < num_coeffs; i++){
    coeffs[i] = sqrt(imodel->coeffs[i]*jmodel->coeffs[i]);
  }
}

//-----------------------------------------

//******************
// Hooke
//******************
void Hooke::Hooke(ContactModel &c){
  contact = c;
  num_coeffs = 2;
  allocate_coeffs();
}

void Hooke::coeffs_to_local(){
  k_norm = coeffs[0];
  damp = coeffs[1];
}

double Hooke::calculate_forces(){
  contact.area = sqrt(contact.dR);
  contact.knfac = k_norm * contact.area;
  return contact.knfac * contact.delta;
}


//******************
// Hertz
//******************
void Hertz::Hertz(ContactModel &c, int mat_flag){
  contact = c;
  material_prop_flag = mat_flag;
  if (material_prop_flag){
    num_coeffs = 3;
  }
  else{
    num_coeffs = 2;
  }
  allocate_coeffs();
}

void Hertz::coeffs_to_local(){
  if (material_prop_flag){
    Emod = coeffs[0];
    poiss = coeffs[1];
    k_norm = 4/3*Emod;
  }
  else{
    k_norm = coeffs[0];
  }
}

double Hertz::calculate_forces(){
  contact.area = sqrt(contact.dR);
  contact.knfac = contact.k_norm * contact.area;
  return contact.knfac * contact.delta;
}

//******************
// DMT
//******************
void DMT::DMT(ContactModel &c){
  contact = c;
  material_prop_flag = 1;
  num_coeffs = 4;
  allocate_coeffs();
}

double DMT::calculate_forces(){
  contact.area = sqrt(contact.dR);
  contact.knfac = contact.k_norm * contact.area;
  double Fne = contact.knfac * contact.delta;
  F_pulloff = 4 * MathConst::MY_PI * cohesion * contact.Reff;
  Fne -= F_pulloff;
  return Fne;
}

void DMT::set_fncrit(){
  contact.Fncrit = fabs(contact.Fne + 2* F_pulloff);
}

//******************
// JKR
//******************
void JKR::JKR(ContactModel &c){
  contact = c;
  material_prop_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
  allocate_coeffs();
}

bool JKR::touch(int touch){
  double Escaled, R2, delta_pulloff, dist_pulloff;
  bool touchflag;

  Escaled = k_norm * THREEQUARTERS;
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

double JKR::calculate_forces(){
  double Escaled, R2, dR2, t0, t1, t2, t3, t4, t5, t6;
    double sqrt1, sqrt2, sqrt3, a2, F_pulloff, Fne;

    Escaled = k_norm * THREEQUARTERS;

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
    a = INVROOT6 * (t6 + sqrt(sqrt3));
    a2 = a * a;
    Fne = Escaled * a * a2 / Reff - MY_2PI * a2 * sqrt(4 * cohesion * Escaled / (MY_PI * a));
    F_pulloff = 3 * MY_PI * cohesion * Reff;

    knfac = Escaled * a;
    return Fne;
}

void JKR::set_fncrit(){
  contact.Fncrit = fabs(contact.Fne + 2 * F_pulloff);
}




void ContactModel::rolling(double *history)
{
  int rhist0, rhist1, rhist2, frameupdate;
  double rolldotn, rollmag, prjmag, magfr, hist_temp[3], temp_dbl, temp_array[3];

  rhist0 = roll_history_index;
  rhist1 = rhist0 + 1;
  rhist2 = rhist1 + 1;

  Frcrit = mu_roll * Fncrit;

  if (history_update) {
    hist_temp[0] = history[rhist0];
    hist_temp[1] = history[rhist1];
    hist_temp[2] = history[rhist2];
    rolldotn = dot3(hist_temp, nx);

    frameupdate = fabs(rolldotn)*k_roll > EPSILON*Frcrit;
    if (frameupdate) { // rotate into tangential plane
      rollmag = len3(hist_temp);
      // projection
      temp_dbl = -rolldotn;
      scale3(temp_dbl, nx, temp_array);
      sub3(hist_temp, temp_array, hist_temp);

      // also rescale to preserve magnitude
      prjmag = len3(hist_temp);
      if (prjmag > 0) temp_dbl = rollmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, hist_temp);
    }
    scale3(dt, vrl, temp_array);
    add3(hist_temp, temp_array, hist_temp);
  }

  scaleadd3(k_roll, hist_temp, gamma_roll, vrl, fr);
  negate3(fr);

  // rescale frictional displacements and forces if needed

  magfr = len3(fr);
  if (magfr > Frcrit) {
    rollmag = len3(hist_temp);
    if (rollmag != 0.0) {
      temp_dbl = -Frcrit / (magfr * k_roll);
      scale3(temp_dbl, fr, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = -gamma_roll/k_roll;
      scale3(temp_dbl, vrl, temp_array);
      add3(hist_temp, temp_array, hist_temp);

      temp_dbl = Frcrit / magfr;
      scale3(temp_dbl, fr);
    } else {
      zero3(fr);
    }
  }

  history[rhist0] = hist_temp[0];
  history[rhist1] = hist_temp[1];
  history[rhist2] = hist_temp[2];
}
