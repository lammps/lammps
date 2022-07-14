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
  material_prop_flag = mat_flag;
  if (material_prop_flag){
    num_coeffs = 3;
  }
  else{
    num_coeffs = 2;
  }
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
  material_prop_flag = 1;
  num_coeffs = 4;
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
  material_prop_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
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



