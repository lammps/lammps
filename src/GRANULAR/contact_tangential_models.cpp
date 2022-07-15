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

#include "tangential_contact_models.h"
#include "math_const.h"
#include "contact.h"

#include <cmath>

using namespace MathConst;

namespace Contact{

// ************************
// Default behaviors where needed
// ************************

//-----------------------------------------

//******************
// Hooke
//******************
void LinearNohistory::calculate_forces(){
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





/* ---------------------------------------------------------------------- */

void ContactModel::tangential_no_history()
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double fsmag, Ft;

  // classic pair gran/hooke (no history)
  fsmag = gamma_scaled * vrel;
  if (vrel != 0.0) Ft = MIN(Fscrit,fsmag) / vrel;
  else Ft = 0.0;

  Ft = -Ft;
  scale3(Ft, vtr, fs);
}

/* ---------------------------------------------------------------------- */

void ContactModel::tangential_history(double *history)
{
  double gamma_scaled = gamma_tang * damp_normal_prefactor;
  double k = k_tang;
  int frame_update = 0;
  double magfs, rsht, shrmag, prjmag, temp_dbl, temp_array[3];

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    frame_update = fabs(rsht) * k > EPSILON * Fscrit;

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    // tangential force
    // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
    temp_dbl = k * dt;
    scale3(temp_dbl, vtr, temp_array);
    sub3(history, temp_array, history);
  }

  // tangential forces = history + tangential velocity damping
  temp_dbl = -gamma_norm;
  scale3(temp_dbl, vtr, fs);

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs, history);
      scale3(gamma_norm, vtr, temp_array);
      add3(history, temp_array, history);
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs);
    } else {
      zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ContactModel::tangential_mindlin(double *history)
{
  double k_scaled, gamma_scaled, magfs, rsht, shrmag, prjmag, temp_dbl;
  double temp_array[3];
  int frame_update = 0;

  gamma_scaled = gamma_tang * damp_normal_prefactor;
  k_scaled = k_tang * a;
  if (mindlin_rescale) {
    // on unloading, rescale the shear displacements/force
    if (a < history[3]) {
      temp_dbl = a / history[3];
      scale3(temp_dbl, history);
    }
  }

  // rotate and update displacements / force.
  // see e.g. eq. 17 of Luding, Gran. Matter 2008, v10,p235
  if (history_update) {
    rsht = dot3(history, nx);
    if (mindlin_force)
      frame_update = fabs(rsht) > EPSILON * Fscrit;
    else
      frame_update = fabs(rsht) * k_scaled > EPSILON * Fscrit;

    if (frame_update) {
      shrmag = len3(history);
      // projection
      scale3(rsht, nx, history);
      // also rescale to preserve magnitude
      prjmag = len3(history);
      if (prjmag > 0) temp_dbl = shrmag / prjmag;
      else temp_dbl = 0;
      scale3(temp_dbl, history);
    }

    // update history
    if (mindlin_force) {
      // tangential force
      // see e.g. eq. 18 of Thornton et al, Pow. Tech. 2013, v223,p30-46
      temp_dbl = -k_scaled * dt;
      scale3(temp_dbl, vtr, temp_array);
    } else {
      scale3(dt, vtr, temp_array);
    }
    add3(history, temp_array, history);

    if (mindlin_rescale) history[3] = a;
  }

  // tangential forces = history + tangential velocity damping
  temp_dbl = -gamma_scaled;
  scale3(temp_dbl, vtr, fs);

  if (! mindlin_force) {
    scale3(k_scaled, history, temp_array);
    add3(fs, temp_array, fs);
  }

  // rescale frictional displacements and forces if needed
  magfs = len3(fs);
  if (magfs > Fscrit) {
    shrmag = len3(history);
    if (shrmag != 0.0) {
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs, history);
      scale3(gamma_tang, vtr, temp_array);
      add3(history, temp_array, history);
      if (! mindlin_force) {
        temp_dbl = -1.0 / k_tang;
        scale3(temp_dbl, history);
      }
      temp_dbl = Fscrit / magfs;
      scale3(temp_dbl, fs);
    } else {
      zero3(fs);
    }
  }
}

/* ---------------------------------------------------------------------- */
