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

#include "damping_contact_models.h"
#include "math_const.h"
#include "contact.h"

#include <cmath>
#include "math_special.h"

using namespace MathConst;
using namespace MathSpecial;

namespace Contact{

void DampingModel::allocate_coeffs(){
  damp = contact.normal_model->coeffs[1];
}

//-----------------------------------------
double Velocity::calculate_forces(){
  return -damp*contact.vnnr;
}

double MassVelocity::calculate_forces(){
  return -damp*contact.meff*contact.vnnr;
}


double ViscoElastic::calculate_forces(){
  return -damp*contact.meff*contact.area*contact.vnnr;
}


void Tsuji::allocate_coeffs(){
  double cor = contact.normal_model->coeffs[1];
  damp = 1.2728-4.2783*cor+11.087*square(cor)-22.348*cube(cor)+
      27.467*powint(cor,4)-18.022*powint(cor,5)+4.8218*powint(cor,6);
}

double Tsuji::calculate_forces(){
  return -damp*sqrt(contact.meff*contact.knfac)*contact.vnnr;
}

