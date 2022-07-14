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

#ifndef DAMPING_CONTACT_MODELS_H_
#define DAMPING_CONTACT_MODELS_H_

#include "contact.h";
#include "sub_model.h"

namespace Contact {

class DampingModel:SubModel{
public:
  DampingModel(){};
  virtual ~DampingModel(){};
  virtual double calculate_forces() = 0;
  virtual void allocate_coeffs();
  double damp;
};

class Velocity:DampingModel{
public:
  double calculate_forces();
};

class MassVelocity:DampingModel{
public:
  double calculate_forces();
};

class ViscoElastic:DampingModel{
public:
  double calculate_forces();
};

class Tsuji:DampingModel{
public:
  double calculate_forces();
  void allocate_coeffs();
};


}

#endif /*DAMPING_CONTACT_MODELS_H_ */

