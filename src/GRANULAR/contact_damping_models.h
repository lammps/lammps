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

class DampingModel:SubModel
{
 public:
  DampingModel();
  virtual ~DampingModel() {};
  virtual double calculate_forces() = 0;
  virtual void coeffs_to_local();
  virtual void mix_coeffs(NormalModel*, NormalModel*);
  double damp;
};

/* ---------------------------------------------------------------------- */

class DampingVelocity:DampingModel
{
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingMassVelocity:DampingModel
{
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingViscoElastic:DampingModel
{
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingTsuji:DampingModel
{
 public:
  double calculate_forces();
  void coeffs_to_local();
};


}

#endif /*DAMPING_CONTACT_MODELS_H_ */

