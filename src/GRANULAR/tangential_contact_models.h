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

#ifndef TANGENTIAL_CONTACT_MODELS_H_
#define TANGENTIAL_CONTACT_MODELS_H_

#include "contact.h";
#include "sub_model.h"

namespace Contact {

class TangentialModel:SubModel{
public:
  TangentialModel(){};
  virtual ~TangentialModel(){};
  virtual double calculate_forces() = 0;
  virtual void coeffs_to_local();
  virtual void mix_coeffs(TangentialModel*, TangentialModel*); //When mixing is needed
  int rescale_flag = 0;
private:
  int beyond_contact = 0;
  int allow_limit_damping = 1;

};

class LinearNohistory:TangentialModel{
public:
  double calculate_forces();
private:
  void allocate_coeffs();
  double k_t, damp;
};

class LinearHistory:TangentialModel{
public:
  LinearHistory(){}
  ~LinearHistory(){};
  void coeffs_to_local();
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss;
};

class Mindlin:TangentialModel{
public:
  Mindlin(int);
  ~Mindlin(){};
  void coeffs_to_local();
  void coeffs_to_local(TangentialModel*, TangentialModel*);
  double calculate_forces();
private:

  double k_norm, damp, Emod, poiss, coh;
};

class MindlinForce:TangentialModel{
public:
  JKR(ContactModel &c);
  ~JKR(){};
  void coeffs_to_local();
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss, coh;

};
}

#endif /*TANGENTIAL_CONTACT_MODELS_H_ */

