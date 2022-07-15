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

#ifndef NORMAL_CONTACT_MODELS_H_
#define NORMAL_CONTACT_MODELS_H_

#include "contact.h";
#include "sub_model.h"

namespace Contact {

class NormalModel:SubModel{
public:
  NormalModel(){};
  virtual ~NormalModel(){};
  virtual bool check_contact();
  virtual void prep_contact();
  virtual void set_fncrit();
  virtual double calculate_forces() = 0;
  virtual void coeffs_to_local();
  virtual void mix_coeffs(NormalModel*, NormalModel*); //When mixing is needed

private:
  int beyond_contact = 0;
  int allow_limit_damping = 1;

};

class Hooke:NormalModel{
public:
  Hooke(ContactModel &c);
  ~Hooke(){};
  void coeffs_to_local();
  double calculate_forces();
private:
  double k_norm, damp;
};

class Hertz:NormalModel{
public:
  Hertz(ContactModel&, int);
  ~Hertz(){};
  void coeffs_to_local();
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss;
};

class DMT:NormalModel{
public:
  DMT(ContactModel &c);
  ~DMT(){};
  void coeffs_to_local();
  void coeffs_to_local(NormalModel*, NormalModel*);
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss, coh;
};

class JKR:NormalModel{
public:
  JKR(ContactModel &c);
  ~JKR(){};
  void coeffs_to_local();
  double calculate_forces();
private:
  double k_norm, damp, Emod, poiss, coh;

};
}

#endif /*NORMAL_CONTACT_MODELS_H_ */

