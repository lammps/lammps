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
  NormalModel();
  virtual ~NormalModel() {};
  virtual bool check_contact();
  virtual void set_fncrit();
  virtual double calculate_forces() = 0;
  virtual void coeffs_to_local();
  virtual void mix_coeffs(NormalModel*, NormalModel*);
  void pulloff_distance(double, double);
  double damp;  // Vestigial argument needed by damping
  double Fncrit, Fne, knfac;
  int material_properties;
};

/* ---------------------------------------------------------------------- */

class NormalHooke:NormalModel
{
 public:
  NormalHooke();
  ~NormalHooke() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  double calculate_forces();
  double Emod, poiss;
};

/* ---------------------------------------------------------------------- */

class NormalHertz:NormalModel
{
 public:
  NormalHertz();
  ~NormalHertz() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  double calculate_forces();
 private:
  double k_norm;
};

/* ---------------------------------------------------------------------- */

class NormalHertzMaterial:NormalHertz
{
 public:
  NormalHertzMaterial();
  ~NormalHertzMaterial() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
 private:
  double k_norm;
};

/* ---------------------------------------------------------------------- */

class NormalDMT:NormalModel
{
 public:
  NormalDMT();
  ~NormalDMT() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  void set_fncrit();
  double calculate_forces();
 private:
  double k_norm, cohesion;
  double F_pulloff;
};

/* ---------------------------------------------------------------------- */

class NormalJKR:NormalModel
{
 public:
  NormalJKR();
  ~NormalJKR() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  void set_fncrit();
  double calculate_forces();
  void pulloff_distance(double, double);
 private:
  double k_norm, cohesion;
  double Escaled, F_pulloff;
};
}

#endif /*NORMAL_CONTACT_MODELS_H_ */

