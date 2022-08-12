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

#ifndef CONTACT_NORMAL_MODELS_H_
#define CONTACT_NORMAL_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class NormalModel : public SubModel {
 public:
  NormalModel();
  ~NormalModel() {};
  virtual void coeffs_to_local() {};
  virtual void mix_coeffs(NormalModel*, NormalModel*) {};
  virtual bool touch();
  virtual double pulloff_distance(double, double);
  virtual double calculate_area();
  virtual double calculate_forces() = 0;
  virtual void set_knfac() = 0;
  virtual void set_fncrit();
  double damp;  // Vestigial argument needed by damping
  double Emod, poiss;
  double Fncrit, Fne, knfac;
  int material_properties;
};

/* ---------------------------------------------------------------------- */

class NormalHooke: public NormalModel {
 public:
  NormalHooke();
  ~NormalHooke() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  double calculate_forces();
  void set_knfac();
 private:
  double k;
};

/* ---------------------------------------------------------------------- */

class NormalHertz: public NormalModel {
 public:
  NormalHertz();
  ~NormalHertz() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  double calculate_forces();
  void set_knfac();
 private:
  double k;
};

/* ---------------------------------------------------------------------- */

class NormalHertzMaterial: public NormalHertz {
 public:
  NormalHertzMaterial();
  ~NormalHertzMaterial() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
 private:
  double k;
};

/* ---------------------------------------------------------------------- */

class NormalDMT: public NormalModel {
 public:
  NormalDMT();
  ~NormalDMT() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  double calculate_forces();
  void set_knfac();
  void set_fncrit();
 private:
  double k, cohesion;
  double F_pulloff;
};

/* ---------------------------------------------------------------------- */

class NormalJKR: public NormalModel {
 public:
  NormalJKR();
  ~NormalJKR() {};
  void coeffs_to_local();
  void mix_coeffs(NormalModel*, NormalModel*);
  bool touch();
  double pulloff_distance(double, double);
  double calculate_area();
  double calculate_forces();
  void set_knfac();
  void set_fncrit();
 private:
  double cohesion;
  double Escaled, F_pulloff;
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_NORMAL_MODELS_H_ */
