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

#ifndef CONTACT_TANGENTIAL_MODELS_H_
#define CONTACT_TANGENTIAL_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class TangentialModel : public SubModel {
 public:
  TangentialModel(class LAMMPS *);
  virtual ~TangentialModel() {};
  virtual void coeffs_to_local() {};
  virtual void mix_coeffs(TangentialModel*, TangentialModel*) {};
  virtual void calculate_forces() = 0;
  void init();
  int rescale_flag;
  double k, damp, mu; // Used by Marshall twisting model
 protected:
  double xt;
};

/* ---------------------------------------------------------------------- */

class TangentialLinearNoHistory : public TangentialModel {
 public:
  TangentialLinearNoHistory(class LAMMPS *);
  void coeffs_to_local();
  void mix_coeffs(TangentialModel*, TangentialModel*);
  void calculate_forces();
};

/* ---------------------------------------------------------------------- */

class TangentialLinearHistory : public TangentialModel {
 public:
  TangentialLinearHistory(class LAMMPS *);
  void coeffs_to_local();
  void mix_coeffs(TangentialModel*, TangentialModel*);
  void calculate_forces();
};

/* ---------------------------------------------------------------------- */

class TangentialMindlin : public TangentialModel {
 public:
  TangentialMindlin(class LAMMPS *);
  void coeffs_to_local();
  void mix_coeffs(TangentialModel*, TangentialModel*);
  void calculate_forces();
 protected:
  int mindlin_rescale, mindlin_force;
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinForce : public TangentialMindlin {
 public:
  TangentialMindlinForce(class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinRescale : public TangentialMindlin {
 public:
  TangentialMindlinRescale(class LAMMPS *);
};

/* ---------------------------------------------------------------------- */

class TangentialMindlinRescaleForce : public TangentialMindlin {
 public:
  TangentialMindlinRescaleForce(class LAMMPS *);
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_TANGENTIAL_MODELS_H_ */
