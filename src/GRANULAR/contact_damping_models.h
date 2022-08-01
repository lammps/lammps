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

#ifndef CONTACT_DAMPING_MODELS_H_
#define CONTACT_DAMPING_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class DampingModel : public SubModel {
 public:
  DampingModel() {};
  ~DampingModel() {};
  virtual void coeffs_to_local();
  virtual void mix_coeffs(DampingModel*, DampingModel*) {};
  virtual double calculate_forces() = 0;
  double damp;
};

/* ---------------------------------------------------------------------- */

class DampingVelocity: public DampingModel {
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingMassVelocity: public DampingModel {
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingViscoelastic: public DampingModel {
 public:
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class DampingTsuji: public DampingModel {
 public:
  void coeffs_to_local();
  double calculate_forces();
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_DAMPING_MODELS_H_ */
