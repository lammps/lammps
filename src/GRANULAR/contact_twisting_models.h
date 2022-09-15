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

#ifndef CONTACT_TWISTING_MODELS_H_
#define CONTACT_TWISTING_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class TwistingModel : public SubModel {
 public:
  TwistingModel(class LAMMPS *);
  virtual ~TwistingModel() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual void calculate_forces() = 0;
};

/* ---------------------------------------------------------------------- */

class TwistingNone : public TwistingModel {
 public:
  TwistingNone(class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class TwistingMarshall : public TwistingModel {
 public:
  TwistingMarshall(class LAMMPS *);
  void calculate_forces();
  void init();
 protected:
  double k_tang, damp_tang, mu_tang;
};

/* ---------------------------------------------------------------------- */

class TwistingSDS : public TwistingModel {
 public:
  TwistingSDS(class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double k, mu, damp;
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_TWISTING_MODELS_H_ */
