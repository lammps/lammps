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

#ifndef CONTACT_ROLLING_MODELS_H_
#define CONTACT_ROLLING_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class RollingModel : public SubModel {
 public:
  RollingModel(class LAMMPS *);
  ~RollingModel() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual void calculate_forces() = 0;
};

/* ---------------------------------------------------------------------- */

class RollingNone : public RollingModel {
 public:
  RollingNone(class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class RollingSDS : public RollingModel {
 public:
  RollingSDS(class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double k, mu, gamma;
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_ROLLING_MODELS_H_ */
