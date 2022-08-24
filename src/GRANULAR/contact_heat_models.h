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

#ifndef CONTACT_HEAT_MODELS_H_
#define CONTACT_HEAT_MODELS_H_

#include "contact_sub_models.h"

namespace LAMMPS_NS {
namespace Contact {

class HeatModel : public SubModel {
 public:
  HeatModel(class LAMMPS *);
  ~HeatModel() {};
  virtual void coeffs_to_local() {};
  virtual void mix_coeffs(HeatModel*, HeatModel*) {};
  virtual double calculate_heat() = 0;
};

/* ---------------------------------------------------------------------- */

class HeatArea : public HeatModel {
 public:
  HeatArea(class LAMMPS *);
  void coeffs_to_local();
  void mix_coeffs(HeatModel*, HeatModel*);
  double calculate_heat();
 protected:
  double conductivity;
};

}    // namespace Contact
}    // namespace LAMMPS_NS

#endif /*CONTACT_HEAT_MODELS_H_ */
