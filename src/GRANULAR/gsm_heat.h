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

#ifdef GSM_CLASS
// clang-format off
GSMStyle(none,
         GSMHeatNone,
         HEAT);

GSMStyle(area,
         GSMHeatArea,
         HEAT);
// clang-format on
#else

#ifndef GSM_HEAT_H_
#define GSM_HEAT_H_

#include "gsm.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMHeat : public GSM {
 public:
  GSMHeat(class GranularModel *, class LAMMPS *);
  ~GSMHeat() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual double calculate_heat() = 0;
};

/* ---------------------------------------------------------------------- */

class GSMHeatNone : public GSMHeat {
 public:
  GSMHeatNone(class GranularModel *, class LAMMPS *);
  double calculate_heat();
};

/* ---------------------------------------------------------------------- */

class GSMHeatArea : public GSMHeat {
 public:
  GSMHeatArea(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  double calculate_heat();
 protected:
  double conductivity;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_HEAT_H_ */
#endif /*GSM_CLASS_H_ */
