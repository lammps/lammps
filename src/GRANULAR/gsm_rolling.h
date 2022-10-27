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

#ifndef GSM_ROLLING_H_
#define GSM_ROLLING_H_

#include "gsm.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GSMRolling : public GSM {
 public:
  GSMRolling(class GranularModel *, class LAMMPS *);
  ~GSMRolling() {};
  virtual void coeffs_to_local() {};
  virtual void init() {};
  virtual void calculate_forces() = 0;
};

/* ---------------------------------------------------------------------- */

class GSMRollingNone : public GSMRolling {
 public:
  GSMRollingNone(class GranularModel *, class LAMMPS *);
  void calculate_forces() {};
};

/* ---------------------------------------------------------------------- */

class GSMRollingSDS : public GSMRolling {
 public:
  GSMRollingSDS(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void calculate_forces();
 protected:
  double k, mu, gamma;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GSM_ROLLING_H_ */
