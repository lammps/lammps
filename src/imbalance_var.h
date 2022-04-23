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

#ifndef LMP_IMBALANCE_VAR_H
#define LMP_IMBALANCE_VAR_H

#include "imbalance.h"

namespace LAMMPS_NS {

class ImbalanceVar : public Imbalance {
 public:
  ImbalanceVar(class LAMMPS *);
  ~ImbalanceVar() override;

 public:
  // parse options. return number of arguments consumed.
  int options(int, char **) override;
  // re-initialize internal data, e.g. variable ID
  void init(int) override;
  // compute per-atom imbalance and apply to weight array
  void compute(double *) override;
  // print information about the state of this imbalance compute (required)
  std::string info() override;

 private:
  char *name;    // variable name
  int id;        // variable index
};

}    // namespace LAMMPS_NS

#endif
