/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_IMBALANCE_GROUP_H
#define LMP_IMBALANCE_GROUP_H

#include "imbalance.h"

namespace LAMMPS_NS {

class ImbalanceGroup : public Imbalance {
 public:
  ImbalanceGroup() : Imbalance() {};
  virtual ~ImbalanceGroup() {};

  // disallow copy constructor and assignment operator
 private:
  ImbalanceGroup(const ImbalanceGroup &) {};
  ImbalanceGroup &operator=(const ImbalanceGroup &) {return *this;};

  // required member functions
 public:
  // parse options. return number of arguments consumed.
  virtual int options(LAMMPS *lmp, int narg, char **arg);
  // compute per-atom imbalance and apply to weight array
  virtual void compute(LAMMPS *lmp, double *weight);
};

}

#endif
