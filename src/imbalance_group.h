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
  ImbalanceGroup(LAMMPS *lmp) : Imbalance(lmp),_num(0),_id(0),_factor(0) {};
  virtual ~ImbalanceGroup() { delete[] _id; delete[] _factor; };

  // internal data members
 private:
  int _num;                     // number of groups with weights
  int *_id;                     // list numerical id's of groups
  double *_factor;              // list if group weight factors

  // required member functions
 public:
  // parse options. return number of arguments consumed.
  virtual int options(int narg, char **arg);
  // compute and apply weight factors to local atom array
  virtual void compute(double *weight);
  // print information about the state of this imbalance compute
  virtual void info(FILE *fp);
};

}

#endif
