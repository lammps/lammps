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
  ImbalanceGroup(class LAMMPS *);
  virtual ~ImbalanceGroup();

  // parse options, return number of arguments consumed
  virtual int options(int, char **);
  // compute and apply weight factors to local atom array
  virtual void compute(double *);
  // print information about the state of this imbalance compute
  virtual void info(FILE *);

 private:
  int num;                     // number of groups with weights
  int *id;                     // numerical IDs of groups
  double *factor;              // group weight factors
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Unknown group in balance weight command

UNDOCUMENTED

*/
