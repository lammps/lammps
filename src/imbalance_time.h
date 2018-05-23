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

#ifndef LMP_IMBALANCE_TIME_H
#define LMP_IMBALANCE_TIME_H

#include "imbalance.h"

namespace LAMMPS_NS {

class ImbalanceTime : public Imbalance {
 public:
  ImbalanceTime(class LAMMPS *);
  virtual ~ImbalanceTime() {}

 public:
  // parse options, return number of arguments consumed
  virtual int options(int, char **);
  // reinitialize internal data
  virtual void init(int);
  // compute and apply weight factors to local atom array
  virtual void compute(double *);
  // print information about the state of this imbalance compute
  virtual void info(FILE *);

 private:
  double factor;               // weight factor for time imbalance
  double last;                 // combined wall time from last call
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Balance weight <= 0.0

UNDOCUMENTED

*/
