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

#ifndef LMP_IMBALANCE_NEIGH_H
#define LMP_IMBALANCE_NEIGH_H

#include "imbalance.h"

namespace LAMMPS_NS {

class ImbalanceNeigh : public Imbalance {
 public:
  ImbalanceNeigh(class LAMMPS *);
  virtual ~ImbalanceNeigh() {}

 public:
  // parse options, return number of arguments consumed
  virtual int options(int, char **);
  // compute and apply weight factors to local atom array
  virtual void compute(double *);
  // print information about the state of this imbalance compute
  virtual void info(FILE *);

 private:
  double factor;               // weight factor for neighbor imbalance
  int did_warn;                // 1 if warned about no suitable neighbor list
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

W: Balance weight neigh skipped b/c no list found

UNDOCUMENTED

E: Balance weight <= 0.0

UNDOCUMENTED

*/
