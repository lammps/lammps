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
  ImbalanceNeigh(LAMMPS *lmp) : Imbalance(lmp), _factor(0.0) {};
  virtual ~ImbalanceNeigh() {};

  // internal data members
 private:
  double _factor;               // weight factor for neighbor imbalance

 public:
  // parse options. return number of arguments consumed
  virtual int options(int narg, char **arg);
  // compute and apply weight factors to local atom array
  virtual void compute(double *weights);
  // print information about the state of this imbalance compute
  virtual void info(FILE *fp);
};

}

#endif
