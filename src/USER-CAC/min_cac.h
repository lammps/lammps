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

#ifndef LMP_CAC_MIN_H
#define LMP_CAC_MIN_H

#include "min.h"

namespace LAMMPS_NS {

class CACMin : public Min {
 public:
  CACMin(class LAMMPS *);
  virtual ~CACMin();
  virtual void init();

 protected:
  class FixCACMinimize *fix_cac_minimize;  // fix that stores auxiliary data
  
  double energy_force(int);
};

}

#endif

/* ERROR/WARNING messages:

*/
