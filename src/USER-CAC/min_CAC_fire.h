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
#ifdef MINIMIZE_CLASS

MinimizeStyle(cac/fire,CACMinFire)

#else

#ifndef LMP_CAC_MIN_FIRE_H
#define LMP_CAC_MIN_FIRE_H

#include "min.h"

namespace LAMMPS_NS {

class CACMinFire : public Min {
 public:
  CACMinFire(class LAMMPS *);
  ~CACMinFire() {}
  void init();
  void setup_style();
  void reset_vectors();
  virtual int iterate(int);
 protected:
  virtual void copy_vectors();
 private:
  double dt,dtmax;
  double alpha;
  bigint last_negative;
};

}

#endif
#endif
