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

/* ----------------------------------------------------------------------
   Contributing author: Markus Hohnerbach (RWTH)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(rebo/intel,PairREBOIntel)

#else

#ifndef LMP_PAIR_REBO_INTEL_H
#define LMP_PAIR_REBO_INTEL_H

#include "pair_airebo_intel.h"

namespace LAMMPS_NS {

class PairREBOIntel : public PairAIREBOIntel {
 public:
  PairREBOIntel(class LAMMPS *);
  virtual void settings(int, char **);
};

}

#endif
#endif
