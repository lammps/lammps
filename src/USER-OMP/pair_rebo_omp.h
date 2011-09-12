/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(rebo/omp,PairREBOOMP)

#else

#ifndef LMP_PAIR_REBO_OMP_H
#define LMP_PAIR_REBO_OMP_H

#include "pair_airebo_omp.h"

namespace LAMMPS_NS {

class PairREBOOMP : public PairAIREBOOMP {
 public:
  PairREBOOMP(class LAMMPS *);
  virtual void settings(int, char **);
};

}

#endif
#endif
