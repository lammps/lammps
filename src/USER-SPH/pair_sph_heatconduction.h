/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sph/heatconduction,PairSPHHeatConduction);
// clang-format on
#else

#ifndef LMP_PAIR_SPH_HEATCONDUCTION_H
#define LMP_PAIR_SPH_HEATCONDUCTION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSPHHeatConduction : public Pair {
 public:
  PairSPHHeatConduction(class LAMMPS *);
  virtual ~PairSPHHeatConduction();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double **cut, **alpha;
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
