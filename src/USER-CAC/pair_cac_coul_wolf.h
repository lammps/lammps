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

#ifdef PAIR_CLASS

PairStyle(cac/coul/wolf,PairCACCoulWolf)

#else

#ifndef LMP_PAIR_COUL_WOLF_CAC_H
#define LMP_PAIR_COUL_WOLF_CAC_H

#include "pair_cac.h"

namespace LAMMPS_NS {

class PairCACCoulWolf : public PairCAC {
 public:
  PairCACCoulWolf(class LAMMPS *);
  virtual ~PairCACCoulWolf();
  
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);

 protected:

  double cut_coul, cut_coulsq, alf;
  
  void allocate();
  void force_densities(int, double, double, double, double, double
    &fx, double &fy, double &fz);
  virtual void settings(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unexpected argument in cac/coul/wolf invocation

Self-explanatory. Check the input script. See the documention
for appropriate syntax.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cac/coul/wolf requires atom attribute q

The atom style defined does not have this attribute.

*/