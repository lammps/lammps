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

#ifndef PAIR_COUL_DEBYE_H
#define PAIR_COUL_DEBYE_H

#include "pair_coul_cut.h"

namespace LAMMPS_NS {

class PairCoulDebye : public PairCoulCut {
 public:
  PairCoulDebye(class LAMMPS *);
  void compute(int, int);
  void settings(int, char **);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 private:
  double kappa;
};

}

#endif
