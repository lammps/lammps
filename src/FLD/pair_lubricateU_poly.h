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

PairStyle(lubricateU/poly,PairLubricateUPoly)

#else

#ifndef LMP_PAIR_LUBRICATEU_POLY_H
#define LMP_PAIR_LUBRICATEU_POLY_H

#include "pair_lubricateU.h"

namespace LAMMPS_NS {

class PairLubricateUPoly : public PairLubricateU {
 public:
  PairLubricateUPoly(class LAMMPS *);
  ~PairLubricateUPoly() {}
  void compute(int, int);
  void settings(int, char **);
  void init_style();

 private:
  void iterate(double **, int);
  void compute_RE(double **);
  void compute_RU(double **);
  void compute_Fh(double **);
};

}

#endif
#endif
