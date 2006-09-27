/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_H
#define PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_H

#include "pair_lj_charmm_coul_charmm.h"

class PairLJCharmmCoulCharmmImplicit : public PairLJCharmmCoulCharmm {
 public:
  void compute(int, int);
  void single(int, int, int, int, double, double, double, int, One &);
};

#endif
