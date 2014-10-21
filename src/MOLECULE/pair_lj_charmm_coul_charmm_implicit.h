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

PairStyle(lj/charmm/coul/charmm/implicit,PairLJCharmmCoulCharmmImplicit)

#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_H
#define LMP_PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_H

#include "pair_lj_charmm_coul_charmm.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulCharmmImplicit : public PairLJCharmmCoulCharmm {
 public:
  PairLJCharmmCoulCharmmImplicit(class LAMMPS *);
  virtual void compute(int, int);
  double single(int, int, int, int, double, double, double, double &);
};

}

#endif
#endif
