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
PairStyle(body/nparticle,PairBodyNparticle);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_NPARTICLE_H
#define LMP_PAIR_BODY_NPARTICLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBodyNparticle : public Pair {
 public:
  PairBodyNparticle(class LAMMPS *);
  ~PairBodyNparticle() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4;

  class AtomVecBody *avec;
  class BodyNparticle *bptr;

  double **discrete;    // list of all sub-particles for all bodies
  int ndiscrete;        // number of discretes in list
  int dmax;             // allocated size of discrete list
  int *dnum;            // number of discretes per line, 0 if uninit
  int *dfirst;          // index of first discrete per each line
  int nmax;             // allocated size of dnum,dfirst vectors

  void allocate();
  void body2space(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
