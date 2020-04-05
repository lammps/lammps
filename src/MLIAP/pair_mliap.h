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

PairStyle(mliap,PairMLIAP)

#else

#ifndef LMP_PAIR_MLIAP_H
#define LMP_PAIR_MLIAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMLIAP : public Pair {
public:
  PairMLIAP(class LAMMPS *);
  ~PairMLIAP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();
  int *map;                     // mapping from atom types to elements

protected:
  virtual void allocate();

  double** beta;                // betas for all atoms in list
  double** descriptors;         // descriptors for all atoms in list
  int ndescriptors;             // number of descriptors 
  int beta_max;                 // number of atoms allocated for beta, descriptors

  class MLIAPModel* model;
  class MLIAPDescriptor* descriptor;
};

}

#endif
#endif

