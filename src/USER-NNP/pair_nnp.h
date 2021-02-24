/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   This file initially came from n2p2 (https://github.com/CompPhysVienna/n2p2)
   Copyright (2018) Andreas Singraber (University of Vienna)

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andreas Singraber
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(nnp,PairNNP)

#else

#ifndef LMP_PAIR_NNP_H
#define LMP_PAIR_NNP_H

#include "pair.h"

namespace nnp {
    class InterfaceLammps;
}

namespace LAMMPS_NS {

class PairNNP : public Pair {

 public:

  PairNNP(class LAMMPS *);
  virtual ~PairNNP();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);

 protected:

  virtual void allocate();
  void transferNeighborList();
  void handleExtrapolationWarnings();

  bool showew;
  bool resetew;
  int showewsum;
  int maxew;
  long numExtrapolationWarningsTotal;
  long numExtrapolationWarningsSummary;
  double cflength;
  double cfenergy;
  double maxCutoffRadius;
  char *directory;
  char *emap;
  nnp::InterfaceLammps *interface;
};

}

#endif
#endif
