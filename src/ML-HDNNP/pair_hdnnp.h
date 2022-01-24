/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
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
// clang-format off
PairStyle(hdnnp,PairHDNNP);
// clang-format on
#else

#ifndef LMP_PAIR_HDNNP_H
#define LMP_PAIR_HDNNP_H

#include "pair.h"

namespace nnp {
class InterfaceLammps;
}

namespace LAMMPS_NS {

class PairHDNNP : public Pair {

 public:
  PairHDNNP(class LAMMPS *);
  virtual ~PairHDNNP();
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
  std::string emap;
  nnp::InterfaceLammps *interface;
};

}    // namespace LAMMPS_NS

#endif
#endif
