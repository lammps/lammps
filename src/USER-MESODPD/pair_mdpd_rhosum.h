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
PairStyle(mdpd/rhosum,PairMDPDRhoSum);
// clang-format on
#else

#ifndef LMP_PAIR_MDPD_RHOSUM_H
#define LMP_PAIR_MDPD_RHOSUM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMDPDRhoSum : public Pair {
 public:
  PairMDPDRhoSum(class LAMMPS *);
  virtual ~PairMDPDRhoSum();
  void init_style();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual double init_one(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

 protected:
  double **cut;
  int first;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
