/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(list,PairList);
// clang-format on
#else

#ifndef LMP_PAIR_LIST_H
#define LMP_PAIR_LIST_H

#include "pair.h"

namespace LAMMPS_NS {

class PairList : public Pair {
 public:
  PairList(class LAMMPS *);
  ~PairList() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

 protected:
  void allocate();

  // potential specific parameters
  struct harm_p {
    double k, r0;
  };
  struct morse_p {
    double d0, alpha, r0;
  };
  struct lj126_p {
    double epsilon, sigma;
  };
  struct quartic_p {
    double k, r0, b1, b2;
  };

  union param_u {
    harm_p harm;
    morse_p morse;
    lj126_p lj126;
    quartic_p quartic;
  };

  struct list_param {
    int style;              // potential style indicator
    tagint id1, id2;        // global atom ids
    double cutsq;           // cutoff**2 for this pair
    double offset;          // energy offset
    union param_u param;    // parameters for style
  };

 protected:
  double cut_global;     // global cutoff distance
  list_param *params;    // lisf of pair interaction parameters
  int npairs;            // # of atom pairs in global list
  int check_flag;        // 1 if checking for missing pairs
};

}    // namespace LAMMPS_NS

#endif
#endif
