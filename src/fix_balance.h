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

#ifdef FIX_CLASS
// clang-format off
FixStyle(balance,FixBalance);
// clang-format on
#else

#ifndef LMP_FIX_BALANCE_H
#define LMP_FIX_BALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class LAMMPS *, int, char **);
  ~FixBalance() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup(int) override;
  void setup_pre_exchange() override;
  void pre_exchange() override;
  void pre_neighbor() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  int nevery, lbstyle, nitermax;
  double thresh, stopthresh;
  char bstr[4];
  int wtflag;               // 1 for weighted balancing
  int sortflag;             // 1 for sorting comm messages

  double imbnow;            // current imbalance factor
  double imbprev;           // imbalance factor before last rebalancing
  double imbfinal;          // imbalance factor after last rebalancing
  double maxloadperproc;    // max load on any processor
  int itercount;            // iteration count of last call to Balance
  int pending;
  bigint lastbalance;       // last timestep balancing was attempted

  class Balance *balance;
  class Irregular *irregular;

  void rebalance();
};

}    // namespace LAMMPS_NS

#endif
#endif
