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
PairStyle(srp,PairSRP);
// clang-format on
#else

#ifndef LMP_PAIR_SRP_H
#define LMP_PAIR_SRP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSRP : public Pair {
 public:
  PairSRP(class LAMMPS *);
  ~PairSRP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

 protected:
  inline void onetwoexclude(int *&, int &, int *&, int *&, int **&);
  inline void remapBonds(int &);
  void allocate();
  void getMinDist(double **&, double &, double &, double &, double &, double &, int &, int &, int &,
                  int &);
  bool min, midpoint;
  double **cut;
  double **a0;
  double **srp;
  double cut_global;
  int bptype;
  int btype;
  class Fix *f_srp;
  char *fix_id;
  int exclude, maxcount;
  int **segment;
};

}    // namespace LAMMPS_NS

#endif
#endif
