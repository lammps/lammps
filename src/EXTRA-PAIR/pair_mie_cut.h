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
PairStyle(mie/cut,PairMIECut);
// clang-format on
#else

#ifndef LMP_PAIR_MIE_CUT_H
#define LMP_PAIR_MIE_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMIECut : public Pair {
 public:
  PairMIECut(class LAMMPS *);
  ~PairMIECut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;

 protected:
  double cut_global;
  double **cut;
  double **epsilon, **sigma;
  double **gamR, **gamA, **Cmie;
  double **mie1, **mie2, **mie3, **mie4, **offset;
  double *cut_respa;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
