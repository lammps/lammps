/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
This pair style is written by Daniele Rapetti (iximiel@gmail.com)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(smatb/single,PairSMATBSingle);
// clang-format on
#else

#ifndef LMP_PAIR_SMATB_SINGLE_H
#define LMP_PAIR_SMATB_SINGLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSMATBSingle : public Pair {
 public:
  PairSMATBSingle(class LAMMPS *);
  ~PairSMATBSingle() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  virtual void allocate();
  // allocated size of per-atom arrays
  int nmax;
  //allocated to store up calculation values
  double *on_eb{nullptr};
  // interaction radius, user-given
  double r0;
  // parameters user-given
  double p;
  double A;
  double q;
  double QSI;
  //cut offs, user given
  double cutOffStart;
  double cutOffEnd;
  //squared cut off end, calculated
  double cutOffEnd2;
  //polynomial for cutoff linking to zero:   Ae^p substitution
  double a3;
  double a4;
  double a5;
  //polynomial for cutoff linking to zero: QSIe^q substitution
  double x3;
  double x4;
  double x5;
};
}    // namespace LAMMPS_NS
#endif
#endif
