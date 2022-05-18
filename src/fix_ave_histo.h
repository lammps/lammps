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

#ifdef FIX_CLASS
// clang-format off
FixStyle(ave/histo,FixAveHisto);
// clang-format on
#else

#ifndef LMP_FIX_AVE_HISTO_H
#define LMP_FIX_AVE_HISTO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveHisto : public Fix {
 public:
  FixAveHisto(class LAMMPS *, int, char **);
  ~FixAveHisto() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_vector(int) override;
  double compute_array(int, int) override;

 protected:
  int me, nvalues;
  int nrepeat, nfreq, irepeat;
  bigint nvalid, nvalid_last;
  int *which, *argindex, *value2index;
  char **ids;
  FILE *fp;
  double lo, hi, binsize, bininv;
  int kind, beyond, overwrite;
  bigint filepos;

  double stats[4], stats_total[4], stats_all[4];
  double **stats_list;

  int nbins;
  double *bin, *bin_total, *bin_all;
  double **bin_list;
  double *coord;

  double *vector;
  int maxatom;

  int ave, nwindow, startstep, mode;
  char *title1, *title2, *title3;
  int iwindow, window_limit;

  void bin_one(double);
  void bin_vector(int, double *, int);
  void bin_atoms(double *, int);
  void options(int, int, char **);
  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
