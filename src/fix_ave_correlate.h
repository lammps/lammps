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
FixStyle(ave/correlate,FixAveCorrelate);
// clang-format on
#else

#ifndef LMP_FIX_AVE_CORRELATE_H
#define LMP_FIX_AVE_CORRELATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveCorrelate : public Fix {
 public:
  FixAveCorrelate(class LAMMPS *, int, char **);
  ~FixAveCorrelate() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_array(int, int) override;

 private:
  int me, nvalues;
  int nrepeat, nfreq;
  bigint nvalid, nvalid_last;
  int *which, *argindex, *value2index;
  char **ids;
  FILE *fp;

  int type, ave, startstep, overwrite;
  double prefactor;
  bigint filepos;

  int firstindex;    // index in values ring of earliest time sample
  int lastindex;     // index in values ring of latest time sample
  int nsample;       // number of time samples in values ring

  int npair;    // number of correlation pairs to calculate
  int *count;
  double **values, **corr;

  int *save_count;    // saved values at Nfreq for output via compute_array()
  double **save_corr;

  void accumulate();
  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
