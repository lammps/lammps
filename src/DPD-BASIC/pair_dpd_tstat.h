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
PairStyle(dpd/tstat,PairDPDTstat);
// clang-format on
#else

#ifndef LMP_PAIR_DPD_TSTAT_H
#define LMP_PAIR_DPD_TSTAT_H

#include "pair_dpd.h"

namespace LAMMPS_NS {

class PairDPDTstat : public PairDPD {
 public:
  PairDPDTstat(class LAMMPS *);
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;

 protected:
  double t_start, t_stop;
};

}    // namespace LAMMPS_NS

#endif
#endif
