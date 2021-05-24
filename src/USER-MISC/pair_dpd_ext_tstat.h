/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dpd/ext/tstat,PairDPDExtTstat);
// clang-format on
#else

#ifndef LMP_PAIR_DPD_EXT_TSTAT_H
#define LMP_PAIR_DPD_EXT_TSTAT_H

#include "pair_dpd_ext.h"

namespace LAMMPS_NS {

class PairDPDExtTstat : public PairDPDExt {
 public:
  PairDPDExtTstat(class LAMMPS *);
  ~PairDPDExtTstat() {}
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);

 protected:
  double t_start, t_stop;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
