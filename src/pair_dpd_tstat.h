/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(dpd/tstat,PairDPDTstat)

#else

#ifndef LMP_PAIR_DPD_TSTAT_H
#define LMP_PAIR_DPD_TSTAT_H

#include "pair_dpd.h"

namespace LAMMPS_NS {

class PairDPDTstat : public PairDPD {
 public:
  PairDPDTstat(class LAMMPS *);
  ~PairDPDTstat() {}
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double t_start,t_stop;
};

}

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
