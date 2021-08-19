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
PairStyle(oxrna2/hbond,PairOxrna2Hbond);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_HBOND_H
#define LMP_PAIR_OXRNA2_HBOND_H

#include "pair_oxdna_hbond.h"

namespace LAMMPS_NS {

class PairOxrna2Hbond : public PairOxdnaHbond {
 public:
  PairOxrna2Hbond(class LAMMPS *);
  virtual ~PairOxrna2Hbond() {}
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
