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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pe/mol/tally,ComputePEMolTally);
// clang-format on
#else

#ifndef LMP_COMPUTE_PE_MOL_TALLY_H
#define LMP_COMPUTE_PE_MOL_TALLY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePEMolTally : public Compute {

 public:
  ComputePEMolTally(class LAMMPS *, int, char **);
  virtual ~ComputePEMolTally();

  void init();
  void compute_vector();

  void pair_setup_callback(int, int);
  void pair_tally_callback(int, int, int, int, double, double, double, double, double, double);

 private:
  bigint did_setup;
  int igroup2, groupbit2;
  double etotal[4];
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
