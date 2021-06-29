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
ComputeStyle(temp/drude,ComputeTempDrude);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_DRUDE_H
#define LMP_COMPUTE_TEMP_DRUDE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempDrude : public Compute {
 public:
  ComputeTempDrude(class LAMMPS *, int, char **);
  ~ComputeTempDrude();
  void init();
  void setup();
  void compute_vector();
  double compute_scalar();
  int modify_param(int, char **);

 private:
  int fix_dof;
  class FixDrude *fix_drude;
  char *id_temp;
  class Compute *temperature;
  bigint dof_core, dof_drude;
  double kineng_core, kineng_drude;
  double temp_core, temp_drude;

  void dof_compute();
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
