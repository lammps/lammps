/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
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
FixStyle(smd/adjust_dt,FixSMDTlsphDtReset);
// clang-format on
#else

#ifndef LMP_FIX_TLSPH_DT_RESET_H
#define LMP_FIX_TLSPH_DT_RESET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMDTlsphDtReset : public Fix {
 public:
  FixSMDTlsphDtReset(class LAMMPS *, int, char **);
  ~FixSMDTlsphDtReset() {}
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void end_of_step();
  double compute_scalar();
  void write_restart(FILE *);
  void restart(char *);

 private:
  double safety_factor;
  double dt, t_elapsed;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

 E: Illegal ... command

 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.

 E: Use of fix dt/reset with undefined lattice

 Must use lattice command with fix dt/reset command if units option is
 set to lattice.

 W: Dump dcd/xtc timestamp may be wrong with fix dt/reset

 If the fix changes the timestep, the dump dcd file will not
 reflect the change.

 */
