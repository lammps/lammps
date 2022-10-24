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
 LAMMPS development team: developers@lammps.org

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

  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  double safety_factor;
  double dt, t_elapsed;
};

}    // namespace LAMMPS_NS

#endif
#endif
