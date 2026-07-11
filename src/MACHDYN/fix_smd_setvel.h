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

/* -*- c++ -*- ----------------------------------------------------------
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
FixStyle(smd/setvel,FixSMDSetVel);
// clang-format on
#else

#ifndef LMP_FIX_SMD_SETVEL_H
#define LMP_FIX_SMD_SETVEL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMDSetVel : public Fix {
 public:
  FixSMDSetVel(class LAMMPS *, int, char **);
  ~FixSMDSetVel() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  //void initial_integrate(int);
  void post_force(int) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  char *idregion;
  class Region *region;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double foriginal[3], foriginal_all[3];
  int force_flag;

  int maxatom;
  double **sforce;
};

}    // namespace LAMMPS_NS

#endif
#endif
