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

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdi/aimd,FixMDIAimd);
// clang-format on
#else

#ifndef LMP_FIX_MDI_AIMD_H
#define LMP_FIX_MDI_AIMD_H

#include "fix.h"
#include <mdi.h>

namespace LAMMPS_NS {

class FixMDIAimd : public Fix {
 public:
  // MDI communicator, public so that LAMMPS can work with a plugin

  MDI_Comm mdicomm;

  FixMDIAimd(class LAMMPS *, int, char **);
  ~FixMDIAimd();
  int setmask();

  void setup(int);
  void setup_pre_reverse(int, int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int nprocs;
  int plugin;

  int eflag_caller;
  double engine_energy;
  int lmpunits;

  // unit conversion factors

  double lmp2mdi_length, mdi2lmp_length;
  double lmp2mdi_energy, mdi2lmp_energy;
  double lmp2mdi_force, mdi2lmp_force;
  double lmp2mdi_pressure, mdi2lmp_pressure;
  double lmp2mdi_velocity, mdi2lmp_velocity;

  // buffers for MDI comm

  int maxbuf;
  double *buf3, *buf3all;

  // methods

  void reallocate();
  void unit_conversions();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
