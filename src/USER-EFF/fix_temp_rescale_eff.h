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

#ifdef FIX_CLASS

FixStyle(temp/rescale/eff,FixTempRescaleEff)

#else

#ifndef LMP_FIX_TEMP_RESCALE_EFF_H
#define LMP_FIX_TEMP_RESCALE_EFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempRescaleEff : public Fix {
 public:
  FixTempRescaleEff(class LAMMPS *, int, char **);
  virtual ~FixTempRescaleEff();
  int setmask();
  void init();
  virtual void end_of_step();
  int modify_param(int, char **);
  void reset_target(double);
  double compute_scalar();

 protected:
  int which;
  double t_start,t_stop,t_window;
  double fraction,energy,efactor;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}

#endif
#endif
