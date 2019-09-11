/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(propel/self,FixPropelSelf)

#else

#ifndef LMP_FIX_PROPEL_SELF_H
#define LMP_FIX_PROPEL_SELF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropelSelf : public Fix {
 public:

  enum operation_modes {
    VELOCITY = 0,
    QUATERNION = 1
  };
	
  FixPropelSelf(class LAMMPS *, int, char **);
  virtual ~FixPropelSelf();
  virtual int setmask();
  virtual void post_force(int);
  // virtual void post_force_respa(int, int, int);

  double memory_usage();

private:
  double magnitude;
  int thermostat_orient;
  int mode;


  // If 0, apply fix to everything in group. If > 0, apply only to those
  // types i for which i <= n_types_filter _and_ apply_to_type[i] == 1:
  int n_types_filter;
  int *apply_to_type; //< Specifies, per type, if the fix applies to it or not.
  

  int verify_atoms_have_quaternion();

  template <int filter_by_type> void post_force_velocity(int);
  template <int filter_by_type> void post_force_quaternion(int);
	
	
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
