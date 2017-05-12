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

FixStyle(nve/spin,FixNVESpin)

#else

#ifndef LMP_FIX_NVE_SPIN_H
#define LMP_FIX_NVE_SPIN_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVESpin : public FixNVE {
	friend class FixSpinDamping;
	
 public:
  FixNVESpin(class LAMMPS *, int, char **);
  virtual ~FixNVESpin() {}
  void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  int extra;
  double dts;
  double alpha_t;
  
 private:
  class FixSpinDamping *lockspindamping;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere requires atom style sphere

Self-explanatory.

E: Fix nve/sphere update dipole requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere requires extended particles

This fix can only be used for particles of a finite size.
 
E: Fix nve/sphere dlm must be used with update dipole
 
The DLM algorithm can only be used in conjunction with update dipole.


*/
