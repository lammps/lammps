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

FixStyle(rigid/meso,FixRigidMeso)

#else

#ifndef LMP_FIX_RIGID_MESO_H
#define LMP_FIX_RIGID_MESO_H

#include "fix_rigid.h"

namespace LAMMPS_NS {

class FixRigidMeso : public FixRigid {
 public:
  FixRigidMeso (class LAMMPS *, int, char **);
  ~FixRigidMeso ();
  int setmask ();
  void setup (int);
  void initial_integrate (int);
  void final_integrate ();
  double compute_scalar () { return 0.0; }
  double compute_array (int, int);

 protected:
  void set_xv ();
  void set_v ();
  double **conjqm;                    // conjugate quaternion momentum
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: fix rigid/meso command requires atom_style with both energy and density

You should use atom_style meso with this fix

E: Can not use thermostat with fix rigid/meso

Self-explanatory

E: Can not use barostat with fix rigid/meso

Self-explanatory

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
