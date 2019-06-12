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

#ifdef COMPUTE_CLASS

ComputeStyle(CAC/nodal_temp,ComputeNodalTemp)

#else

#ifndef LMP_COMPUTE_NODAL_TEMP_H
#define LMP_COMPUTE_NODAL_TEMP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNodalTemp : public Compute {
 public:
  ComputeNodalTemp(class LAMMPS *, int, char **);
  virtual ~ComputeNodalTemp();
  void init() {}
  void setup();
  virtual double compute_scalar();

 protected:
  double tfactor;

  virtual void dof_compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
