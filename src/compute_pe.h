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

#ifdef COMPUTE_CLASS

ComputeStyle(pe,ComputePE)

#else

#ifndef LMP_COMPUTE_PE_H
#define LMP_COMPUTE_PE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePE : public Compute {
 public:
  ComputePE(class LAMMPS *, int, char **);
  ~ComputePE() {}
  void init() {}
  double compute_scalar();

 private:
  int pairflag,bondflag,angleflag,dihedralflag,improperflag,kspaceflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pe must use group all

Energies computed by potentials (pair, bond, etc) are computed on all
atoms.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
