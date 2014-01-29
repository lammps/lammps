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

ComputeStyle(vacf,ComputeVACF)

#else

#ifndef LMP_COMPUTE_VACF_H
#define LMP_COMPUTE_VACF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVACF : public Compute {
 public:
  ComputeVACF(class LAMMPS *, int, char **);
  ~ComputeVACF();
  void init();
  virtual void compute_vector();

 protected:
  bigint nvacf;
  char *id_fix;
  class FixStore *fix;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute vacf fix ID

UNDOCUMENTED

U: Could not find compute msd fix ID

Self-explanatory.

*/
