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

ComputeStyle(heat/flux,ComputeHeatFlux)

#else

#ifndef LMP_COMPUTE_HEAT_FLUX_H
#define LMP_COMPUTE_HEAT_FLUX_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHeatFlux : public Compute {
 public:
  ComputeHeatFlux(class LAMMPS *, int, char **);
  ~ComputeHeatFlux();
  void init();
  void compute_vector();

 private:
  char *id_ke,*id_pe,*id_stress;
  class Compute *c_ke,*c_pe,*c_stress;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute heat/flux compute ID

Self-explanatory.

E: Compute heat/flux compute ID does not compute ke/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute pe/atom

Self-explanatory.

E: Compute heat/flux compute ID does not compute stress/atom

Self-explanatory.

*/
