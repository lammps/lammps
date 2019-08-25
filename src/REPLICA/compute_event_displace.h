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

ComputeStyle(event/displace,ComputeEventDisplace)

#else

#ifndef LMP_COMPUTE_EVENT_DISPLACE_H
#define LMP_COMPUTE_EVENT_DISPLACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEventDisplace : public Compute {
 public:
  ComputeEventDisplace(class LAMMPS *, int, char **);
  ~ComputeEventDisplace();
  void init();
  double compute_scalar();

  int all_events();
  void reset_extra_compute_fix(const char *);


 private:
  int triclinic;
  double displace_distsq;
  char *id_event;
  class FixEvent *fix_event;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Distance must be > 0 for compute event/displace

Self-explanatory.

E: Could not find compute event/displace fix ID

Self-explanatory.

E: Compute event/displace has invalid fix event assigned

This is an internal LAMMPS error.  Please report it to the
developers.

*/
