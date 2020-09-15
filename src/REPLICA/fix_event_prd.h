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

FixStyle(EVENT/PRD,FixEventPRD)

#else

#ifndef LMP_FIX_EVENT_PRD_H
#define LMP_FIX_EVENT_PRD_H

#include "fix_event.h"

namespace LAMMPS_NS {

class FixEventPRD : public FixEvent {
 public:
  int event_number;      // event counter
  bigint event_timestep; // timestep of last event on any replica
  bigint clock;          // total elapsed timesteps across all replicas
  int replica_number;    // replica where last event occurred
  int correlated_event;  // 1 if last event was correlated, 0 otherwise
  int ncoincident;       // # of simultaneous events on different replicas

  FixEventPRD(class LAMMPS *, int, char **);
  ~FixEventPRD() {}

  void write_restart(FILE *);
  void restart(char *);

  // methods specific to FixEventPRD, invoked by PRD

  void store_event_prd(bigint, int);

 private:

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
