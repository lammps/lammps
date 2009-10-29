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

#ifndef FIX_EVENT_H
#define FIX_EVENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvent : public Fix {
 public:
  int event_number;      // event counter
  int event_timestep;    // timestep of last event on any replica
  int clock;             // total elapsed timesteps across all replicas
  int replica_number;    // replica where last event occured
  int correlated_event;  // 1 if last event was correlated, 0 otherwise

  FixEvent(class LAMMPS *, int, char **);
  ~FixEvent();
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  void write_restart(FILE *);
  void restart(char *);

  // methods specific to FixEvent, invoked by PRD

  void store_event(int, int);
  void store_state();
  void restore_state();

 private:
  double **xevent;       // atom coords at last event
  double **xold;         // atom coords for reset/restore
  int *imageold;         // image flags for reset/restore
};

}

#endif
