/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(EVENT/PRD,FixEventPRD);
// clang-format on
#else

#ifndef LMP_FIX_EVENT_PRD_H
#define LMP_FIX_EVENT_PRD_H

#include "fix_event.h"

namespace LAMMPS_NS {

class FixEventPRD : public FixEvent {
 public:
  int event_number;         // event counter
  bigint event_timestep;    // timestep of last event on any replica
  bigint clock;             // total elapsed timesteps across all replicas
  int replica_number;       // replica where last event occurred
  int correlated_event;     // 1 if last event was correlated, 0 otherwise
  int ncoincident;          // # of simultaneous events on different replicas

  FixEventPRD(class LAMMPS *, int, char **);

  void write_restart(FILE *) override;
  void restart(char *) override;

  // methods specific to FixEventPRD, invoked by PRD

  void store_event_prd(bigint, int);

 private:
};

}    // namespace LAMMPS_NS

#endif
#endif
