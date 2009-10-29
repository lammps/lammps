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

#ifndef COMPUTE_EVENT_DISPLACE_H
#define COMPUTE_EVENT_DISPLACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEventDisplace : public Compute {
 public:
  ComputeEventDisplace(class LAMMPS *, int, char **);
  ~ComputeEventDisplace();
  void init();
  double compute_scalar();
  void reset_extra_compute_fix(char *);

 private:
  int triclinic;
  double displace_distsq;
  char *id_event;
  class Fix *fix;
};

}

#endif
