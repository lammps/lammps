/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(event/displace,ComputeEventDisplace);
// clang-format on
#else

#ifndef LMP_COMPUTE_EVENT_DISPLACE_H
#define LMP_COMPUTE_EVENT_DISPLACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEventDisplace : public Compute {
 public:
  ComputeEventDisplace(class LAMMPS *, int, char **);
  ~ComputeEventDisplace() override;
  void init() override;
  double compute_scalar() override;

  int all_events();
  void reset_extra_compute_fix(const char *) override;

 private:
  int triclinic;
  double displace_distsq;
  char *id_event;
  class FixEvent *fix_event;
};

}    // namespace LAMMPS_NS

#endif
#endif
