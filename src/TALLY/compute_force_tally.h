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
ComputeStyle(force/tally,ComputeForceTally);
// clang-format on
#else

#ifndef LMP_COMPUTE_FORCE_TALLY_H
#define LMP_COMPUTE_FORCE_TALLY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeForceTally : public Compute {

 public:
  ComputeForceTally(class LAMMPS *, int, char **);
  ~ComputeForceTally() override;

  void init() override;

  double compute_scalar() override;
  void compute_peratom() override;

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

  void pair_setup_callback(int, int) override;
  void pair_tally_callback(int, int, int, int, double, double, double, double, double,
                           double) override;

 private:
  bigint did_setup;
  int nmax, igroup2, groupbit2;
  double **fatom;
  double ftotal[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
