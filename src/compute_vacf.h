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
ComputeStyle(vacf,ComputeVACF);
// clang-format on
#else

#ifndef LMP_COMPUTE_VACF_H
#define LMP_COMPUTE_VACF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVACF : public Compute {
 public:
  ComputeVACF(class LAMMPS *, int, char **);
  ~ComputeVACF() override;
  void init() override;
  void compute_vector() override;
  void set_arrays(int) override;

 protected:
  bigint nvacf;
  char *id_fix;
  class FixStorePeratom *fix;
};

}    // namespace LAMMPS_NS

#endif
#endif
