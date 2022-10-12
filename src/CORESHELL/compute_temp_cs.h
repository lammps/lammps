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
ComputeStyle(temp/cs,ComputeTempCS);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_CS_H
#define LMP_COMPUTE_TEMP_CS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempCS : public Compute {
 public:
  ComputeTempCS(class LAMMPS *, int, char **);
  ~ComputeTempCS() override;
  void init() override;
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;
  double memory_usage() override;

  void remove_bias(int, double *) override;
  void remove_bias_all() override;
  void reapply_bias_all() override;
  void restore_bias(int, double *) override;
  void restore_bias_all() override;

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int groupbit_c, groupbit_s;
  int nshells;
  int firstflag;
  int maxatom;
  int cgroup, sgroup;

  double tfactor;
  double **vint;

  char *id_fix;
  class FixStorePeratom *fix;

  void dof_compute();
  void vcm_pairs();
};

}    // namespace LAMMPS_NS

#endif
#endif
