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

#ifdef FIX_CLASS
// clang-format off
FixStyle(numdiff,FixNumDiff);
// clang-format on
#else

#ifndef LMP_FIX_NUMDIFF_H
#define LMP_FIX_NUMDIFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNumDiff : public Fix {
 public:
  FixNumDiff(class LAMMPS *, int, char **);
  ~FixNumDiff() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double memory_usage() override;

 private:
  double delta;
  int maxatom;
  int ilevel_respa;

  int pair_compute_flag;      // 0 if pair->compute is skipped
  int kspace_compute_flag;    // 0 if kspace->compute is skipped

  char *id_pe;
  class Compute *pe;

  double **numdiff_forces;    // finite diff forces
  double **temp_x;            // original coords
  double **temp_f;            // original forces

  void calculate_forces();
  void displace_atoms(int, int, int);
  void restore_atoms(int, int);
  double update_energy();
  void force_clear(double **);
  void reallocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
