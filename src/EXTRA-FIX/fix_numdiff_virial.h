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
FixStyle(numdiff/virial,FixNumDiffVirial);
// clang-format on
#else

#ifndef LMP_FIX_NUMDIFF_VIRIAL_H
#define LMP_FIX_NUMDIFF_VIRIAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNumDiffVirial : public Fix {
 public:
  FixNumDiffVirial(class LAMMPS *, int, char **);
  ~FixNumDiffVirial() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double memory_usage() override;

 private:
  static constexpr int NDIR_VIRIAL = 6;    // dimension of virial and strain vectors
  double delta;                            // strain magnitude
  int maxatom;                             // allocated size of atom arrays
  int ilevel_respa;

  int pair_compute_flag;      // 0 if pair->compute is skipped
  int kspace_compute_flag;    // 0 if kspace->compute is skipped

  char *id_pe;          // name of energy compute
  class Compute *pe;    // pointer to energy compute

  double virial[NDIR_VIRIAL];     // finite diff virial components (Voigt order)
  double **temp_x;                // original coords
  double **temp_f;                // original forces
  double fixedpoint[3];           // define displacement field origin
  int dirlist[NDIR_VIRIAL][2];    // strain cartesian indices (Voigt order)

  double compute_vector(int) override;      // access function for virial
  void calculate_virial();                  // virial calculation
  void displace_atoms(int, int, double);    // apply displacement field
  void restore_atoms(int, int);             // restore original positions
  double update_energy();                   // calculate new energy
  void virial_clear();                      // set virial to zero
  void reallocate();                        // grow the atom arrays
};

}    // namespace LAMMPS_NS

#endif
#endif
