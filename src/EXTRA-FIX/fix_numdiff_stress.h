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
FixStyle(numdiff/stress,FixNumDiffStress);
// clang-format on
#else

#ifndef LMP_FIX_NUMDIFF_STRESS_H
#define LMP_FIX_NUMDIFF_STRESS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNumDiffStress : public Fix {
 public:
  FixNumDiffStress(class LAMMPS *, int, char **);
  ~FixNumDiffStress();
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double memory_usage() override;

 private:
  static const int NDIR_STRESS = 6; // dimension of stress and strain vectors 
  double delta;              // strain magnitude
  int maxatom;               // allocated size of atom arrays
  int ilevel_respa;

  int pair_compute_flag;      // 0 if pair->compute is skipped
  int kspace_compute_flag;    // 0 if kspace->compute is skipped

  char *id_pe;                // name of energy compute 
  class Compute *pe;          // pointer to energy compute

  double stress[NDIR_STRESS]; // finite diff stress components (Voigt order)
  double **temp_x;            // original coords
  double **temp_f;            // original forces
  double fixedpoint[3];       // define displacement field origin
  int dirlist[NDIR_STRESS][2];// strain cartesian indices (Voigt order)
  
  double compute_vector(int) override;  // access function for stress
  void calculate_stress();              // stress calculation
  void displace_atoms(int, int, double);// apply displacement field
  void restore_atoms(int, int);         // restore original positions
  double update_energy();               // calculate new energy
  void stress_clear();                  // set stress to zero
  void reallocate();                    // grow the atom arrays
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix numdiff/stress does not exist

Self-explanatory.

*/
