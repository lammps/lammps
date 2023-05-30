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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(spin/dipole/cut,PairSpinDipoleCut);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_DIPOLE_CUT_H
#define LMP_PAIR_SPIN_DIPOLE_CUT_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinDipoleCut : public PairSpin {
 public:
  double cut_coul;
  double **sigma;

  PairSpinDipoleCut(LAMMPS *);
  ~PairSpinDipoleCut() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  void compute_dipolar(int, int, double *, double *, double *, double *, double);
  void compute_dipolar_mech(int, int, double *, double *, double *, double *, double);

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  double cut_spin_long_global;    // global long cutoff distance

 protected:
  double hbar;            // reduced Planck's constant
  double mub;             // Bohr's magneton
  double mu_0;            // vacuum permeability
  double mub2mu0;         // prefactor for mech force
  double mub2mu0hbinv;    // prefactor for mag force

  double **cut_spin_long;    // cutoff distance long

  double g_ewald;
  int ewald_order;

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
