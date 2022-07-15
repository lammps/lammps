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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(spin/dmi,PairSpinDmi);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_DMI_H
#define LMP_PAIR_SPIN_DMI_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinDmi : public PairSpin {
 public:
  PairSpinDmi(LAMMPS *lmp) : PairSpin(lmp) {}
  ~PairSpinDmi() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  void compute_dmi(int, int, double *, double *, double *);
  void compute_dmi_mech(int, int, double, double *, double *, double *, double *);

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  double cut_spin_dmi_global;    // short range pair cutoff

 protected:
  double **DM;                                     // dmi coeff in eV
  double **v_dmx, **v_dmy, **v_dmz;                // dmi direction
  double **vmech_dmx, **vmech_dmy, **vmech_dmz;    // dmi mech direction
  double **cut_spin_dmi;                           // cutoff distance dmi

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
