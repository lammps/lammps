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
PairStyle(spin/magelec,PairSpinMagelec);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_MAGELEC_H
#define LMP_PAIR_SPIN_MAGELEC_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinMagelec : public PairSpin {
 public:
  PairSpinMagelec(LAMMPS *lmp) : PairSpin(lmp) {}
  ~PairSpinMagelec() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  void compute_magelec(int, int, double *, double *, double *);
  void compute_magelec_mech(int, int, double *, double *, double *);

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  double cut_spin_magelec_global;    // global me cutoff

 protected:
  double **ME, **ME_mech;              // magelec coeff in eV
  double **v_mex, **v_mey, **v_mez;    // magelec direction
  double **cut_spin_magelec;           // magelec cutoff distance

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
