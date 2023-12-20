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
PairStyle(spin/exchange,PairSpinExchange);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_EXCHANGE_H
#define LMP_PAIR_SPIN_EXCHANGE_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinExchange : public PairSpin {
 public:
  PairSpinExchange(class LAMMPS *);
  ~PairSpinExchange() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  void compute_exchange(int, int, double, double *, double *);
  void compute_exchange_mech(int, int, double, double *, double *, double *, double *);
  double compute_energy(int, int, double, double *, double *);

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  double cut_spin_exchange_global;    // global exchange cutoff distance

 protected:
  int e_offset;                  // apply energy offset
  double **J1_mag;               // exchange coeffs in eV
  double **J1_mech;              // mech exchange coeffs in
  double **J2, **J3;             // J1 in eV, J2 adim, J3 in Ang
  double **cut_spin_exchange;    // cutoff distance exchange

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
