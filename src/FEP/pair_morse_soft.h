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
/* ----------------------------------------------------------------------
   Contributing author: Stefan Paquay (Eindhoven University of Technology)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(morse/soft,PairMorseSoft);
// clang-format on
#else

#ifndef LMP_PAIR_MORSE_SOFT_H
#define LMP_PAIR_MORSE_SOFT_H

#include "pair_morse.h"

namespace LAMMPS_NS {

class PairMorseSoft : public PairMorse {
 public:
  PairMorseSoft(class LAMMPS *lmp) :
      PairMorse(lmp), lambda(nullptr), nlambda(0), shift_range(1.0){};
  ~PairMorseSoft() override;
  void compute(int, int) override;

  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double **lambda;

  int nlambda;
  double shift_range;

  void allocate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
