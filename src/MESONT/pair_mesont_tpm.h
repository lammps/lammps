/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mesont/tpm,PairMESONTTPM);
// clang-format on
#else

#ifndef LMP_PAIR_MESONT_TPM_H
#define LMP_PAIR_MESONT_TPM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMESONTTPM : public Pair {
 public:
  PairMESONTTPM(class LAMMPS *);
  ~PairMESONTTPM() override;
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
  void init_style() override;

  double energy_s;                        // accumulated energies for stretching
  double energy_b;                        // accumulated energies for bending
  double energy_t;                        // accumulated energies for tube-tube interaction
  double *eatom_s, *eatom_b, *eatom_t;    // accumulated per-atom values

 protected:
  int BendingMode, TPMType;
  char *tab_path;
  int tab_path_length;
  double cut_global;
  double **cut;
  static int instance_count;
  int nmax;

  virtual void allocate();
  void *extract(const char *, int &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
