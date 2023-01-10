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
PairStyle(quip,PairQUIP);
// clang-format on
#else

#ifndef LMP_PAIR_QUIP_H
#define LMP_PAIR_QUIP_H

#include "pair.h"

extern "C" {
int quip_lammps_api_version();
void quip_lammps_wrapper(int *, int *, int *, int *, int *, int *, int *, int *, int *, double *,
                         int *, int *, double *, double *, double *, double *, double *, double *);
void quip_lammps_potential_initialise(int *, int *, double *, char *, int *, char *, int *);
}

namespace LAMMPS_NS {

class PairQUIP : public Pair {
 public:
  PairQUIP(class LAMMPS *);
  ~PairQUIP() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void allocate();

 private:
  double cutoff;
  int *quip_potential;
  int n_quip_potential;
  int *map;           // mapping from atom types to elements
  char *quip_file;    // mapping from atom types to elements
  int n_quip_file;
  char *quip_string;    // mapping from atom types to elements
  int n_quip_string;
};

}    // namespace LAMMPS_NS

#endif
#endif
