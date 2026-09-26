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
PairStyle(bpm/peri,PairBPMPeri);
// clang-format on
#else

#ifndef LMP_PAIR_BPM_PERI_H
#define LMP_PAIR_BPM_PERI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBPMPeri : public Pair {
 public:
  PairBPMPeri(class LAMMPS *);
  ~PairBPMPeri() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double **kspring;    // short-range contact stiffness per type pair
  double **cut;        // horizon (also the neighbor cutoff) per type pair
  int index_vfrac;     // per-atom d_vfrac (nodal volume) index

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
