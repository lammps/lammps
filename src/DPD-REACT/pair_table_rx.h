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
PairStyle(table/rx,PairTableRX);
// clang-format on
#else

#ifndef LMP_PAIR_TABLE_RX_H
#define LMP_PAIR_TABLE_RX_H

#include "pair_table.h"

namespace LAMMPS_NS {

class PairTableRX : public PairTable {
 public:
  PairTableRX(class LAMMPS *);
  ~PairTableRX() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  int nspecies;
  char *site1, *site2;
  int isite1, isite2;
  void getMixingWeights(int, double &, double &, double &, double &);
  bool fractionalWeighting;
};

}    // namespace LAMMPS_NS

#endif
#endif
