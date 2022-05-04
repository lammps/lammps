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
PairStyle(vashishta/table,PairVashishtaTable);
// clang-format on
#else

#ifndef LMP_PAIR_VASHISHITA_TABLE_H
#define LMP_PAIR_VASHISHITA_TABLE_H

#include "pair_vashishta.h"

namespace LAMMPS_NS {

class PairVashishtaTable : public PairVashishta {
 public:
  PairVashishtaTable(class LAMMPS *);
  ~PairVashishtaTable() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  double memory_usage() override;

 protected:
  int ntable;
  double deltaR2;
  double oneOverDeltaR2;
  double ***forceTable;        // table of forces per element pair
  double ***potentialTable;    // table of potential energies

  void twobody_table(const Param &, double, double &, int, double &);
  void setup_params() override;
  void create_tables();
};

}    // namespace LAMMPS_NS

#endif
#endif
