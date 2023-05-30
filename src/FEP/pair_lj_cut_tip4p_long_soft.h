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
PairStyle(lj/cut/tip4p/long/soft,PairLJCutTIP4PLongSoft);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_SOFT_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_SOFT_H

#include "pair_lj_cut_coul_long_soft.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongSoft : public PairLJCutCoulLongSoft {
 public:
  PairLJCutTIP4PLongSoft(class LAMMPS *);
  ~PairLJCutTIP4PLongSoft() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart_settings(FILE *fp) override;
  void read_restart_settings(FILE *fp) override;
  void *extract(const char *, int &) override;
  double memory_usage() override;

 protected:
  int typeH, typeO;    // atom types of TIP4P water H and O atoms
  int typeA, typeB;    // angle and bond types of TIP4P water
  double alpha;        // geometric constraint parameter for TIP4P

  int nmax;            // info on off-oxygen charge sites
  int **hneigh;        // 0,1 = indices of 2 H associated with O
                       // 2 = 0 if site loc not yet computed, 1 if yes
  double **newsite;    // locations of charge sites

  void compute_newsite(double *, double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
